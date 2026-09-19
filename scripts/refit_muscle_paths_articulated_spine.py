"""Give the articulated-spine variant's muscles moment arms about its new joints.

    .venv/bin/python -m scripts.refit_muscle_paths_articulated_spine \
        [--work DIR] [--threads N] [--reuse] [--install]

WHY THIS EXISTS, and it is not the reason the variant's own document first gave.

`docs/ARTICULATED_SPINE.md` measured that no muscle has a moment arm about any of
the fifteen coordinates `data/models/articulated_spine_v1` adds -- exactly zero,
through the engine -- and attributed it to the muscle paths being fitted
polynomials in the coordinates that existed when they were fitted.  That is true,
and it is only half of a diagnosis, because it does not separate two cases that
need opposite fixes:

  REPRESENTATION  the muscle's real geometry DOES cross the joint, and the fitted
                  polynomial cannot express the arm because the joint's
                  coordinate is not one of its arguments.  A refit recovers it.
  TOPOLOGY        no muscle in this model has a path point on either side of the
                  joint.  The arm is zero in the geometry itself, the refit will
                  agree, and only a NEW MUSCLE can change it.

This script measures which each coordinate is, by sampling the model's OWN
`GeometryPath`s -- wrap objects included -- through
`scripts/native_polynomial_path_fit sample`, and comparing that ground truth to
the shipped `FunctionBasedPathSet` evaluated on the same model at the same poses.
Then it refits, and gates the refit against the same ground truth.

WHAT THE DIAGNOSIS FOUND (2026-09-18, reproduced by `--reuse` from the report):

  subtalar_angle_{l,r}   REPRESENTATION.  11 muscles per side carry a real arm,
                         up to 32.7 mm (perbrev).  The fitted set reports 0.
  mtp_angle_{l,r}        REPRESENTATION, AND IT IS NOT NEW.  4 muscles per side,
                         up to 25.3 mm (ehl).  These coordinates are in the BASE
                         plant, free and unlocked, and have read zero since the
                         engine was built -- because the shipped path set was
                         fitted upstream on a body with the toes WELDED
                         (`exampleMocoInverse.cpp:51`,
                         `ModOpReplaceJointsWithWelds({"mtp_r","mtp_l"})`, applied
                         BEFORE `ModOpReplacePathsWithFunctionBasedPaths`).
                         `scripts/native_mechanical_stream.cpp` does not weld
                         them.  So the base plant has carried two free toe hinges
                         no muscle can control since before this variant existed.
  the nine spine/neck    TOPOLOGY.  The GeometryPath truth is zero too.  This
  and the four wrist     model holds no muscle with a path point above the torso
  coordinates            or beyond the radius.  No fit can help.

HOW THE FIX WORKS.  `PolynomialPathFitter` does not take a coordinate list: it
calls `AbstractGeometryPath::findIndependentCoordinates` per path and discovers
the coordinates from the geometry (`PolynomialPathFitter.cpp:1151`).  Refitting
on a model whose subtalar is a PinJoint therefore picks subtalar up by itself.
This was already demonstrated before this script existed and is on disk:
`scripts/verify_path_refitting.py`'s gated refit of the UNMODIFIED source model
recovered `mtp_angle` for exactly edl/ehl/fdl/fhl on both sides
(`data/derived/path-refitting/verify/base-a/`), which the shipped set omits.

WHAT IS AND IS NOT REFITTED.  The refit covers all 98 `PathActuator`s, but the
installed set is FILTERED to the 80 source muscles -- the same 80 the shipped set
covers.  The 18 added muscles (`arm26_*`, `gait2392_*`) keep their real
`GeometryPath`s, which are exact, and this variant does not change them.  That
keeps the change to exactly what the diagnosis says needs changing.

OUTCOME, 2026-09-18 -- read this before the rest.  Two REFITS were pre-registered
and both FAILED the same no-regression gate (knee_angle_r: gasmed_r 1.66 mm, then
gaslat_r 4.50 mm, against a 1.36 mm bar).  The fitter's own log showed why the
second was worse: two of its four fits discarded 706 and 1,170 of 1,485 sample
rows as NaN, and even the clean fits put the right gastrocnemius knee arm at
either ~0.2 or ~1.65 mm.  The refit route is closed here.  What is INSTALLED
(attempt 3, pre-registered after attempt 2 failed) is not a fit at all: the
shipped set with the 22 subtalar/mtp-crossing muscles removed, so the engine
runs those 22 on the model's own GeometryPaths.  Every verdict and both
pre-registrations are in the report this script writes
(data/models/articulated_spine_v1/muscle_paths_report.json).  The name of this
script is kept because the diagnosis and the two failed refits are what it is
for.

THE FITTER IS NOT DETERMINISTIC (`LatinHypercubeDesign` seeds `std::mt19937` from
`std::random_device`), so the installed XML is a committed artefact, not something
`build_articulated_spine.py` regenerates.  The builder stays a function of
committed inputs and only hashes this file.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import math
import shutil
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ihm.native import path_refitting as pf  # noqa: E402

RAW = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
PATHSET = 'subject_walk_scaled_FunctionBasedPathSet.xml'
MODEL = ROOT / 'data/models/articulated_spine_v1/model.osim'
BASE_MODEL = ROOT / 'data/models/engineering_stance_v1/model.osim'
SHIPPED = ROOT / RAW / PATHSET
COORDINATES = ROOT / RAW / 'coordinates.sto'
INSTALLED = ROOT / 'data/models/articulated_spine_v1/muscle_paths.xml'
REPORT = ROOT / 'data/models/articulated_spine_v1/muscle_paths_report.json'

#: The fifteen coordinates the variant adds, plus the two the BASE model already
#: had and the fitted set already could not see.  `mtp` is in this list because
#: the diagnosis below is about it too, and leaving it out would have hidden the
#: finding that the defect predates the variant.
NEW = ('thoracic_extension', 'thoracic_bending', 'thoracic_rotation',
       'neck_extension', 'neck_bending', 'neck_rotation',
       'head_extension', 'head_bending', 'head_rotation',
       'wrist_flex_r', 'wrist_dev_r', 'wrist_flex_l', 'wrist_dev_l',
       'subtalar_angle_r', 'subtalar_angle_l')
INHERITED = ('mtp_angle_r', 'mtp_angle_l')
#: Controls: coordinates the shipped fit DOES carry.  A run in which these read
#: zero is a broken instrument, not a finding.
CONTROLS = ('ankle_angle_r', 'ankle_angle_l', 'knee_angle_r', 'lumbar_extension')

#: A moment arm below this is "no arm".  It is not a measurement of anything in
#: the body: it is two orders of magnitude below the smallest real arm the
#: diagnosis finds (ehl_r about subtalar, 4.2 mm) and four above the sampler's
#: own assembly noise, so every muscle falls unambiguously on one side of it.
ARM_FLOOR_M = 1e-4

#: Sampling the fit explores for the two foot coordinates the reference
#: trajectory does not carry.  (0.5 rad, 10 rad/s, 0 phase) is the UPSTREAM
#: example's own treatment of `mtp_angle` and is reused unchanged for
#: `subtalar_angle`.  IT IS NOT A MEASURED TRAJECTORY and is not claimed to be
#: one: what actually sets the fit's coverage is the +/-30 degree global sampling
#: bound around each frame, and the achieved subtalar span is measured back out
#: of the fitter's own sampled-values file and reported as
#: `subtalar_sampled_span_rad`.
FOOT_SINE = (0.5, 10.0, 0.0)
#: The held-out evaluation drives the same two coordinates with a DIFFERENT sine,
#: so the evaluation poses are not poses the fit centred a sample cloud on.  An
#: engineering choice, declared here rather than buried.
EVAL_SINE = (0.45, 7.0, 1.1)
#: Rows 5, 55, 105, ... -- `verify_path_refitting.py`'s own held-out selection.
EVAL_STRIDE, EVAL_OFFSET = 50, 5

#: ATTEMPT 2 -- written 2026-09-18 AFTER attempt 1 failed and BEFORE any attempt-2
#: fit existed.  Attempt 1's verdict is not changed by anything below.
AVERAGE_FITS = 4
ATTEMPT_2_PREREGISTRATION = {
    'written': '2026-09-18, after attempt 1 FAILED, before any attempt-2 fit ran',
    'attempt_1_verdict': ('FAIL: the_coordinates_that_already_worked_still_work, on '
                          'knee_angle_r. Worst muscle gasmed_r, peak-arm error 1.66 mm '
                          'against a bar of 1.5 x the shipped set\'s worst (0.91 mm) '
                          '= 1.36 mm. RECORDED, NOT RESCORED.'),
    'diagnosed_cause': ('Draw noise in a nondeterministic fitter, on three pieces of '
                        'evidence: gasmed_l, the mirror image fitted in the SAME run '
                        'with the SAME added coordinate, is not among the six worst '
                        'left knee muscles (so under 0.42 mm); recfem_r, '
                        'which crosses no new joint at all, moved from 0.54 to 1.42 mm; '
                        'and the pooled knee RMS over every knee muscle is unchanged '
                        '(0.64 -> 0.65 mm). None of this is a rescore -- it is why the '
                        'instrument changes, not a reason to move the bar.'),
    'instrument_change': ('(1) Muscles whose fitted coordinate set did not change keep '
                          'the SHIPPED coefficients byte for byte, so the refit cannot '
                          'move them at all; they are the paths the identified plant '
                          'runs on. (2) Muscles whose set grew get the coefficient-wise '
                          'mean of %d NEW refits; attempt 1\'s fit is excluded because it '
                          'is the draw that suggested the change. Averaging divides the '
                          'draw variance by %d and involves no selection on the outcome.'
                          % (AVERAGE_FITS, AVERAGE_FITS)),
    'why_%d' % AVERAGE_FITS: ('Cost, not measurement: each fit is ~3.3 min at 4 threads '
                              'on this shared machine, and 4 halves the draw sd.'),
    'gates': 'IDENTICAL code and bars to attempt 1 (function score()).',
    'if_it_fails': ('Recorded as FAILED, the path set is not installed, and subtalar '
                    'stays undriven in this variant.'),
}

#: ATTEMPT 3 -- written 2026-09-18 AFTER attempt 2 FAILED, BEFORE attempt 3 was
#: scored and before any plant carrying it was integrated.  Attempts 1 and 2 stay
#: FAILED; the refit route is closed in this session, as attempt 2's
#: pre-registration said it would be.
ATTEMPT_3_PREREGISTRATION = {
    'written': '2026-09-18, after attempt 2 FAILED, before attempt 3 was scored',
    'attempt_2_verdict': ('FAIL, same gate, worse: knee_angle_r worst gaslat_r 4.50 mm '
                          'against the 1.36 mm bar. RECORDED, NOT RESCORED.'),
    'what_attempt_2_revealed': ('The fitter\'s own log: of five fits, two removed 706 '
                                'and 1,170 of 1,485 sample rows as NaN (the one that '
                                'wrecked the mean fitted on ~21% of its samples, in '
                                '51 s against ~190 s). The other three removed none. '
                                'A seeded probe of 675 fitter-like poses through the '
                                'model\'s own GeometryPaths produced 0 NaN, so the '
                                'cause is specific to the fitter\'s own design and was '
                                'not isolated here. Among the three clean fits the '
                                'right gastrocnemius knee error is still bimodal '
                                '(1.66, 1.64, 0.21 mm). The fitter is not a stable '
                                'instrument for these paths, and a third refit would '
                                'be a third draw.'),
    'instrument': ('NO FIT. The shipped path set with the 22 muscles whose geometry '
                   'crosses subtalar or mtp REMOVED from it. The engine\'s '
                   'replacePathsWithFunctionBasedPaths replaces only the paths the '
                   'set names, so those 22 run on the model\'s own GeometryPaths -- '
                   'the representation this plant already uses for its 18 arm and '
                   'trunk muscles -- and the other 58 keep the shipped coefficients. '
                   '18 of the 22 are plain via-point paths; the 4 gastrocnemii carry '
                   'two wrap objects each.'),
    'what_is_NOT_evidence': ('The accuracy and recovery gates below pass BY '
                             'CONSTRUCTION for the 22: the candidate IS the ground '
                             'truth they are scored against. They are run and '
                             'reported so the file is complete, and they are not '
                             'counted as a result.'),
    'what_can_fail': ('(a) the topology gate -- nothing may acquire an arm about '
                      'a spine or wrist coordinate; (b) the engine must load and '
                      'integrate it; (c) the DYNAMIC gate G-S in '
                      'scripts/verify_articulated_spine.py, pre-registered in the '
                      'same edit: under docs/NATIVE_JOINT_LIMITS.md\'s protocol '
                      'ankle_angle_{l,r} and subtalar_angle_{l,r} must each stay '
                      'within the F2 bar, the BASE model\'s own worst excursion '
                      'measured in the same run. The wall clock per advance is '
                      'measured and reported with no bar, because no bar for it '
                      'has provenance here.'),
    'install_rule': ('Installed as muscle_paths.xml if (a) and (b) hold. G-S is a '
                     'measurement of whether the ankle regression is repaired and is '
                     'recorded either way; it does not un-install the arms, which '
                     'exist whether or not the tonic plant is quieter.'),
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def model_coordinates(path):
    """[(coordinate name, full /jointset path, default value)] in model order."""
    out = []
    for joint in ET.parse(path).getroot().find('.//JointSet/objects'):
        for coordinate in joint.iter('Coordinate'):
            out.append((coordinate.get('name'),
                        '/jointset/%s/%s' % (joint.get('name'), coordinate.get('name')),
                        float(coordinate.findtext('default_value'))))
    return out


def path_point_bodies(path):
    """{muscle: [body per path point]} straight from the model's GeometryPaths."""
    out = {}
    for force in ET.parse(path).getroot().find('.//ForceSet/objects'):
        geometry = force.find('.//GeometryPath')
        if geometry is None or force.get('name') is None:
            continue
        bodies = []
        for point in geometry.iter():
            if point.tag in ('PathPoint', 'ConditionalPathPoint', 'MovingPathPoint'):
                bodies.append(point.findtext('socket_parent_frame').split('/')[-1])
        out[force.get('name')] = bodies
    return out


def joint_topology(path):
    """{child body: (joint name, parent body, [coordinate names])}."""
    out = {}
    for joint in ET.parse(path).getroot().find('.//JointSet/objects'):
        frames = {f.get('name'): f.findtext('socket_parent').split('/')[-1]
                  for f in joint.iter('PhysicalOffsetFrame')}
        resolve = lambda s: frames.get(s.split('/')[-1], s.split('/')[-1])  # noqa: E731
        out[resolve(joint.findtext('socket_child_frame'))] = (
            joint.get('name'), resolve(joint.findtext('socket_parent_frame')),
            [c.get('name') for c in joint.iter('Coordinate')])
    return out


def chain_to_ground(topology, body):
    chain = []
    while body in topology:
        name, parent, _coordinates = topology[body]
        chain.append(name)
        body = parent
    return chain


def spans(topology, bodies, coordinate):
    """Does a path with these path-point bodies cross the joint owning `coordinate`?

    Purely topological, read off the .osim: TRUE when the path has a point on a
    body distal to the joint and a point that is not.  This is computed without
    any solver so it can be compared against the measured arm; the two agreeing
    is what makes the REPRESENTATION/TOPOLOGY split a measurement and not a
    reading of the same number twice.
    """
    joint = next(name for _child, (name, _p, cs) in topology.items() if coordinate in cs)
    distal = [b for b in set(bodies) if joint in chain_to_ground(topology, b)]
    return bool(distal) and len(set(distal)) != len(set(bodies))


def write_trajectory(destination, model_path, *, source, stride, offset, sine):
    """The reference trajectory, extended to every coordinate the variant holds.

    Columns the recorded trajectory does not carry are filled: the two foot
    coordinates with `sine`, everything else at the model's own declared default.
    Which is which is returned, so no coordinate is silently defaulted.
    """
    labels, rows = pf.read_sto(source)
    rows = rows[offset::stride]
    if not rows:
        raise ValueError('Empty trajectory selection')
    present = {label.rsplit('/', 2)[-2]: index for index, label in enumerate(labels[1:], 1)
               if label.endswith('/value')}
    amplitude, omega, phase = sine
    driven, defaulted = [], []
    for name, full, default in model_coordinates(model_path):
        if name in present:
            continue
        labels.append(full + '/value')
        if name in INHERITED + ('subtalar_angle_r', 'subtalar_angle_l'):
            driven.append(name)
            for row in rows:
                row.append(amplitude * math.sin(omega * row[0] + phase))
        else:
            defaulted.append(name)
            for row in rows:
                row.append(default)
    pf.write_sto(destination, labels, rows)
    return dict(frames=len(rows), driven=sorted(driven), defaulted=sorted(defaulted))


def arms(series, coordinate):
    """{muscle: peak |moment arm| about `coordinate|} from a sampler CSV dict."""
    out = {}
    for (actuator, quantity, column), values in series.items():
        if quantity != 'moment_arm' or column.rsplit('/', 1)[-1] != coordinate:
            continue
        out[actuator.rsplit('/', 1)[-1]] = max(abs(v) for v in values)
    return out


def load_csv(path):
    series = collections.defaultdict(list)
    with open(path) as handle:
        next(handle)
        for line in handle:
            _row, _time, actuator, quantity, column, value = line.rstrip('\n').split(',')
            series[(actuator, quantity, column)].append(float(value))
    if not series:
        raise RuntimeError('Sampler produced no rows for %s' % path)
    return series


def filter_pathset(fitted, destination, keep):
    """Write out only the paths for muscles in `keep`, and report what was dropped."""
    tree = ET.parse(fitted)
    # The fitter writes <OpenSimDocument><Set_FunctionBasedPath_><objects>; the
    # shipped file uses the same layout.  Find the one <objects> that holds paths.
    containers = [c for c in tree.getroot().iter('objects')
                  if c.find('FunctionBasedPath') is not None]
    if len(containers) != 1:
        raise RuntimeError('Expected one <objects> holding FunctionBasedPaths, found %d'
                           % len(containers))
    container = containers[0]
    dropped = []
    for element in list(container):
        name = element.get('name', '').rsplit('/', 1)[-1]
        if name not in keep:
            dropped.append(name)
            container.remove(element)
    ET.indent(tree, space='\t')
    tree.write(destination, encoding='UTF-8', xml_declaration=True)
    return sorted(dropped)


def pathset_coordinates(path):
    """{muscle: (coordinate names)} of a FunctionBasedPathSet."""
    return {element.get('name').rsplit('/', 1)[-1]:
            tuple(x.rsplit('/', 1)[-1] for x in element.findtext('coordinate_paths').split())
            for element in ET.parse(path).getroot().iter('FunctionBasedPath')}


def sampled_span(directory, coordinate):
    """The range the fitter's own Latin hypercube actually explored, in radians.

    Read out of the fitter's `*_coordinate_values_sampled.sto`, so the coverage
    claim is measured from the run rather than inferred from the bounds asked for.
    """
    files = sorted(Path(directory).glob('*_coordinate_values_sampled.sto'))
    if len(files) != 1:
        return None
    labels, rows = pf.read_sto(files[0])
    matches = [i for i, label in enumerate(labels)
               if label.endswith('/%s/value' % coordinate)]
    if len(matches) != 1:
        return None
    column = [row[matches[0]] for row in rows]
    return [min(column), max(column)]


def run(work, threads, reuse, install):
    work = Path(work)
    work.mkdir(parents=True, exist_ok=True)
    report = {'schema': 'ihm.articulated-spine-muscle-paths.v1',
              'model': str(MODEL.relative_to(ROOT)),
              'model_sha256': sha(MODEL),
              'shipped_pathset': str(SHIPPED.relative_to(ROOT)),
              'shipped_pathset_sha256': sha(SHIPPED),
              'arm_floor_m': ARM_FLOOR_M}

    topology = joint_topology(MODEL)
    bodies = path_point_bodies(MODEL)

    # ---- trajectories ----------------------------------------------------
    fit_trajectory = work / 'fit_coordinates.sto'
    fit_fill = write_trajectory(fit_trajectory, MODEL, source=COORDINATES,
                                stride=1, offset=0, sine=FOOT_SINE)
    evaluation = work / 'held_out.sto'
    eval_fill = write_trajectory(evaluation, MODEL, source=COORDINATES,
                                 stride=EVAL_STRIDE, offset=EVAL_OFFSET, sine=EVAL_SINE)
    report['trajectories'] = {
        'source': str(COORDINATES.relative_to(ROOT)),
        'fit': fit_fill, 'held_out': eval_fill,
        'fit_sine_amplitude_rad_omega_rad_s_phase': list(FOOT_SINE),
        'held_out_sine_amplitude_rad_omega_rad_s_phase': list(EVAL_SINE),
        'note': ('The reference trajectory carries 31 of this model\'s 48 '
                 'coordinates. The two foot coordinates it omits are driven by a '
                 'sine -- the upstream example\'s own treatment of mtp_angle, '
                 'reused for subtalar_angle and NOT a measured trajectory -- and '
                 'the thirteen spine and wrist coordinates are held at the '
                 'model\'s declared default, which is honest because no muscle '
                 'in this model crosses any of them (see classification).')}

    def do_sample(name, pathset=None, coordinates=None):
        out = work / (name + '.csv')
        if not (reuse and out.exists()):
            pf.sample(MODEL, coordinates or evaluation, out, pathset=pathset,
                      log=work / (name + '.samplelog'))
        return load_csv(out)

    # ---- 1. diagnosis: GeometryPath truth vs the shipped fitted paths -----
    truth = do_sample('truth-geometry')
    shipped = do_sample('shipped', SHIPPED)

    classification = {}
    for coordinate in NEW + INHERITED + CONTROLS:
        real = arms(truth, coordinate)
        fitted = arms(shipped, coordinate)
        crossing = sorted(m for m, value in real.items() if value > ARM_FLOOR_M)
        expressed = sorted(m for m, value in fitted.items() if value > ARM_FLOOR_M)
        topological = sorted(m for m, points in bodies.items()
                             if spans(topology, points, coordinate))
        if not crossing:
            verdict = 'topology'
        elif not expressed:
            verdict = 'representation'
        else:
            verdict = 'expressed'
        classification[coordinate] = {
            'verdict': verdict,
            'muscles_crossing_in_geometry': crossing,
            'muscles_with_a_fitted_arm': expressed,
            'muscles_spanning_the_joint_by_path_point_body': topological,
            # Two independent routes to the same set: path-point bodies read off
            # the XML, and the arm measured through OpenSim.  They must agree.
            'topology_agrees_with_measurement': set(topological) == set(crossing),
            'spans_but_no_measured_arm': sorted(set(topological) - set(crossing)),
            'measured_arm_but_does_not_span': sorted(set(crossing) - set(topological)),
            'peak_geometrypath_arm_m': {m: round(real[m], 6) for m in crossing},
            'peak_fitted_arm_m': max(fitted.values()) if fitted else 0.0}
    report['classification'] = classification
    report['diagnosis'] = {
        'representation': sorted(k for k, v in classification.items()
                                 if v['verdict'] == 'representation'),
        'topology': sorted(k for k, v in classification.items()
                           if v['verdict'] == 'topology'),
        'expressed': sorted(k for k, v in classification.items()
                            if v['verdict'] == 'expressed')}

    # ---- 2. attempt 1: one refit, installed whole -------------------------
    settings = dict(pf.UPSTREAM_SETTINGS)
    settings['fill_sine'] = {}
    source_muscles = set(pathset_coordinates(SHIPPED))

    def one_fit(name):
        directory = work / name
        produced = sorted(directory.glob('*_FunctionBasedPathSet.xml')) if directory.exists() else []
        if reuse and len(produced) == 1:
            return produced[0]
        started = time.time()
        result = pf.fit(MODEL, fit_trajectory, directory, settings=settings,
                        threads=threads, log=work / (name + '.log'))
        report.setdefault('fit_seconds', {})[name] = round(time.time() - started, 1)
        return result

    fitted = one_fit('refit')
    attempt1 = work / 'attempt1_muscle_paths.xml'
    dropped = filter_pathset(fitted, attempt1, source_muscles)
    covered = pathset_coordinates(attempt1)
    shipped_coordinates = pathset_coordinates(SHIPPED)
    gained = sorted(m for m in covered if set(covered[m]) != set(shipped_coordinates[m]))
    report['refit'] = {
        'settings': {k: v for k, v in settings.items() if k != 'fill_sine'},
        'paths_fitted_total': len(pathset_coordinates(fitted)),
        'paths_installed': len(covered),
        'paths_dropped_keeping_their_GeometryPath': dropped,
        'subtalar_sampled_span_rad': sampled_span(work / 'refit', 'subtalar_angle_r'),
        'muscles_whose_coordinate_set_changed': gained,
        'muscles_that_gained_subtalar_r': sorted(
            m for m, cs in covered.items() if 'subtalar_angle_r' in cs),
        'muscles_that_gained_mtp_r': sorted(
            m for m, cs in covered.items() if 'mtp_angle_r' in cs),
        'coordinates_appearing_in_any_installed_path': sorted(
            {c for cs in covered.values() for c in cs})}

    report['attempt_1'] = score(attempt1, 'refit', do_sample, truth, shipped,
                                classification, eval_fill)
    report['attempt_1']['instrument'] = (
        'ONE refit of all 80 source paths, installed whole.')

    # ---- 3. attempt 2: PRE-REGISTERED after attempt 1 FAILED ----------------
    report['attempt_2_preregistration'] = ATTEMPT_2_PREREGISTRATION
    averaged = work / 'refit_mean.xml'
    fits = [one_fit('refit-mean-%d' % k) for k in range(1, AVERAGE_FITS + 1)]
    average_pathsets(fits, averaged, source_muscles)
    attempt2 = work / 'attempt2_muscle_paths.xml'
    hybrid_pathset(SHIPPED, averaged, attempt2, keep_shipped=set(source_muscles) - set(gained))
    report['attempt_2'] = score(attempt2, 'attempt2', do_sample, truth, shipped,
                                classification, eval_fill)
    report['attempt_2']['instrument'] = (
        'HYBRID. The %d muscles whose fitted coordinate set did not change keep '
        'the shipped coefficients byte for byte; the %d whose set grew get the '
        'coefficient-wise mean of %d NEW refits (attempt 1 excluded).'
        % (len(source_muscles) - len(gained), len(gained), AVERAGE_FITS))
    # The run-to-run spread of the statistic that failed, measured over every
    # independent fit this run made.  This is the noise measurement the
    # attempt-1 diagnosis rests on, and it is reported whatever it says.
    spread = {}
    for index, path in enumerate([fitted] + fits):
        single = work / ('single_%d.xml' % index)
        filter_pathset(path, single, source_muscles)
        series = do_sample('single-%d' % index, single)
        spread['fit_%d' % index] = worst_errors(series, truth, ('knee_angle_r', 'knee_angle_l'))
    report['fitter_spread_on_the_failed_statistic'] = spread

    # ---- 4. attempt 3: PRE-REGISTERED after attempt 2 FAILED -- no fit ------
    report['attempt_3_preregistration'] = ATTEMPT_3_PREREGISTRATION
    attempt3 = work / 'attempt3_muscle_paths.xml'
    filter_pathset(SHIPPED, attempt3, set(source_muscles) - set(gained))
    report['attempt_3'] = score(attempt3, 'attempt3', do_sample, truth, shipped,
                                classification, eval_fill)
    report['attempt_3']['instrument'] = (
        'NO FIT: the shipped set minus the %d muscles whose geometry crosses '
        'subtalar or mtp, which therefore run on the model\'s own GeometryPaths.'
        % len(gained))
    report['attempt_3']['muscles_on_their_own_GeometryPath'] = gained
    report['attempt_3']['accuracy_gates_are_evidence'] = False
    report['attempt_3']['gate_that_can_fail_here'] = 'the_topology_class_is_still_zero'

    report['passed'] = report['attempt_3']['gates']['the_topology_class_is_still_zero']['passed']
    report['verdict'] = {
        'attempt_1': 'PASS' if report['attempt_1']['passed'] else 'FAIL',
        'attempt_2': 'PASS' if report['attempt_2']['passed'] else 'FAIL',
        'attempt_3': ('PASS (accuracy gates by construction; see preregistration)'
                      if report['attempt_3']['passed'] else 'FAIL')}

    if install:
        if not report['passed']:
            raise SystemExit('Refusing to install a path set that failed its gates')
        INSTALLED.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(attempt3, INSTALLED)
        report['installed'] = {'path': str(INSTALLED.relative_to(ROOT)),
                               'sha256': sha(INSTALLED), 'attempt': 3}
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=2, default=str) + '\n')
    return report


def worst_errors(series, truth, coordinates):
    """{coordinate: (worst muscle, |peak refit arm - peak true arm|)}."""
    out = {}
    for coordinate in coordinates:
        real, got = arms(truth, coordinate), arms(series, coordinate)
        crossing = [m for m, v in real.items() if v > ARM_FLOOR_M]
        worst = max(((abs(got.get(m, 0.) - real[m]), m) for m in crossing), default=(0., None))
        out[coordinate] = {'muscle': worst[1], 'error_m': worst[0],
                           'gasmed_error_m': abs(got.get('gasmed_' + coordinate[-1], 0.)
                                                 - real.get('gasmed_' + coordinate[-1], 0.))}
    return out


def average_pathsets(files, destination, keep):
    """Coefficient-wise mean of several fitted path sets.

    Exact, not an approximation: every path here is ONE order-5
    MultivariatePolynomialFunction in its length, moment arms are its
    derivatives, and both are linear in the coefficients -- so the mean of the
    coefficient vectors IS the mean polynomial.  Refuses unless every file
    agrees on each muscle's coordinate list, dimension and order.
    """
    trees = [ET.parse(f) for f in files]
    index = []
    for tree in trees:
        index.append({e.get('name').rsplit('/', 1)[-1]: e
                      for e in tree.getroot().iter('FunctionBasedPath')})
    names = set(index[0]) & set(keep)
    for other in index[1:]:
        if set(other) & set(keep) != names:
            raise ValueError('Fits disagree on which muscles they cover')
    for name in names:
        signatures = {(e.findtext('coordinate_paths').split().__repr__(),
                       e.findtext('length_function/MultivariatePolynomialFunction/dimension'),
                       e.findtext('length_function/MultivariatePolynomialFunction/order'))
                      for e in (d[name] for d in index)}
        if len(signatures) != 1:
            raise ValueError('Fits disagree on the polynomial structure of %s' % name)
        vectors = [[float(x) for x in d[name].findtext(
            'length_function/MultivariatePolynomialFunction/coefficients').split()]
            for d in index]
        if len({len(v) for v in vectors}) != 1:
            raise ValueError('Coefficient count differs for %s' % name)
        mean = [sum(column) / len(column) for column in zip(*vectors)]
        index[0][name].find('length_function/MultivariatePolynomialFunction/coefficients').text = \
            ' '.join(repr(v) for v in mean)
    filter_pathset_tree(trees[0], destination, keep)


def filter_pathset_tree(tree, destination, keep):
    containers = [c for c in tree.getroot().iter('objects')
                  if c.find('FunctionBasedPath') is not None]
    if len(containers) != 1:
        raise RuntimeError('Expected one <objects> holding FunctionBasedPaths')
    for element in list(containers[0]):
        if element.get('name', '').rsplit('/', 1)[-1] not in keep:
            containers[0].remove(element)
    ET.indent(tree, space='\t')
    tree.write(destination, encoding='UTF-8', xml_declaration=True)


def hybrid_pathset(shipped, fitted, destination, *, keep_shipped):
    """Shipped paths for `keep_shipped`, fitted paths for every other muscle.

    The shipped elements are copied as parsed XML, so their coefficient TEXT is
    carried unchanged (indentation and XML comments are not); the verifier
    re-reads both files and checks every kept path's coordinate list and
    coefficient string for string equality.
    """
    tree = ET.parse(fitted)
    shipped_paths = {e.get('name').rsplit('/', 1)[-1]: e
                     for e in ET.parse(shipped).getroot().iter('FunctionBasedPath')}
    containers = [c for c in tree.getroot().iter('objects')
                  if c.find('FunctionBasedPath') is not None]
    container = containers[0]
    for position, element in enumerate(list(container)):
        name = element.get('name').rsplit('/', 1)[-1]
        if name in keep_shipped:
            container.remove(element)
            container.insert(position, shipped_paths[name])
    missing = set(keep_shipped) - {e.get('name').rsplit('/', 1)[-1] for e in container}
    if missing:
        raise ValueError('Hybrid lost muscles: %r' % sorted(missing))
    ET.indent(tree, space='\t')
    tree.write(destination, encoding='UTF-8', xml_declaration=True)


def score(candidate, name, do_sample, truth, shipped, classification, eval_fill):
    """The four gates, unchanged between attempts.  Written before attempt 1 ran."""
    refit = do_sample(name, candidate)
    length = pf.compare(refit, truth, 'length')
    shipped_length = pf.compare(shipped, truth, 'length')
    arm_error = pf.compare(refit, truth, 'moment_arm')
    shipped_arm_error = pf.compare(shipped, truth, 'moment_arm')
    gates = {}
    gates['accuracy_against_the_models_own_GeometryPaths'] = {
        'note': ('Held-out frames, refit and shipped set scored against the SAME '
                 'truth. A raw RMS alone would say nothing; the shipped set is '
                 'the baseline and the refit must not be worse than 1.5x it.'),
        'held_out_frames': eval_fill['frames'],
        'refit_length_rms_m': length['rms_m'],
        'shipped_length_rms_m': shipped_length['rms_m'],
        'refit_moment_arm_rms_m': arm_error['rms_m'],
        'shipped_moment_arm_rms_m': shipped_arm_error['rms_m'],
        'refit_worst_actuator': length['worst_actuator'],
        'refit_worst_actuator_rms_m': length['max_actuator_rms_m'],
        'passed': (length['rms_m'] <= 1.5 * shipped_length['rms_m']
                   and arm_error['rms_m'] <= 1.5 * shipped_arm_error['rms_m'])}

    recovered = {}
    for coordinate in ('subtalar_angle_r', 'subtalar_angle_l',
                       'mtp_angle_r', 'mtp_angle_l'):
        real = arms(truth, coordinate)
        got = arms(refit, coordinate)
        crossing = classification[coordinate]['muscles_crossing_in_geometry']
        worst = max(((abs(got.get(m, 0.0) - real[m]), m) for m in crossing),
                    default=(0.0, None))
        recovered[coordinate] = {
            'muscles_expected': len(crossing),
            'muscles_now_nonzero': sum(1 for m in crossing if got.get(m, 0.0) > ARM_FLOOR_M),
            'worst_peak_arm_error_m': worst[0], 'worst_muscle': worst[1],
            'per_muscle_m': {m: {'geometrypath': round(real[m], 6),
                                 'refit': round(got.get(m, 0.0), 6)} for m in crossing}}
    gates['the_representation_class_is_recovered'] = {
        'note': ('Every muscle the GeometryPath truth says crosses subtalar or '
                 'mtp must now carry a nonzero fitted arm, and its peak must sit '
                 'within 5 mm of the truth. 5 mm is not a physiological bar: it '
                 'is the scale at which "the arm is there" stops being the claim '
                 'and "the arm is right" starts, and the measured worst is '
                 'reported beside it.'),
        'coordinates': recovered,
        'passed': all(v['muscles_now_nonzero'] == v['muscles_expected']
                      and v['worst_peak_arm_error_m'] < 5e-3 for v in recovered.values())}

    unchanged = {}
    for coordinate in CONTROLS:
        real, got, was = arms(truth, coordinate), arms(refit, coordinate), arms(shipped, coordinate)
        crossing = [m for m, v in real.items() if v > ARM_FLOOR_M]
        worst = max(((abs(got.get(m, 0.) - real[m]), m) for m in crossing), default=(0., None))
        unchanged[coordinate] = {
            'muscles': len(crossing),
            'refit_worst_error_m': worst[0], 'refit_worst_muscle': worst[1],
            'shipped_worst_error_m': max((abs(was.get(m, 0.) - real[m]) for m in crossing), default=0.)}
    gates['the_coordinates_that_already_worked_still_work'] = {
        'note': ('A refit that recovered subtalar by wrecking the ankle would '
                 'pass the gate above. These are coordinates the shipped set '
                 'already carried; the refit must not be worse than it was.'),
        'coordinates': unchanged,
        'passed': all(v['refit_worst_error_m'] <= 1.5 * max(v['shipped_worst_error_m'], 1e-4)
                      for v in unchanged.values())}

    topology = [c for c in NEW if classification[c]['verdict'] == 'topology']
    gates['the_topology_class_is_still_zero'] = {
        'note': ('The control that can fail in the other direction. The nine '
                 'spine and four wrist coordinates have NO muscle crossing them '
                 'in the geometry, so a refit that reported an arm about one of '
                 'them would be fitting noise, and every number in this file '
                 'would be suspect.'),
        'coordinates': {c: {'peak_refit_arm_m': max(arms(refit, c).values(), default=0.0),
                            'peak_truth_arm_m': max(arms(truth, c).values(), default=0.0)}
                        for c in topology},
        'passed': all(max(arms(refit, c).values(), default=0.0) <= ARM_FLOOR_M
                      for c in topology)}
    return {'gates': gates, 'passed': all(g['passed'] for g in gates.values()),
            'failed_gates': sorted(k for k, g in gates.items() if not g['passed'])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--work', default=str(ROOT / 'data/derived/articulated-spine-muscle-paths/work'))
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--reuse', action='store_true')
    parser.add_argument('--install', action='store_true')
    args = parser.parse_args()
    report = run(args.work, args.threads, args.reuse, args.install)
    print(json.dumps({k: v for k, v in report.items() if k != 'classification'},
                     indent=2, default=str))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
