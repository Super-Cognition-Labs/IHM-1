"""Derive data/models/thoracic_drive_v1 from data/models/shoulder_girdle_v1.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_thoracic_drive
    OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_thoracic_drive --check

WHAT THIS CHANGES, and it is a small list.

`shoulder_girdle_v1` put EVERY torso-side attachment of the sixty Seth 2019
donor muscles on the body `torso`, because the girdle map is expressed in the
torso frame.  `articulated_spine_v1` had already split the rib cage off `torso`
into a separate `thorax` body carrying the `thoracic` joint.  So the girdle's
rib-cage attachments ended up on the lumbar-side body, and the rib cage's own
joint had no muscle crossing it -- which is the state
docs/ARTICULATED_SPINE.md records as BLOCKED.

`scripts/measure_thoracic_reassignability.py` asked, under a rule fixed before
it ran, which of those attachments belong on the rib cage.  Ten muscles, five
per side, passed both halves: their nearest anatomical surface in this body's
own meshes is a member of the thorax partition, and they sit more than
sqrt(2) x 62 mm above the thoracic joint centre.  This builder moves exactly
those ten, and the rib-cage wrap surface they all use.

  1. The ten muscles' torso attachments move from `torso` to `thorax`.
  2. The `Thorax_r` and `Thorax_l` wrap ellipsoids move from `torso` to
     `thorax`.

Nothing else changes.  No mass moves, no body is added, no coordinate is added,
no muscle parameter is touched, and no number is introduced that is not already
in the base model.

WHY THE ARITHMETIC IS EXACT.  The `thoracic` joint is `torso -> thorax`.  Its
parent frame is a PURE TRANSLATION on `torso` with zero orientation and its child
frame is the `thorax` origin, so at `thoracic_* = 0` a point at `p` in torso
coordinates is at `p - t` in thorax coordinates.  The builder asserts both
frames before using them.  There is no registration, no similarity factor and
nothing fitted here: the known answer is that every path length is UNCHANGED at
the rest pose, and `scripts/verify_thoracic_drive.py` measures it through the
engine instead of trusting this paragraph.

WHY THE WRAP MOVES TOO, and what it costs.  `Thorax_{r,l}` is the rib-cage
ellipsoid.  Leaving it on `torso` while its attachments ride `thorax` would put
the rib cage's SURFACE and the rib cage's ATTACHMENTS on two different bodies,
and at the declared thoracic range limit the two would disagree by about 44 mm
(the ellipsoid centre sits ~0.17 m from the joint centre and the range is
0.262 rad).  The repository has already made this call once in the same
direction: `hat_ribs_scap.vtp` was moved to `thorax` by
`scripts/build_articulated_spine.py`.

  The cost, stated because it is not small: NINETEEN muscles per side wrap on
  `Thorax_{r,l}`, not five.  The other fourteen -- trapezius, rhomboid, levator
  scapulae, latissimus dorsi, teres major, supraspinatus, pectoralis major
  clavicular -- keep their attachments on `torso` (the measurement refused them:
  their endpoints are on VERTEBRAE, which stayed in the torso residual) and yet
  acquire a dependence on `thoracic_*` through the wrap surface alone.  That is
  a real mechanism and a much weaker one than an attachment, and
  docs/WRIST_AND_THORAX.md reports the two populations separately rather than
  summing them into one "muscles that drive the thoracic joint" count.

RANGES AND STOPS.  No coordinate is added, so nothing new is needed.
`thoracic_extension`, `thoracic_bending` and `thoracic_rotation` already carry a
declared `<range>` and an `articulated_stop_*` `CoordinateLimitForce` at exactly
that range, inherited unchanged from `articulated_spine_v1`.  The builder
asserts all three are still present and still at their declared range, because
docs/NATIVE_JOINT_LIMITS.md says nothing else enforces a range here.

WHAT THIS DOES NOT DO.  It does not make the thoracic joint controllable: five
muscles per side, all of them anterior chest muscles, is not a bidirectional
trunk.  The measured bound is in docs/WRIST_AND_THORAX.md and it is a capacity
at one pose, not behaviour.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import pathlib
import shutil
import sys
import xml.etree.ElementTree as ET

ROOT = pathlib.Path(__file__).resolve().parents[1]

BASE_DIR = ROOT / 'data/models/shoulder_girdle_v1'
OUT_DIR = ROOT / 'data/models/thoracic_drive_v1'
REPORT = ROOT / 'out/thoracic_reassignability.json'

#: The rib-cage wrap ellipsoids.  Named here and asserted to exist; the builder
#: refuses to run if either is missing or has moved.
RIB_CAGE_WRAPS = ('Thorax_r', 'Thorax_l')

THORACIC_COORDINATES = ('thoracic_extension', 'thoracic_bending', 'thoracic_rotation')

CARRIED = ('catalog.json', 'equilibrium_excitations.json', 'initial_pose.json')


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def triple(text: str):
    return [float(x) for x in text.split()]


def write_triple(values):
    return ' '.join(repr(float(v)) for v in values)


def thoracic_offset(root: ET.Element):
    """The pure translation that carries a torso point into the thorax frame.

    Asserted, not assumed: both offset frames must have zero orientation and the
    child frame must sit at the thorax origin, or the subtraction below is not
    the right transform and the builder refuses to write anything.
    """
    for element in root.iter():
        if not (element.tag.endswith('Joint') and element.get('name') == 'thoracic'):
            continue
        parent, child = None, None
        for frame in element.iter('PhysicalOffsetFrame'):
            body = frame.findtext('socket_parent').rstrip('/').split('/')[-1]
            entry = (triple(frame.findtext('translation')),
                     triple(frame.findtext('orientation')))
            if body == 'torso':
                parent = entry
            elif body == 'thorax':
                child = entry
        if parent is None or child is None:
            raise SystemExit('the thoracic joint does not connect torso to thorax')
        if max(abs(v) for v in parent[1]) != 0.0:
            raise SystemExit('the thoracic parent frame is rotated; a translation '
                             'is not the right transform')
        if max(abs(v) for v in child[0]) != 0.0 or max(abs(v) for v in child[1]) != 0.0:
            raise SystemExit('the thoracic child frame is not the thorax origin')
        coordinates = [c.get('name') for c in element.iter('Coordinate')]
        if sorted(coordinates) != sorted(THORACIC_COORDINATES):
            raise SystemExit('unexpected thoracic coordinates: %r' % coordinates)
        return parent[0]
    raise SystemExit('no thoracic joint in the base model')


def assert_stops(root: ET.Element):
    """Every thoracic coordinate still declares a range and still carries a stop
    at exactly that range.  docs/NATIVE_JOINT_LIMITS.md: nothing else enforces it."""
    ranges = {}
    for coordinate in root.iter('Coordinate'):
        if coordinate.get('name') in THORACIC_COORDINATES:
            ranges[coordinate.get('name')] = triple(coordinate.findtext('range'))
    stops = {}
    for force in root.iter():
        if not force.tag.endswith('CoordinateLimitForce'):
            continue
        name = force.findtext('coordinate')
        if name in THORACIC_COORDINATES:
            stops[name] = (float(force.findtext('lower_limit')),
                           float(force.findtext('upper_limit')))
    # `CoordinateLimitForce` reads its limits in OpenSim's DEGREE convention
    # while `<range>` is in radians, and `scripts/native_mechanical_stream.cpp`
    # converts them on the way in. docs/ARTICULATED_SPINE.md records the same
    # convention (and the damping term that is NOT converted, which is a
    # separate and unfixed disagreement). Compare in degrees or the assertion
    # fires on a stop that is correctly placed.
    to_degrees = 180.0 / 3.141592653589793
    for name in THORACIC_COORDINATES:
        if name not in ranges:
            raise SystemExit('%s declares no range' % name)
        if name not in stops:
            raise SystemExit('%s carries no CoordinateLimitForce' % name)
        lower, upper = (v * to_degrees for v in ranges[name])
        if abs(stops[name][0] - lower) > 1e-9 or abs(stops[name][1] - upper) > 1e-9:
            raise SystemExit('%s stop %r deg is not at its declared range %r deg'
                             % (name, stops[name], (lower, upper)))
    return {k: list(v) for k, v in ranges.items()}


def reassignable_muscles():
    """The list is READ from the pre-registered measurement, never typed here."""
    if not REPORT.exists():
        raise SystemExit('run scripts.measure_thoracic_reassignability first; '
                         'this builder does not decide which muscles move')
    report = json.loads(REPORT.read_text())
    if report.get('schema') != 'ihm.thoracic-reassignability.v1':
        raise SystemExit('unexpected report schema %r' % report.get('schema'))
    if report.get('verdict') != 'UNBLOCKED':
        raise SystemExit('the measurement says %r; nothing to build'
                         % report.get('verdict'))
    verdicts = {m['muscle']: m['verdict'] for m in report['muscles']}
    muscles = sorted(report['reassignable'])
    if sorted(m for m, v in verdicts.items() if v == 'REASSIGNABLE') != muscles:
        raise SystemExit('the report disagrees with itself about which muscles '
                         'are reassignable')
    return muscles, report


def build(root: ET.Element, muscles, offset):
    moved_points, moved_wraps = [], []
    targets = set(muscles)

    for muscle in root.iter():
        if muscle.get('name') not in targets or not muscle.tag.endswith('Muscle'):
            continue
        for point in muscle.iter('PathPoint'):
            parent = point.findtext('socket_parent_frame')
            if not (parent and parent.rstrip('/').endswith('torso')):
                continue
            before = triple(point.findtext('location'))
            after = [before[i] - offset[i] for i in range(3)]
            point.find('socket_parent_frame').text = '/bodyset/thorax'
            point.find('location').text = write_triple(after)
            moved_points.append({'muscle': muscle.get('name'),
                                 'point': point.get('name'),
                                 'from_body': 'torso', 'to_body': 'thorax',
                                 'location_torso_m': before,
                                 'location_thorax_m': after})

    torso_set = thorax_set = None
    for body in root.iter('Body'):
        if body.get('name') == 'torso':
            torso_set = body.find('WrapObjectSet').find('objects')
        elif body.get('name') == 'thorax':
            thorax_set = body.find('WrapObjectSet').find('objects')
    if torso_set is None or thorax_set is None:
        raise SystemExit('torso or thorax has no WrapObjectSet')

    for name in RIB_CAGE_WRAPS:
        found = [w for w in list(torso_set) if w.get('name') == name]
        if len(found) != 1:
            raise SystemExit('expected exactly one %s on torso, found %d'
                             % (name, len(found)))
        wrap = found[0]
        before = triple(wrap.findtext('translation'))
        after = [before[i] - offset[i] for i in range(3)]
        wrap.find('translation').text = write_triple(after)
        torso_set.remove(wrap)
        thorax_set.append(wrap)
        moved_wraps.append({'wrap_object': name, 'from_body': 'torso',
                            'to_body': 'thorax',
                            'translation_torso_m': before,
                            'translation_thorax_m': after})

    if len(moved_wraps) != len(RIB_CAGE_WRAPS):
        raise SystemExit('not every rib-cage wrap moved')
    return moved_points, moved_wraps


def rewrite_catalog(muscles, moved_points, report_sha):
    rows = json.loads((BASE_DIR / 'catalog.json').read_text())
    touched = {m['muscle'] for m in moved_points}
    for row in rows:
        if row.get('id') not in touched:
            continue
        row['attachment_bodies'] = sorted(
            'thorax' if b == 'torso' else b for b in row['attachment_bodies'])
        crossed = set(row.get('crosses_joints') or [])
        crossed.add('thoracic')
        row['crosses_joints'] = sorted(crossed)
        row['reassignment_basis'] = (
            'Torso attachment moved to the thorax body by '
            'scripts/build_thoracic_drive.py. The transform is the thoracic '
            'joint offset subtracted exactly; no registration and no fitted '
            'quantity. Membership of the thorax was decided by '
            'scripts/measure_thoracic_reassignability.py under a rule fixed '
            'before it ran (nearest anatomical surface in this body\'s own '
            'meshes is a thorax-partition structure, AND the attachment sits '
            'more than sqrt(2) x 62 mm above the thoracic joint centre). '
            'Report sha256 %s. The muscle\'s donor parameters are unchanged.'
            % report_sha)
    return rows, sorted(touched)


def emit(out_dir: pathlib.Path, tree: ET.ElementTree, rows, registration):
    out_dir.mkdir(parents=True, exist_ok=True)
    tree.write(out_dir / 'model.osim', encoding='UTF-8', xml_declaration=True)
    (out_dir / 'model.osim').write_bytes(
        (out_dir / 'model.osim').read_bytes().rstrip() + b'\n')
    (out_dir / 'catalog.json').write_text(json.dumps(rows, indent=2) + '\n')
    for name in CARRIED:
        if name == 'catalog.json':
            continue
        shutil.copyfile(BASE_DIR / name, out_dir / name)
    shutil.copyfile(REPORT, out_dir / 'thoracic_reassignability.json')
    registration['model_sha256'] = sha256(out_dir / 'model.osim')
    registration['catalog_sha256'] = sha256(out_dir / 'catalog.json')
    (out_dir / 'registration.json').write_text(json.dumps(registration, indent=2) + '\n')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true',
                        help='build into a scratch directory and require the '
                             'bytes to equal what is already installed')
    args = parser.parse_args(argv)

    muscles, report = reassignable_muscles()
    report_sha = sha256(REPORT)

    base_model = BASE_DIR / 'model.osim'
    tree = ET.parse(base_model)
    root = tree.getroot()

    offset = thoracic_offset(root)
    ranges = assert_stops(root)

    # The round-trip reference: the base model parsed and written back out with
    # no change at all.  Everything that is not an intended change is then a
    # byte-for-byte match against it, and `diff` shows exactly the edit.
    reference = copy.deepcopy(tree)

    moved_points, moved_wraps = build(root, muscles, offset)
    if not moved_points:
        raise SystemExit('nothing moved; the base model is not what was expected')
    rows, touched = rewrite_catalog(muscles, moved_points, report_sha)

    registration = {
        'schema': 'ihm.thoracic-drive-variant.v1',
        'model_path': 'data/models/thoracic_drive_v1/model.osim',
        'model_sha256': None,
        'catalog_path': 'data/models/thoracic_drive_v1/catalog.json',
        'catalog_sha256': None,
        'base_model_path': 'data/models/shoulder_girdle_v1/model.osim',
        'base_model_sha256': sha256(base_model),
        'sources': {
            'data/models/shoulder_girdle_v1/model.osim': sha256(base_model),
            'data/models/shoulder_girdle_v1/catalog.json': sha256(BASE_DIR / 'catalog.json'),
            'data/models/shoulder_girdle_v1/registration.json':
                sha256(BASE_DIR / 'registration.json'),
            'scripts/build_thoracic_drive.py': sha256(pathlib.Path(__file__)),
            'scripts/measure_thoracic_reassignability.py':
                sha256(ROOT / 'scripts/measure_thoracic_reassignability.py'),
            'out/thoracic_reassignability.json': report_sha,
        },
        'muscle_count': 158,
        'body_count': 29,
        'default_enabled': False,
        'native_acceptance_complete': False,
        'native_verified': False,
        'donor': {
            'note': 'No new donor. Every muscle, parameter and path point here is '
                    'already in data/models/shoulder_girdle_v1 and carries the '
                    'provenance that build recorded: Seth, Dong, Matias and Delp '
                    '2019, Front. Neurorobot. 13:90, CC BY 4.0 from the authors on '
                    'SimTK and Apache-2.0 as a file of opensim-core. '
                    'docs/SHOULDER_GIRDLE.md section 1. NOTHING here comes from '
                    'MoBL-ARMS 4.1 or from the Gonzalez/Buchanan/Delp wrist model, '
                    'directly or indirectly.',
        },
        'changes': {
            'muscle_attachments_moved_torso_to_thorax': moved_points,
            'wrap_objects_moved_torso_to_thorax': moved_wraps,
            'thoracic_joint_offset_m': list(offset),
            'transform': 'p_thorax = p_torso - thoracic_joint_offset. The thoracic '
                         'joint parent frame is a pure translation on torso with '
                         'zero orientation and its child frame is the thorax '
                         'origin, both asserted by the builder, so this is exact '
                         'at thoracic_* = 0 and every path length is unchanged '
                         'at the rest pose.',
            'muscles_reassigned': touched,
            'mass_moved_kg': 0.0,
            'bodies_added': 0,
            'coordinates_added': 0,
            'muscle_parameters_changed': 0,
        },
        'thoracic_ranges_rad': ranges,
        'thoracic_stops': 'inherited unchanged from articulated_spine_v1; one '
                          'CoordinateLimitForce per thoracic coordinate, at '
                          'exactly its declared range, asserted by the builder',
        'selection': "NativeMechanicalStream(..., augmented_registration="
                     "'data/models/thoracic_drive_v1/registration.json')",
        'limitations': [
            'Five muscles per side is not a controllable trunk. All ten are '
            'anterior chest muscles; the measured torque bound and its sign are '
            'in docs/WRIST_AND_THORAX.md and are a capacity at one pose, not '
            'behaviour.',
            'Nineteen muscles per side wrap on the rib-cage ellipsoid that moved, '
            'so fourteen per side acquire a dependence on thoracic_* through the '
            'WRAP SURFACE alone while their attachments stay on torso. That is a '
            'weaker mechanism than an attachment and is reported separately.',
            'The thoracic joint centre is inherited from a registration with a '
            '62 mm RMS proxy residual, and two of the five reassignments '
            '(PectoralisMajorThorax_M, SerratusAnterior_M) won their nearest-'
            'surface test by 8-15 mm, well inside that residual. They passed the '
            'rule as written and were not rescored; the margins are reported.',
            'No linearization.npz, no stance acceptance, no controller, and no '
            'actuator port on any thoracic coordinate.',
        ],
    }

    if args.check:
        scratch = ROOT / 'data/derived/thoracic-drive-check'
        if scratch.exists():
            shutil.rmtree(scratch)
        emit(scratch, tree, rows, dict(registration))
        bad = []
        for name in sorted(p.name for p in scratch.iterdir()):
            a, b = scratch / name, OUT_DIR / name
            if not b.exists() or a.read_bytes() != b.read_bytes():
                bad.append(name)
        shutil.rmtree(scratch)
        if bad:
            print('DIFFERS from the installed build: %s' % ', '.join(bad))
            return 1
        print('--check: the builder is a function; %d files identical'
              % len(list(OUT_DIR.iterdir())))
        return 0

    reference_path = ROOT / 'data/derived/thoracic-drive-reference.osim'
    reference_path.parent.mkdir(parents=True, exist_ok=True)
    reference.write(reference_path, encoding='UTF-8', xml_declaration=True)
    reference_path.write_bytes(reference_path.read_bytes().rstrip() + b'\n')

    emit(OUT_DIR, tree, rows, registration)
    print('wrote %s' % OUT_DIR.relative_to(ROOT))
    print('  muscles reassigned : %d (%s)' % (len(touched), ', '.join(touched)))
    print('  path points moved  : %d' % len(moved_points))
    print('  wrap objects moved : %d' % len(moved_wraps))
    print('  thoracic offset    : %r' % (offset,))
    print('  round-trip reference for diffing: %s'
          % reference_path.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
