"""What data/models/shoulder_girdle_v1 is, checked rather than asserted.

    OPENBLAS_NUM_THREADS=1 nice -n 10 \\
        .venv/bin/python -m unittest scripts.verify_shoulder_girdle -v

Two native sessions dominate the wall clock. ONE at a time, and the stream
prlimit-caps itself at 4 GB.

The checks are in three groups, and the ones that matter are the ones that CAN
fail for the reason they are about:

  Partition  the torso was repartitioned and not added to, and the five parts
             reassemble the one body they came from.
  File       every new coordinate declares a range, its rest value is inside it,
             and a CoordinateLimitForce sits at exactly that range -- and no
             coordinate that already existed gained one.
  Native     the engine's inventory; the arm did NOT move when it was re-parented
             onto the scapula; each new joint moves only its distal bodies; the
             left girdle is the exact z-mirror of the right; and the measured
             moment arms of the new muscles about the new coordinates, with the
             two controls this repository already uses -- a repeat call must be
             bit-equal, and soleus_r about ankle_angle_r must still read
             -0.0497 m.
"""
from __future__ import annotations

import json
import shutil
import sys
import unittest
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ihm.assembly.cervical_inertia import combine_bodies  # noqa: E402
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402
from scripts.build_shoulder_girdle import (  # noqa: E402
    BASE_MODEL, GIRDLE_RANGE_HALFWIDTH_RAD, MASS_SCALE, OUT,
    STOP_STIFFNESS_NM_PER_RAD, build, inertia_matrix, numbers, text_of,
)

TARGET_MASS_KG = 77.6122029
DEGREES = 180.0 / np.pi
MODEL = ROOT / OUT / 'model.osim'
REGISTRATION = OUT + '/registration.json'
RECORD = json.loads((ROOT / OUT / 'registration.json').read_text())
NEW_RANGES = {k: tuple(v) for k, v in RECORD['new_coordinates'].items()}
NEW_BODIES = ('clavicle_r', 'scapula_r', 'clavicle_l', 'scapula_l')
MIRROR = np.diag([1.0, 1.0, -1.0])

#: What each new joint is ALLOWED to move, and bodies it must not. The second
#: half is what makes this a control that can fail.
DISTAL = {
    'clav_elev_r': (('clavicle_r',), ('torso', 'pelvis', 'clavicle_l', 'thorax')),
    'scapula_upward_rot_r': (('scapula_r', 'humerus_r', 'hand_r'),
                             ('torso', 'pelvis', 'scapula_l', 'humerus_l')),
    'scapula_abduction_l': (('scapula_l', 'humerus_l'),
                            ('torso', 'pelvis', 'scapula_r', 'humerus_r')),
}


def coordinates(path):
    found = {}
    for element in ET.parse(path).getroot().iter('Coordinate'):
        text = element.findtext('range')
        found[element.get('name')] = (
            None if text is None else tuple(float(v) for v in text.split()),
            float(element.findtext('default_value')))
    return found


def masses(path):
    root = ET.parse(path).getroot()
    return {b.get('name'): float(b.findtext('mass'))
            for b in root.findall('.//BodySet/objects/Body')}


def body_records(path, names):
    out = {}
    for body in ET.parse(path).getroot().findall('.//BodySet/objects/Body'):
        if body.get('name') in names:
            out[body.get('name')] = {
                'mass_kg': float(text_of(body, 'mass')),
                'center_m': numbers(body, 'mass_center').tolist(),
                'inertia_kg_m2': inertia_matrix(numbers(body, 'inertia')).tolist()}
    return out


class Partition(unittest.TestCase):
    """The girdle came OUT of the torso."""

    def test_total_mass_conserved(self):
        base = masses(ROOT / BASE_MODEL)
        variant = masses(MODEL)
        self.assertEqual(len(base), 25)
        self.assertEqual(len(variant), 29)
        self.assertLess(abs(sum(variant.values()) - sum(base.values())), 1e-9)

    def test_only_torso_changed(self):
        base = masses(ROOT / BASE_MODEL)
        variant = masses(MODEL)
        for name, value in base.items():
            if name == 'torso':
                continue
            self.assertEqual(value, variant[name], name)
        debited = base['torso'] - variant['torso']
        added = sum(variant[n] for n in NEW_BODIES)
        self.assertLess(abs(debited - added), 1e-12)

    def test_the_partition_recombines_to_the_original_torso(self):
        """The strongest known answer available: mass, COM and the FULL inertia
        tensor of the five parts must reassemble the one body they came from.

        The new bodies' records are read back out of the written model file and
        carried into the torso frame by their own joints' declared placement, so
        this is a check on the FILE, not a replay of the builder's arithmetic.
        """
        base = body_records(ROOT / BASE_MODEL, {'torso'})['torso']
        variant = body_records(MODEL, set(NEW_BODIES) | {'torso'})
        parts = [variant['torso']]
        frames = self.placement()
        for name in NEW_BODIES:
            rotation, origin = frames[name]
            record = variant[name]
            centre = rotation @ np.asarray(record['center_m'], float) + origin
            parts.append({'mass_kg': record['mass_kg'], 'center_m': centre.tolist(),
                          'inertia_kg_m2': (rotation @ np.asarray(record['inertia_kg_m2'],
                                                                 float) @ rotation.T).tolist()})
        total = combine_bodies(parts)
        self.assertLess(abs(total['mass_kg'] - base['mass_kg']), 1e-12)
        self.assertLess(np.abs(np.asarray(total['center_m'])
                               - np.asarray(base['center_m'])).max(), 1e-12)
        self.assertLess(np.linalg.norm(np.asarray(total['inertia_kg_m2'])
                                       - np.asarray(base['inertia_kg_m2'])), 1e-12)

    def placement(self):
        """Each new body's (rotation, origin) in the torso frame at the rest pose.

        Read from the measured donor frames and the registration's own map, which
        is what the model file was written from; the NATIVE mirror test below is
        the independent check that the file agrees with the engine.
        """
        frames = json.loads((ROOT / OUT / 'donor_reference_frames.json').read_text())
        out = {}
        for name in NEW_BODIES:
            donor_body = name.split('_')[0]
            transform = np.asarray(frames['frames'][donor_body]['body_to_thorax_reference'],
                                   float)
            side = 1 if name.endswith('_r') else -1
            rotation = transform[:3, :3]
            origin = (RECORD['map']['s_lat'] * transform[:3, 3]
                      + np.asarray(RECORD['map']['delta_torso_frame_m'], float))
            if side < 0:
                rotation = MIRROR @ rotation @ MIRROR
                origin = MIRROR @ origin
            out[name] = (rotation, origin)
        return out

    def test_the_new_body_inertias_are_physically_admissible(self):
        """The donor's own clavicle and scapula tensors are NOT (the second
        moment has a negative eigenvalue), which is why they were replaced.
        Check that what was written instead is."""
        for name, record in body_records(MODEL, set(NEW_BODIES)).items():
            tensor = np.asarray(record['inertia_kg_m2'], float)
            second = np.trace(tensor) / 2 * np.eye(3) - tensor
            self.assertGreater(np.linalg.eigvalsh(tensor).min(), 0, name)
            self.assertGreater(np.linalg.eigvalsh(second).min(), -1e-12, name)

    def test_the_donor_tensors_really_are_inadmissible(self):
        """The control for the one above: it would pass trivially if the donor's
        tensors had been fine all along."""
        donor = ET.parse(ROOT / RECORD['donor']['path']).getroot().find('Model')
        broken = []
        for body in donor.find('BodySet/objects'):
            tensor = inertia_matrix(numbers(body, 'inertia'))
            second = np.trace(tensor) / 2 * np.eye(3) - tensor
            if np.linalg.eigvalsh(second).min() < -1e-12:
                broken.append(body.get('name'))
        self.assertEqual(sorted(broken), ['clavicle', 'radius', 'scapula'])


class File(unittest.TestCase):
    """What the model file declares about the twelve new coordinates."""

    def setUp(self):
        self.declared = coordinates(MODEL)
        self.base = coordinates(ROOT / BASE_MODEL)

    def test_every_new_coordinate_declares_a_range_and_rests_inside_it(self):
        self.assertEqual(len(NEW_RANGES), 12)
        for name, (lower, upper) in NEW_RANGES.items():
            bounds, default = self.declared[name]
            self.assertEqual(bounds, (lower, upper), name)
            self.assertGreater(default, lower, name)
            self.assertLess(default, upper, name)
            self.assertAlmostEqual(upper - lower, 2 * GIRDLE_RANGE_HALFWIDTH_RAD, 12, name)

    def test_every_new_coordinate_carries_a_stop_at_that_range(self):
        stops = {}
        for element in ET.parse(MODEL).getroot().iter('CoordinateLimitForce'):
            stops[element.findtext('coordinate').strip()] = element
        for name, (lower, upper) in NEW_RANGES.items():
            self.assertIn(name, stops, name)
            stop = stops[name]
            self.assertAlmostEqual(float(stop.findtext('lower_limit')) / DEGREES, lower, 12)
            self.assertAlmostEqual(float(stop.findtext('upper_limit')) / DEGREES, upper, 12)
            for side in ('lower_stiffness', 'upper_stiffness'):
                self.assertAlmostEqual(float(stop.findtext(side)) * DEGREES,
                                       STOP_STIFFNESS_NM_PER_RAD, 9)

    def test_no_coordinate_that_already_existed_gained_a_stop(self):
        """The base model's own recorded gates were measured without one; giving
        an old coordinate a stop here would change the plant under them."""
        before = {e.findtext('coordinate').strip()
                  for e in ET.parse(ROOT / BASE_MODEL).getroot().iter('CoordinateLimitForce')}
        after = {e.findtext('coordinate').strip()
                 for e in ET.parse(MODEL).getroot().iter('CoordinateLimitForce')}
        self.assertEqual(after - before, set(NEW_RANGES))
        self.assertEqual(before - after, set())

    def test_the_donor_declares_no_usable_range(self):
        """Why the ranges above are an engineering bound and not a transfer.
        This is the evidence for that sentence, not a restatement of it."""
        donor = ET.parse(ROOT / RECORD['donor']['path']).getroot().find('Model')
        girdle = {'clav_prot', 'clav_elev', 'scapula_abduction', 'scapula_elevation',
                  'scapula_upward_rot', 'scapula_winging'}
        seen = {}
        for element in donor.iter('Coordinate'):
            if element.get('name') in girdle:
                seen[element.get('name')] = tuple(
                    round(float(v), 6) for v in element.findtext('range').split())
        self.assertEqual(len(seen), 6)
        for name, (lower, upper) in seen.items():
            self.assertIn((abs(lower), abs(upper)),
                          [(3.141593, 3.141593), (1.570796, 1.570796)], name)

    def test_the_builder_is_a_function(self):
        files = [ROOT / OUT / n for n in ('model.osim', 'catalog.json', 'registration.json')]
        before = [f.read_bytes() for f in files]
        build(ROOT)
        self.assertEqual([f.read_bytes() for f in files], before)


class Native(unittest.TestCase):
    """The engine's verdict, not the file's."""

    def setUp(self):
        self.stream = None
        self.output = None

    def tearDown(self):
        if self.stream is not None:
            self.stream.close()
        if self.output is not None and (ROOT / self.output).exists():
            shutil.rmtree(ROOT / self.output)

    def open(self, *, environment='free', initial_pose=None, registration=None):
        self.output = 'data/derived/verify-shoulder-girdle-' + uuid.uuid4().hex[:10]
        self.stream = NativeMechanicalStream(
            ROOT, ROOT / self.output, environment=environment,
            target_mass_kg=TARGET_MASS_KG,
            augmented_registration=registration or REGISTRATION,
            initial_pose=initial_pose)
        return self.stream

    def drop(self):
        self.stream.close()
        self.stream = None
        shutil.rmtree(ROOT / self.output)
        self.output = None

    #: A station off every rotation axis and off every body origin.
    STATION_M = (0.05, 0.05, 0.05)

    def stations(self, stream, bodies, station=None):
        station = self.STATION_M if station is None else station
        return {b: np.asarray(stream.body_point(body=b, station_m=station)
                              ['point_source_m'], float) for b in bodies}

    def test_the_model_loads_and_the_engine_agrees_on_its_inventory(self):
        state = self.open().snapshot()
        self.assertEqual(len(state['bodies']), 29)
        self.assertEqual(len(state['muscles']), 158)
        self.assertEqual(len(state['coordinates']), 60)
        self.assertLess(abs(state['mass_kg'] - TARGET_MASS_KG), 1e-9)
        for name, (lower, upper) in NEW_RANGES.items():
            self.assertEqual(state['coordinates'][name]['unit'], 'rad', name)
            value = state['coordinates'][name]['value']
            self.assertGreater(value, lower, name)
            self.assertLess(value, upper, name)

    def test_the_arm_did_not_move_when_it_was_re_parented(self):
        """THE check on the map. acromial_{r,l} now hangs off the scapula instead
        of the torso, through a joint whose parent frame carries the inverse of
        the scapula's reference rotation. If any factor of the map were wrong the
        arm would be somewhere else, and every existing arm result with it.

        Measured against the BASE plant in the same protocol, not against
        arithmetic: four stations per arm body, so orientation is tested and not
        only position."""
        bodies = ('humerus_r', 'ulna_r', 'radius_r', 'hand_r',
                  'humerus_l', 'ulna_l', 'radius_l', 'hand_l', 'torso', 'pelvis')
        probes = ((0, 0, 0), (0.1, 0, 0), (0, 0.1, 0), (0, 0, 0.1))
        stream = self.open()
        variant = {p: self.stations(stream, bodies, p) for p in probes}
        self.drop()
        stream = self.open(registration='data/models/articulated_spine_v1/registration.json')
        base = {p: self.stations(stream, bodies, p) for p in probes}
        self.drop()
        worst = max(float(np.abs(variant[p][b] - base[p][b]).max())
                    for p in probes for b in bodies)
        self.assertLess(worst, 1e-9, 'the arm moved by %.3e m' % worst)

    def test_the_left_girdle_is_the_exact_mirror_of_the_right(self):
        """docs/UPPER_BODY_ACTUATION.md section 8.5: a reflection prior is not
        licensed by a name suffix. The donor has one arm, so everything on the
        left is the right reflected in z, and this is the measurement that says
        whether the reflection came out right -- position AND orientation, from
        four stations per body."""
        pairs = (('clavicle_r', 'clavicle_l'), ('scapula_r', 'scapula_l'),
                 ('humerus_r', 'humerus_l'), ('hand_r', 'hand_l'))
        # Each probe and its z-mirror, so a body mirrored in POSITION but not in
        # ORIENTATION fails: the left body's station at (x, y, -z) must be the
        # z-mirror of the right body's station at (x, y, z).
        probes = ((0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.0, 0.1, 0.0),
                  (0.0, 0.0, 0.1), (0.0, 0.0, -0.1))
        bodies = [b for pair in pairs for b in pair]
        stream = self.open()
        measured = {p: self.stations(stream, bodies, p) for p in probes}
        self.drop()
        worst = 0.0
        for right, left in pairs:
            for probe in probes:
                mirrored = (probe[0], probe[1], -probe[2])
                expected = MIRROR @ measured[probe][right]
                worst = max(worst, float(np.abs(measured[mirrored][left] - expected).max()))
        self.assertLess(worst, 1e-9, 'left/right mirror residual %.3e m' % worst)

    def test_a_single_girdle_coordinate_cannot_be_set_on_its_own(self):
        """Not a defect -- the closed loop, and the engine's guard working.

        Each girdle is `torso -> clavicle` (2 coordinates) plus
        `torso -> scapula` (4) closed by a 3-equation acromioclavicular
        PointConstraint, so only 3 of the 6 are free and no one of them can be
        given a value the other five must then respect. The engine refuses with
        `assembly changed requested initial pose` rather than silently returning
        a different pose, which is exactly what it should do. Recorded here so
        that a caller who tries it finds the reason next to the refusal."""
        for name in ('clav_elev_r', 'scapula_upward_rot_r', 'scapula_abduction_l'):
            lower, upper = NEW_RANGES[name]
            angle = 0.5 * (lower + upper) + 0.4 * (upper - lower) / 2
            with self.assertRaises(RuntimeError, msg=name):
                self.open(initial_pose={name: angle})
            if self.stream is not None:
                self.stream = None
            if self.output is not None and (ROOT / self.output).exists():
                shutil.rmtree(ROOT / self.output)
            self.output = None

    #: joint path -> coordinate, for the sampler's `/value` column labels.
    JOINT_OF = {'clav_elev_r': 'sternoclavicular_r',
                'scapula_upward_rot_r': 'scapulothoracic_r',
                'scapula_abduction_l': 'scapulothoracic_l'}

    def test_each_new_joint_moves_only_its_distal_entities(self):
        """Through OpenSim's own assembler, because of the test above.

        `initial_pose` cannot be used here: the loop means one coordinate cannot
        be set alone. So the perturbation goes through the path-fitting
        sampler, which sets the coordinate and then calls `model.assemble`,
        letting the other girdle coordinates take whatever values the constraint
        needs. Bodies are measured in the TORSO frame, by trilateration from
        four torso stations, so `torso` is the reference and every other body is
        a candidate.

        The half that can fail: `hand` and `humerus` on the driven side MUST
        move, because the arm now hangs off the scapula. The half that can also
        fail: pelvis, thorax, the contralateral girdle and the contralateral
        arm must not.

        THE BAR, and where it comes from. Measured populations over the three
        perturbations, maximum station displacement in the torso frame:

            distal    3.10e-03 m (clavicle_r) to 5.54e-02 m (hand_r)
            proximal  0 to 3.44e-06 m -- and all but one are below 3.2e-11 m.
                      The one is the CONTRALATERAL HAND, 2.47e-06 and 3.44e-06 m,
                      while the contralateral humerus moves 3e-16: the assembler
                      re-solves the whole body (two knee CoordinateCouplerConstraints
                      and two acromioclavicular PointConstraints) and nudges the
                      far wrist.

        The two populations are 900x apart and the bar sits between them. A bar
        at 1e-6, which is what this test was first written with, is BELOW the
        assembly noise and fails on it -- it did, and it read as a topology
        error. Same shape as the bar docs/ARTICULATED_SPINE.md had to move for
        the same reason, and it is an instrument floor, not a pre-registered
        threshold."""
        from scripts.measure_shoulder_donor_frames import (
            ANCHORS, AXIS_LENGTH, STATIONS, probe_model, read_lengths, trilaterate)
        import subprocess
        import tempfile

        watched = ('clavicle_r', 'scapula_r', 'humerus_r', 'hand_r',
                   'clavicle_l', 'scapula_l', 'humerus_l', 'hand_l',
                   'pelvis', 'thorax', 'head', 'femur_r')
        work = Path(tempfile.mkdtemp(prefix='girdle-distal-', dir=ROOT / 'data/derived'))
        try:
            probe = work / 'probe.osim'
            probe_model(ROOT, probe, str(MODEL), 'torso', watched)

            def poses(columns):
                sto = work / 'pose.sto'
                labels = ['time'] + list(columns)
                values = ['0'] + [repr(v) for v in columns.values()]
                sto.write_text('pose\nversion=1\nnRows=1\nnColumns=%d\n'
                               'inDegrees=no\nendheader\n%s\n%s\n'
                               % (len(labels), '\t'.join(labels), '\t'.join(values)))
                out = work / 'lengths.csv'
                subprocess.run([str(ROOT / 'data/runtime/opensim/native_polynomial_path_fit'),
                                'sample', '--model', str(probe), '--coordinates', str(sto),
                                '--output', str(out), '--moment-arm-coordinates',
                                '/jointset/back/lumbar_extension'],
                               cwd=ROOT, check=True, capture_output=True)
                lengths = read_lengths(out)
                result = {}
                for body in watched:
                    points = [trilaterate(ANCHORS,
                                          [lengths['ruler_%s_%d_%d' % (body, a, s)]
                                           for a in range(len(ANCHORS))])[0]
                              for s in range(len(STATIONS))]
                    result[body] = np.concatenate(points)
                return result

            rest = poses({})
            report = {}
            for name, joint in self.JOINT_OF.items():
                lower, upper = NEW_RANGES[name]
                angle = 0.5 * (lower + upper) + 0.4 * (upper - lower) / 2
                moved = poses({'/jointset/%s/%s/value' % (joint, name): angle})
                displacement = {b: float(np.abs(moved[b] - rest[b]).max()) for b in watched}
                report[name] = displacement
                side = name[-1]
                distal = tuple('%s_%s' % (b, side) for b in
                               ('clavicle', 'scapula', 'humerus', 'hand'))
                other = 'l' if side == 'r' else 'r'
                proximal = ('pelvis', 'thorax', 'head', 'femur_r') + tuple(
                    '%s_%s' % (b, other) for b in ('scapula', 'humerus', 'hand'))
                for body in distal:
                    self.assertGreater(displacement[body], 1e-3,
                                       '%s should move %s (%.3e m)'
                                       % (name, body, displacement[body]))
                for body in proximal:
                    self.assertLess(displacement[body], 1e-5,
                                    '%s must not move %s (%.3e m)'
                                    % (name, body, displacement[body]))
            (ROOT / OUT / 'distal_only_report.json').write_text(
                json.dumps({'protocol': 'one girdle coordinate set through OpenSim\'s own '
                                        'assembler; maximum station displacement in the '
                                        'torso frame, four stations per body',
                            'displacement_m': report}, indent=2) + '\n')
        finally:
            shutil.rmtree(work, ignore_errors=True)

    def test_moment_arms_of_the_new_muscles_about_the_new_coordinates(self):
        """The point of the exercise, with this repository's two standing
        controls: the query must be a FUNCTION (call it twice, compare), and a
        coordinate whose answer is known must still print it (soleus_r about
        ankle_angle_r, -0.0497 m).

        The half that can fail on the girdle itself: a scapulohumeral muscle
        must read ZERO about every girdle coordinate, because it spans scapula
        to humerus and crosses no girdle joint; and a thoracoscapular muscle
        must read nonzero, because it spans torso to scapula and crosses one.
        A build that put the muscles on the wrong bodies would break one of
        those two, and a build that put them all on one body would break both."""
        stream = self.open()
        girdle = sorted(NEW_RANGES)
        thoracoscapular = ['seth_SerratusAnterior_M_r', 'seth_TrapeziusScapula_S_r',
                           'seth_Rhomboideus_I_r', 'seth_PectoralisMinor_r',
                           'seth_LevatorScapulae_r']
        scapulohumeral = ['seth_DeltoideusScapula_M_r', 'seth_Supraspinatus_A_r',
                          'seth_Subscapularis_M_r', 'seth_TeresMinor_r']
        muscles = thoracoscapular + scapulohumeral + ['soleus_r']
        query = girdle + ['ankle_angle_r', 'arm_flex_r', 'arm_add_r']
        first = stream.moment_arms(muscles=muscles, coordinates=query)['moment_arms_m']
        second = stream.moment_arms(muscles=muscles, coordinates=query)['moment_arms_m']
        self.assertEqual(first, second, 'moment_arms is not a function of its input')
        self.assertAlmostEqual(first['soleus_r']['ankle_angle_r'], -0.0497, 4)
        for muscle in muscles:
            for name in girdle:
                if muscle in thoracoscapular:
                    continue
                self.assertLess(abs(first[muscle][name]), 1e-9,
                                '%s must not cross %s' % (muscle, name))
        for muscle in thoracoscapular:
            crossed = [abs(first[muscle][n]) for n in girdle if n.endswith('_r')]
            self.assertGreater(max(crossed), 5e-3, muscle)
            for name in girdle:
                if name.endswith('_l'):
                    self.assertLess(abs(first[muscle][name]), 1e-9,
                                    '%s must not cross the LEFT %s' % (muscle, name))
        for muscle in scapulohumeral:
            self.assertGreater(max(abs(first[muscle]['arm_flex_r']),
                                   abs(first[muscle]['arm_add_r'])), 5e-3, muscle)
        self.drop()


if __name__ == '__main__':
    unittest.main()
