"""Gates for data/models/thoracic_drive_v1.

    OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
        .venv/bin/python -m unittest scripts.verify_thoracic_drive -v

ONE engine session at a time; every test opens its own and closes it.

THE TWO HALVES.

  What must NOT have changed.  The build moves ten attachments and two wrap
  surfaces from `torso` to `thorax` by subtracting the thoracic joint offset,
  and the thoracic joint rests at zero.  So at the rest pose the model is the
  SAME MECHANICAL OBJECT as `shoulder_girdle_v1`: every path length, every
  moment arm about every coordinate that is not thoracic, and every body's mass
  and inertia must be unchanged.  That is the known answer, and it is the half
  that catches a wrong sign, a wrong frame or a dropped component.

  What must have changed.  The ten muscles must now read a NONZERO moment arm
  about `thoracic_*`, where `docs/ARTICULATED_SPINE.md` measured every one of
  the plant's 148 muscles below 1e-9 m.

  AND THE CONTROL THAT CAN FAIL.  The second half is worthless on its own,
  because "the new model reads nonzero" is also what you would see if the query
  were broken.  So the same test re-measures the BASE model in the same run and
  requires it to still read zero.  A test that would pass if the base had been
  fine all along is not a test -- IHM-1 CLAUDE.md, "a control that cannot fail
  has not passed".

THIS REPOSITORY'S STANDING CONTROLS, in `test_moment_arms_*`:
  * `moment_arms` called twice at the identical state must compare bit-equal --
    a native query that is not a function of its arguments produces confident
    wrong conclusions.
  * `soleus_r` about `ankle_angle_r` must still read -0.0497080 m.
  * `gait2392_ercspn/intobl/extobl` about `lumbar_extension` must still
    reproduce 42.69 / -52.82 / -62.79 mm.
"""
from __future__ import annotations

import json
import math
import pathlib
import shutil
import subprocess
import sys
import unittest
import uuid
import xml.etree.ElementTree as ET

import numpy as np

from ihm.native.mechanical_stream import NativeMechanicalStream

ROOT = pathlib.Path(__file__).resolve().parents[1]

BASE = 'data/models/shoulder_girdle_v1/registration.json'
VARIANT = 'data/models/thoracic_drive_v1/registration.json'
BASE_MODEL = ROOT / 'data/models/shoulder_girdle_v1/model.osim'
VARIANT_MODEL = ROOT / 'data/models/thoracic_drive_v1/model.osim'

TARGET_MASS_KG = 77.6122029

THORACIC = ('thoracic_extension', 'thoracic_bending', 'thoracic_rotation')

#: The ten the pre-registered measurement selected. Read from the report that
#: travels with the model, never typed here.
REPORT = json.loads(
    (ROOT / 'data/models/thoracic_drive_v1/thoracic_reassignability.json').read_text())
REASSIGNED = sorted(REPORT['reassignable'])

#: Distal of the thoracic joint, and everything else.
DISTAL = ('thorax', 'cervical', 'head')
NOT_DISTAL = ('pelvis', 'torso', 'femur_r', 'tibia_r', 'calcn_r',
              'clavicle_r', 'scapula_r', 'humerus_r', 'radius_r', 'hand_r',
              'clavicle_l', 'scapula_l', 'humerus_l', 'hand_l')

STATION_M = (0.05, 0.05, 0.05)


def elements(model: pathlib.Path, tag: str):
    return {e.get('name'): e for e in ET.parse(model).getroot().iter(tag)}


def triple(text):
    return [float(x) for x in text.split()]


class Files(unittest.TestCase):
    """What the file says, before any engine runs."""

    def test_no_body_mass_mass_centre_or_inertia_changed(self):
        base = elements(BASE_MODEL, 'Body')
        variant = elements(VARIANT_MODEL, 'Body')
        self.assertEqual(sorted(base), sorted(variant))
        for name in base:
            for field in ('mass', 'mass_center', 'inertia'):
                self.assertEqual(base[name].findtext(field),
                                 variant[name].findtext(field),
                                 '%s %s changed' % (name, field))

    def test_the_rib_cage_wrap_did_not_move_in_the_world(self):
        """It changed body AND coordinates; the two must cancel exactly."""
        offset = np.asarray(
            json.loads((ROOT / 'data/models/thoracic_drive_v1/registration.json')
                       .read_text())['changes']['thoracic_joint_offset_m'], float)
        base = elements(BASE_MODEL, 'WrapEllipsoid')
        variant = elements(VARIANT_MODEL, 'WrapEllipsoid')
        for name in ('Thorax_r', 'Thorax_l'):
            before = np.asarray(triple(base[name].findtext('translation')), float)
            after = np.asarray(triple(variant[name].findtext('translation')), float)
            self.assertLess(float(np.abs((after + offset) - before).max()), 1e-15, name)
            self.assertEqual(base[name].findtext('xyz_body_rotation'),
                             variant[name].findtext('xyz_body_rotation'), name)
            self.assertEqual(base[name].findtext('dimensions'),
                             variant[name].findtext('dimensions'), name)

    def test_the_wrap_moved_body_and_the_references_still_resolve(self):
        root = ET.parse(VARIANT_MODEL).getroot()
        homes = {}
        for body in root.iter('Body'):
            for wrap in body.iter():
                if wrap.tag.startswith('Wrap') and wrap.get('name'):
                    homes[wrap.get('name')] = body.get('name')
        self.assertEqual(homes.get('Thorax_r'), 'thorax')
        self.assertEqual(homes.get('Thorax_l'), 'thorax')
        referenced = {w.findtext('wrap_object') for w in root.iter('PathWrap')}
        self.assertTrue(referenced <= set(homes),
                        'a PathWrap names a wrap object no body owns: %r'
                        % sorted(referenced - set(homes)))

    def test_exactly_the_ten_measured_muscles_changed_body(self):
        def attachments(model):
            out = {}
            for muscle in ET.parse(model).getroot().iter():
                if not (muscle.tag.endswith('Muscle') and muscle.get('name')):
                    continue
                out[muscle.get('name')] = sorted(
                    p.findtext('socket_parent_frame') for p in muscle.iter('PathPoint'))
            return out
        base, variant = attachments(BASE_MODEL), attachments(VARIANT_MODEL)
        self.assertEqual(sorted(base), sorted(variant))
        changed = sorted(n for n in base if base[n] != variant[n])
        self.assertEqual(changed, REASSIGNED)

    def test_every_thoracic_coordinate_declares_a_range_and_a_stop_at_it(self):
        """docs/NATIVE_JOINT_LIMITS.md: nothing else enforces a declared range.
        The stop is written in OpenSim's DEGREE convention and the range in
        radians, so the comparison converts."""
        root = ET.parse(VARIANT_MODEL).getroot()
        ranges = {c.get('name'): triple(c.findtext('range'))
                  for c in root.iter('Coordinate') if c.get('name') in THORACIC}
        defaults = {c.get('name'): float(c.findtext('default_value'))
                    for c in root.iter('Coordinate') if c.get('name') in THORACIC}
        stops = {}
        for force in root.iter():
            if force.tag.endswith('CoordinateLimitForce') \
                    and force.findtext('coordinate') in THORACIC:
                stops[force.findtext('coordinate')] = (
                    float(force.findtext('lower_limit')),
                    float(force.findtext('upper_limit')))
        self.assertEqual(sorted(ranges), sorted(THORACIC))
        self.assertEqual(sorted(stops), sorted(THORACIC))
        for name in THORACIC:
            lower, upper = ranges[name]
            self.assertAlmostEqual(stops[name][0], math.degrees(lower), places=9, msg=name)
            self.assertAlmostEqual(stops[name][1], math.degrees(upper), places=9, msg=name)
            self.assertGreater(defaults[name], lower, name)
            self.assertLess(defaults[name], upper, name)

    def test_the_builder_is_a_function(self):
        result = subprocess.run(
            [sys.executable, '-m', 'scripts.build_thoracic_drive', '--check'],
            cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


class Engine(unittest.TestCase):
    """The engine's verdict, not the file's. One session at a time."""

    def setUp(self):
        self.stream = None
        self.output = None

    def tearDown(self):
        self.drop()

    def open(self, registration=VARIANT, initial_pose=None):
        self.output = 'data/derived/verify-thoracic-drive-' + uuid.uuid4().hex[:10]
        self.stream = NativeMechanicalStream(
            ROOT, ROOT / self.output, environment='free',
            target_mass_kg=TARGET_MASS_KG,
            augmented_registration=registration, initial_pose=initial_pose)
        return self.stream

    def drop(self):
        if self.stream is not None:
            self.stream.close()
            self.stream = None
        if self.output is not None:
            shutil.rmtree(ROOT / self.output, ignore_errors=True)
            self.output = None

    def stations(self, stream, bodies):
        return {b: np.asarray(stream.body_point(body=b, station_m=list(STATION_M))
                              ['point_source_m'], float) for b in bodies}

    def test_the_model_loads_and_the_engine_agrees_on_its_inventory(self):
        state = self.open().snapshot()
        self.assertEqual(len(state['bodies']), 29)
        self.assertEqual(len(state['muscles']), 158)
        self.assertEqual(len(state['coordinates']), 60)
        self.assertLess(abs(state['mass_kg'] - TARGET_MASS_KG), 1e-9)

    def test_every_path_length_is_unchanged_at_the_rest_pose(self):
        """THE known answer. The reassignment is the thoracic joint offset
        subtracted exactly and the joint rests at zero, so no path may move."""
        variant = {k: v['path_length_m']
                   for k, v in self.open().snapshot()['muscles'].items()}
        self.drop()
        base = {k: v['path_length_m']
                for k, v in self.open(BASE).snapshot()['muscles'].items()}
        self.drop()
        self.assertEqual(sorted(base), sorted(variant))
        worst, name = 0.0, None
        for key in base:
            delta = abs(base[key] - variant[key])
            if delta > worst:
                worst, name = delta, key
        self.assertLess(worst, 1e-12,
                        'worst path-length change %.3e m at %s' % (worst, name))
        self.assertGreater(len(base), 0)

    def displacement(self, registration):
        rest = self.open(registration)
        before = self.stations(rest, DISTAL + NOT_DISTAL)
        self.drop()
        driven = self.open(registration, initial_pose={'thoracic_extension': 0.2})
        after = self.stations(driven, DISTAL + NOT_DISTAL)
        self.drop()
        return {b: float(np.linalg.norm(after[b] - before[b]))
                for b in DISTAL + NOT_DISTAL}

    @unittest.expectedFailure
    def test_FAILED_distal_only_unpaired_bar_1e9(self):
        """RECORDED FAILED, 18 September 2026, and NEVER RESCORED.

        Written before the run: setting `thoracic_extension = 0.2` must move the
        three bodies distal of the joint and nothing else, with "nothing else"
        at 1e-9 m.

        It FAILED: `tibia_r` moved **2.949e-07 m**, 295x the bar. The bar stays
        where it was written and this test stays marked as failing.

        The bar was the mistake, not the build, and the mistake was choosing it
        without measuring the instrument's own noise floor.
        docs/SHOULDER_GIRDLE.md section 6 records the identical experience one
        variant earlier -- "the distal-only bar had to be measured, and the first
        one was below the noise ... a bar at 1e-6 sits below that and fails on
        it" -- because this plant's assembler re-solves two knee couplers and two
        acromioclavicular point constraints whenever any coordinate is set, and
        that re-solve nudges far bodies by ~1e-7 m whatever was perturbed.

        The question the gate was FOR is answered instead by
        `test_the_reassignment_adds_no_motion_the_base_plant_does_not_have`,
        which is a different quantity -- a paired difference against the base
        plant under the identical perturbation -- and not this threshold moved.
        """
        moved = self.displacement(VARIANT)
        for body in NOT_DISTAL:
            self.assertLess(moved[body], 1e-9, '%s moved %.3e m' % (body, moved[body]))

    def test_the_reassignment_adds_no_motion_the_base_plant_does_not_have(self):
        """The paired version of the gate above, and the one that measures the
        treatment rather than the fixture.

        Identical perturbation on the base plant and on this one. The base has
        the same joint, the same couplers and the same constraints, and NO
        muscle crossing the thoracic joint, so whatever it does to a distant
        body is the assembler and not the reassignment. Only the PAIRED
        DIFFERENCE is attributable here -- IHM-1 CLAUDE.md, "subtract a run of
        the same fixture with the treatment off"."""
        variant = self.displacement(VARIANT)
        base = self.displacement(BASE)
        paired = {b: abs(variant[b] - base[b]) for b in variant}
        for body in DISTAL:
            self.assertGreater(variant[body], 1e-3,
                               '%s did not move: %r' % (body, variant))
            self.assertGreater(base[body], 1e-3, body)
        worst = max(NOT_DISTAL, key=lambda b: paired[b])
        self.assertLess(paired[worst], 1e-9,
                        'the reassignment moved %s by %.3e m more than the base '
                        'plant does under the same perturbation (variant %.3e, '
                        'base %.3e)' % (worst, paired[worst], variant[worst],
                                        base[worst]))
        (ROOT / 'data/models/thoracic_drive_v1/distal_only_report.json').write_text(
            json.dumps({
                'protocol': 'thoracic_extension set to 0.2 rad through the '
                            'engine\'s own initial_pose; station (0.05,0.05,0.05) '
                            'per body; displacement in the native source ground '
                            'frame; the same perturbation run on the base plant '
                            'and reported as a paired difference',
                'unpaired_bar_1e-9_verdict': 'FAILED, recorded, never rescored; '
                                             'see test_FAILED_distal_only_unpaired_bar_1e9',
                'variant_displacement_m': variant,
                'base_displacement_m': base,
                'paired_difference_m': paired,
            }, indent=2) + '\n')

    def test_moment_arms_the_ten_now_drive_the_thoracic_joint_and_the_base_does_not(self):
        """Both halves in one run, with this repository's standing controls.

        The base half is the one that can fail for the right reason: if the
        query were broken, or if the thoracic arms had always been nonzero, the
        variant reading nonzero would mean nothing."""
        probe = REASSIGNED + ['soleus_r', 'gait2392_ercspn_r', 'gait2392_intobl_r',
                              'gait2392_extobl_r', 'seth_TrapeziusScapula_M_r',
                              'seth_DeltoideusScapula_M_r']
        elsewhere = ['ankle_angle_r', 'lumbar_extension', 'arm_flex_r',
                     'scapula_upward_rot_r', 'elbow_flex_r']
        query = list(THORACIC) + elsewhere

        stream = self.open()
        first = stream.moment_arms(muscles=probe, coordinates=query)['moment_arms_m']
        second = stream.moment_arms(muscles=probe, coordinates=query)['moment_arms_m']
        self.assertEqual(first, second, 'moment_arms is not a function of its input')
        self.drop()

        stream = self.open(BASE)
        base = stream.moment_arms(muscles=probe, coordinates=query)['moment_arms_m']
        self.drop()

        # Standing controls, on the variant.
        self.assertAlmostEqual(first['soleus_r']['ankle_angle_r'], -0.0497080, places=6)
        self.assertAlmostEqual(first['gait2392_ercspn_r']['lumbar_extension'],
                               0.04269, places=4)
        self.assertAlmostEqual(first['gait2392_intobl_r']['lumbar_extension'],
                               -0.05282, places=4)
        self.assertAlmostEqual(first['gait2392_extobl_r']['lumbar_extension'],
                               -0.06279, places=4)

        # THE CONTROL THAT CAN FAIL: in the base, nothing crosses the thoracic
        # joint at all, which is what docs/ARTICULATED_SPINE.md measured.
        for muscle in probe:
            for name in THORACIC:
                self.assertLess(abs(base[muscle][name]), 1e-9,
                                'BASE %s already had a %s arm of %.3e m'
                                % (muscle, name, base[muscle][name]))

        # What must have changed: each of the ten now crosses the thoracic joint.
        for muscle in REASSIGNED:
            crossed = max(abs(first[muscle][n]) for n in THORACIC)
            self.assertGreater(crossed, 1e-3,
                               '%s still does not cross the thoracic joint' % muscle)

        # What must NOT have changed: every arm about every other coordinate.
        worst, where = 0.0, None
        for muscle in probe:
            for name in elsewhere:
                delta = abs(first[muscle][name] - base[muscle][name])
                if delta > worst:
                    worst, where = delta, (muscle, name)
        self.assertLess(worst, 1e-12,
                        'a non-thoracic moment arm moved by %.3e m at %r'
                        % (worst, where))


if __name__ == '__main__':
    unittest.main()
