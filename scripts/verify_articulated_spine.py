"""Known-answer checks for data/models/articulated_spine_v1.

Every test here has an answer that is known before the run, and at least one of
them can FAIL for the reason it exists (IHM-1 CLAUDE.md, "a control that cannot
fail has not passed"): the distal-only test asserts both that the entities below
a new joint MOVE and that the entities above it do not.  A joint that did nothing
at all would pass the second half alone.

  A  total mass is conserved to 1e-9 against engineering_stance_v1, and the
     four-body partition recombines to that model's ONE torso body exactly
  B  every new rotational coordinate declares a <range> that contains zero
  C  every new rotational coordinate carries a CoordinateLimitForce whose limits
     ARE the declared range, and OpenSim's degree convention is honoured
  D  the model loads in the native engine at 25 bodies / 98 muscles / 48
     coordinates and the engine's own total mass is the target mass
  E  each new joint moves its distal entities and nothing proximal
  F1 the rest pose is strictly inside every new declared range, so no new stop
     is active on a body that has not moved
  F2 PRE-REGISTERED GATE, AND IT FAILS: under docs/NATIVE_JOINT_LIMITS.md's own
     protocol no new coordinate may leave its range by more than the base model
     already leaves its own.  thoracic_extension reaches 0.394 rad against a
     0.220 rad bar.  Recorded as FAILED; the threshold is not moved.
  H  no muscle has a moment arm about any new coordinate -- a known answer about
     what this variant IS, with a control (soleus_r about ankle_angle_r) that
     comes back nonzero
  G  the builder is a function: run twice, identical bytes

Bounded run:
  cd <repo> && OPENBLAS_NUM_THREADS=1 nice -n 10 \
      .venv/bin/python -m unittest scripts.verify_articulated_spine -v

Native sessions are opened ONE AT A TIME and closed in tearDown; the engine is
prlimit-capped at 4 GB by NativeMechanicalStream itself.
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
from scripts.build_articulated_spine import (  # noqa: E402
    BASE, CERVICAL_INERTIA, HEAD_RANGE, MASS_SCALE, NECK_RANGE, OUT,
    STOP_STIFFNESS_NM_PER_RAD, SUBTALAR_RANGE, THORACIC_RANGE, WRIST_RANGE,
    build, read_partition,
)

TARGET_MASS_KG = 77.6122029
DEGREES = 180.0 / np.pi

#: Every coordinate this variant adds, and the range it must declare.
NEW_RANGES = {}
NEW_RANGES.update(THORACIC_RANGE)
NEW_RANGES.update(NECK_RANGE)
NEW_RANGES.update(HEAD_RANGE)
for _side in ('r', 'l'):
    NEW_RANGES['wrist_flex_%s' % _side] = WRIST_RANGE['wrist_flex']
    NEW_RANGES['wrist_dev_%s' % _side] = WRIST_RANGE['wrist_dev']
    NEW_RANGES['subtalar_angle_%s' % _side] = SUBTALAR_RANGE

#: What each new joint is ALLOWED to move, and one body it must not.  The second
#: half is what makes this a control that can fail.
DISTAL = {
    'thoracic_extension': (('thorax', 'cervical', 'head'), ('torso', 'pelvis', 'humerus_r')),
    'neck_extension': (('cervical', 'head'), ('thorax', 'torso', 'humerus_r')),
    'head_rotation': (('head',), ('cervical', 'thorax', 'torso')),
    'wrist_flex_r': (('hand_r',), ('radius_r', 'ulna_r', 'hand_l', 'torso')),
    'subtalar_angle_r': (('calcn_r', 'toes_r'), ('talus_r', 'tibia_r', 'calcn_l')),
}

MODEL = ROOT / OUT / 'model.osim'
REGISTRATION = OUT + '/registration.json'


def coordinates(path):
    root = ET.parse(path).getroot()
    found = {}
    for element in root.iter('Coordinate'):
        text = element.findtext('range')
        found[element.get('name')] = None if text is None else tuple(
            float(v) for v in text.split())
    return found


def masses(path):
    root = ET.parse(path).getroot()
    return {b.get('name'): float(b.findtext('mass'))
            for b in root.findall('.//BodySet/objects/Body')}


class Partition(unittest.TestCase):
    """A -- the torso was repartitioned, never added to."""

    def test_total_mass_conserved(self):
        base = masses(ROOT / BASE / 'model.osim')
        variant = masses(MODEL)
        self.assertEqual(len(base), 22)
        self.assertEqual(len(variant), 25)
        self.assertLess(abs(sum(variant.values()) - sum(base.values())), 1e-9)

    def test_only_torso_changed(self):
        base = masses(ROOT / BASE / 'model.osim')
        variant = masses(MODEL)
        for name, value in base.items():
            if name == 'torso':
                continue
            self.assertEqual(value, variant[name], name)
        debited = base['torso'] - variant['torso']
        added = sum(variant[n] for n in ('thorax', 'cervical', 'head'))
        self.assertLess(abs(debited - added), 1e-9)

    def test_partition_recombines_to_the_original_torso(self):
        """The strongest known answer available: mass, COM and the full inertia
        tensor of the four parts must reassemble the one body they came from."""
        inertia = json.loads((ROOT / CERVICAL_INERTIA).read_text())
        _, _, partition = read_partition(ROOT)
        total = combine_bodies([partition[k] for k in ('head', 'cervical', 'thorax', 'torso')])
        original = inertia['current_torso']
        self.assertLess(abs(total['mass_kg'] - original['mass_kg']), 1e-12)
        self.assertLess(np.abs(np.asarray(total['center_m'])
                               - np.asarray(original['center_m'])).max(), 1e-12)
        self.assertLess(np.linalg.norm(np.asarray(total['inertia_kg_m2'])
                                       - np.asarray(original['inertia_kg_m2'])), 1e-12)

    def test_new_body_inertias_are_physically_admissible(self):
        """The engine inscribes a contact sphere per body from the inertia
        ellipsoid and throws if the triangle inequality fails; check it here, in
        the file, rather than discovering it in a native traceback."""
        root = ET.parse(MODEL).getroot()
        for body in root.findall('.//BodySet/objects/Body'):
            moments = np.asarray([float(v) for v in body.findtext('inertia').split()][:3])
            mass = float(body.findtext('mass'))
            for k in range(3):
                slack = moments[(k + 1) % 3] + moments[(k + 2) % 3] - moments[k]
                self.assertGreater(5 * slack / (2 * mass), 0,
                                   '%s axis %d' % (body.get('name'), k))


class Ranges(unittest.TestCase):
    """B and C -- a range, and a stop AT that range, on every new coordinate."""

    def setUp(self):
        self.declared = coordinates(MODEL)
        self.stops = {}
        for force in ET.parse(MODEL).getroot().iter('CoordinateLimitForce'):
            self.stops[force.findtext('coordinate')] = force

    def test_every_new_coordinate_is_present_and_declares_its_range(self):
        base = coordinates(ROOT / BASE / 'model.osim')
        self.assertEqual(set(self.declared) - set(base), set(NEW_RANGES))
        for name, bounds in NEW_RANGES.items():
            self.assertIsNotNone(self.declared[name], name)
            self.assertAlmostEqual(self.declared[name][0], bounds[0], places=12, msg=name)
            self.assertAlmostEqual(self.declared[name][1], bounds[1], places=12, msg=name)

    def test_zero_is_inside_every_new_range(self):
        for name, bounds in NEW_RANGES.items():
            self.assertLess(bounds[0], 0.0, name)
            self.assertGreater(bounds[1], 0.0, name)

    def test_every_new_coordinate_has_a_stop_at_its_declared_range(self):
        self.assertEqual(set(self.stops), set(NEW_RANGES))
        for name, bounds in NEW_RANGES.items():
            force = self.stops[name]
            self.assertAlmostEqual(float(force.findtext('lower_limit')),
                                   bounds[0] * DEGREES, places=9, msg=name)
            self.assertAlmostEqual(float(force.findtext('upper_limit')),
                                   bounds[1] * DEGREES, places=9, msg=name)
            for side in ('lower_stiffness', 'upper_stiffness'):
                self.assertAlmostEqual(float(force.findtext(side)) * DEGREES,
                                       STOP_STIFFNESS_NM_PER_RAD, places=9, msg=name)

    def test_no_existing_coordinate_gained_a_stop(self):
        """A variant that silently stopped the OLD coordinates would be a
        different plant for a reason nobody asked for."""
        base = set(coordinates(ROOT / BASE / 'model.osim'))
        self.assertEqual(set(self.stops) & base, set())

    def test_the_welds_are_gone(self):
        root = ET.parse(MODEL).getroot()
        names = {j.get('name'): j.tag for j in root.find('.//JointSet/objects')}
        self.assertEqual(names['radius_hand_r'], 'UniversalJoint')
        self.assertEqual(names['radius_hand_l'], 'UniversalJoint')
        self.assertEqual(names['subtalar_r'], 'PinJoint')
        self.assertEqual(names['subtalar_l'], 'PinJoint')
        self.assertNotIn('WeldJoint', names.values())


class Idempotent(unittest.TestCase):
    """G -- the builder is a function of its inputs, so it returns the same
    bytes when called twice.  This is not a known answer; it tests whether the
    instrument is a function at all."""

    def test_rebuilding_reproduces_the_same_model(self):
        before = MODEL.read_bytes()
        registration = json.loads((ROOT / OUT / 'registration.json').read_text())
        build(ROOT)
        self.assertEqual(MODEL.read_bytes(), before)
        self.assertEqual(json.loads((ROOT / OUT / 'registration.json').read_text()),
                         registration)


class Native(unittest.TestCase):
    """D, E, F -- the engine's verdict, not the file's."""

    def setUp(self):
        self.stream = None
        self.output = None

    def tearDown(self):
        if self.stream is not None:
            self.stream.close()
        if self.output is not None and (ROOT / self.output).exists():
            shutil.rmtree(ROOT / self.output)

    def open(self, *, environment='free', initial_pose=None, registration=None):
        self.output = 'data/derived/verify-articulated-spine-' + uuid.uuid4().hex[:10]
        self.stream = NativeMechanicalStream(
            ROOT, ROOT / self.output, environment=environment,
            target_mass_kg=TARGET_MASS_KG,
            augmented_registration=registration or REGISTRATION,
            initial_pose=initial_pose)
        return self.stream

    def drop(self):
        """Close and delete the current session.  ONE native session at a time."""
        self.stream.close()
        self.stream = None
        shutil.rmtree(ROOT / self.output)
        self.output = None

    def tonic_excursions(self, registration, model):
        """The docs/NATIVE_JOINT_LIMITS.md protocol: supine, 0.02 tonic on every
        muscle, 2 s, worst excursion past each coordinate's DECLARED range."""
        declared = {k: v for k, v in coordinates(model).items() if v is not None}
        stream = self.open(environment='supine', registration=registration)
        excitation = {m: 0.02 for m in stream.snapshot()['muscles']}
        worst = {}
        for _ in range(200):
            state = stream.advance(0.01, actuation=excitation)
            for name, bounds in declared.items():
                value = state['coordinates'][name]['value']
                worst[name] = max(worst.get(name, -1e9),
                                  bounds[0] - value, value - bounds[1])
        self.drop()
        return worst

    def test_the_model_loads_and_the_engine_agrees_on_its_inventory(self):
        state = self.open().snapshot()
        self.assertEqual(len(state['bodies']), 25)
        self.assertEqual(len(state['muscles']), 98)
        self.assertEqual(len(state['coordinates']), 48)
        self.assertLess(abs(state['mass_kg'] - TARGET_MASS_KG), 1e-9)
        for name in NEW_RANGES:
            self.assertEqual(state['coordinates'][name]['unit'], 'rad', name)
            self.assertEqual(state['coordinates'][name]['value'], 0, name)

    #: A station OFF every rotation axis and off the body origin.  The origin of
    #: thorax, cervical and head IS its proximal joint centre, so probing the
    #: origin would report 8e-12 m for a joint that had just rotated 15 degrees
    #: -- a pivot does not move under its own rotation.  That is exactly the
    #: shape of a control that cannot fail, and it did fail here first.
    STATION_M = (0.05, 0.05, 0.05)

    def origins(self, stream, bodies):
        return {b: np.asarray(stream.body_point(body=b, station_m=self.STATION_M)
                              ['point_source_m'], float) for b in bodies}

    def test_each_new_joint_moves_only_its_distal_entities(self):
        every = sorted({b for pair in DISTAL.values() for group in pair for b in group})
        rest = self.origins(self.open(), every)
        self.drop()
        for name, (distal, proximal) in DISTAL.items():
            angle = 0.5 * min(abs(NEW_RANGES[name][0]), abs(NEW_RANGES[name][1]))
            moved = self.origins(self.open(initial_pose={name: angle}), every)
            displacement = {b: float(np.linalg.norm(moved[b] - rest[b])) for b in every}
            for body in distal:
                # The half that can fail: a joint that did nothing would pass the
                # proximal half and fail here.
                self.assertGreater(displacement[body], 1e-3,
                                   '%s should move %s' % (name, body))
            for body in proximal:
                # NOT zero, and the reason is worth writing down: this model
                # carries CoordinateCouplerConstraints on both knees, so OpenSim
                # re-assembles the WHOLE body whenever any coordinate is set and
                # the solver's answer depends on where it started.  Measured
                # here: 8e-12 m for the three spine joints, and up to 1.46e-6 m
                # for the wrist and subtalar.  The smallest distal displacement
                # in the same sweep is 9.2e-3 m, so the two populations are
                # 6,300x apart and the bar below sits between them.  A bar at
                # 1e-12 would have failed on assembly noise and read as a
                # topology error.
                self.assertLess(displacement[body], 1e-5,
                                '%s must not move %s' % (name, body))
            self.drop()

    def test_the_rest_pose_is_strictly_inside_every_new_declared_range(self):
        """F1, and this is the check docs/NATIVE_JOINT_LIMITS.md actually asks
        for.  That document's finding is that the source model's declared range
        and its own passive stop DISAGREE -- `hip_rotation_l`'s stop sits 0.222
        rad wider than its range, so a body doing nothing already violates the
        range and no guard at the declared range can be used.

        Here the stop IS the declared range by construction, so the equivalent
        failure is a rest pose that is not inside it.  Read the engine's own
        coordinate values at t=0, before any integration."""
        state = self.open(environment='upright').snapshot()
        for name, (lower, upper) in NEW_RANGES.items():
            value = state['coordinates'][name]['value']
            self.assertGreater(value, lower, name)
            self.assertLess(value, upper, name)
            self.assertEqual(value, 0, name)

    def test_no_muscle_has_a_moment_arm_about_any_new_coordinate(self):
        """Not a pass/fail on the model -- a KNOWN ANSWER that pins down what
        this variant is.  The 80 source muscle paths are fitted polynomials in
        the coordinates that existed when they were fitted
        (subject_walk_scaled_FunctionBasedPathSet.xml, substituted in by
        ModelFactory::replacePathsWithFunctionBasedPaths), so un-welding a joint
        does NOT give them a moment arm about it.  The control that can fail is
        the same query at a coordinate that WAS in the fit: soleus_r about
        ankle_angle_r must come back nonzero."""
        stream = self.open()
        muscles = ['soleus_r', 'gasmed_r', 'tibant_r', 'gait2392_ercspn_r',
                   'arm26_BIClong_r']
        arms = stream.moment_arms(muscles=muscles,
                                  coordinates=sorted(NEW_RANGES) + ['ankle_angle_r',
                                                                    'lumbar_extension'])
        self.assertLess(abs(arms['moment_arms_m']['soleus_r']['ankle_angle_r']), 0.06)
        self.assertGreater(abs(arms['moment_arms_m']['soleus_r']['ankle_angle_r']), 0.01)
        self.assertGreater(
            abs(arms['moment_arms_m']['gait2392_ercspn_r']['lumbar_extension']), 0.01)
        for muscle in muscles:
            for name in NEW_RANGES:
                self.assertEqual(arms['moment_arms_m'][muscle][name], 0,
                                 '%s about %s' % (muscle, name))

    def test_gate_new_coordinates_stay_within_the_base_models_worst_excursion(self):
        """F2, PRE-REGISTERED GATE -- AND IT FAILS.  Recorded, not rescored.

        Protocol, taken unchanged from docs/NATIVE_JOINT_LIMITS.md: a body doing
        nothing, 0.02 tonic excitation on every muscle, 2 s, then the worst
        excursion past each declared range.  Bar: no NEW coordinate may exceed
        the worst excursion the BASE model already shows on its OWN coordinates
        under the identical protocol, measured in the same run.

        Measured 2026-09-18:
            base   worst 0.2202 rad (knee_angle_r); knee_angle_l 0.2195,
                   pro_sup 0.1809, ankle_angle 0.1229
            new    thoracic_extension 0.3940, thoracic_bending 0.2099,
                   neck_extension 0.2084, subtalar_angle_l 0.1716,
                   subtalar_angle_r 0.1263, head_extension 0.1261,
                   neck_bending 0.0829; the other eight stay inside
        VERDICT: FAIL on thoracic_extension, 1.79x the bar.

        Two things the number is NOT, and both were checked before writing this:
        the stop is present and it is at the declared range (Ranges, above), and
        the rest pose is inside the range (F1, above).  A CoordinateLimitForce is
        a soft stop -- 30 N*m/rad is the only stiffness this repository has
        measured, and that same sweep reports 0.220 rad of residual excursion for
        it, so 0.394 is in family rather than anomalous.  CONFOUND, stated
        because it is not separable here: the supine environment gives every body
        a posterior contact sphere inscribed in its inertia ellipsoid, so the
        three new bodies bring three new balls and part of this moment is the
        contact proxy rather than the stop.  What would settle it is muscles
        across these joints (docs/WORKBENCH_AUTHENTICITY.md item 0.1), not a
        stiffer number."""
        base = self.tonic_excursions(BASE + '/registration.json',
                                     ROOT / BASE / 'model.osim')
        variant = self.tonic_excursions(REGISTRATION, MODEL)
        # The bar is the BASE MODEL's own worst, measured on the base model.
        # Taking it from the variant's non-new coordinates would be the wrong
        # population: the variant's ankle folds to 1.36 rad past its range once
        # the subtalar weld is gone, which would raise the bar 6x and hand this
        # gate a pass it has not earned.
        bar = max(base.values())
        over = {k: round(variant[k], 4) for k in NEW_RANGES if variant[k] > bar}
        drift = {k: round(variant[k] - base[k], 4) for k in base
                 if variant[k] - base[k] > 0.05}
        print('\n  base worst excursion past a declared range: %.4f rad (%s)'
              % (bar, max(base, key=base.get)))
        print('  new coordinates over that bar: %r' % over)
        print('  EXISTING coordinates the variant made worse by >0.05 rad: %r' % drift)
        self.assertEqual(over, {})


if __name__ == '__main__':
    unittest.main(verbosity=2)
