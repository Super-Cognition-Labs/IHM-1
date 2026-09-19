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

Added 2026-09-18, making the new joints drivable (the F2 gate above is NOT
re-scored; it still runs on the kinematic registration and still fails):

  M1 WHY the arms were zero, measured through OpenSim on the model's OWN
     GeometryPaths against the shipped fitted set at the same poses: subtalar and
     mtp are REPRESENTATION (11 and 4 muscles cross them, the fit reads 0); the
     nine spine and four wrist coordinates are TOPOLOGY (the geometry reads 0
     too).  mtp is in the base plant and has read 0 since the engine existed.
  M2 the installed foot path set keeps the shipped coefficients for every path it
     names, and names none of the 22 muscles that cross subtalar or mtp
  M3 with it, those 22 carry their true arms and nothing acquires an arm about a
     spine or wrist coordinate
  N1 the cervical recipe's frames ARE the variant's joint centres
  N2 every one of the donor's 78 neck muscles is accounted for: 50 transferred,
     10 excluded on the girdle, 18 internal to the lumped cervical body
  N3 KNOWN ANSWER: every transferred path length equals the donor's own, evaluated
     by OpenSim on the donor model, at the reference pose
  N4 MEASURED, NOT GATED: moment-arm signs against the donor's generalized arms
  N5 the transfer is a function: run twice, identical bytes
  D  the muscled plant in the engine: 148 muscles, arms about the neck, head and
     subtalar, the soleus and lumbar controls, bit-equal on a repeat call
  G-S PRE-REGISTERED GATE on the foot-paths plant (see its docstring)

Bounded run:
  cd <repo> && OPENBLAS_NUM_THREADS=1 nice -n 10 \
      .venv/bin/python -m unittest scripts.verify_articulated_spine -v

Native sessions are opened ONE AT A TIME and closed in tearDown; the engine is
prlimit-capped at 4 GB by NativeMechanicalStream itself.
"""
from __future__ import annotations

import json
import math
import re
import shutil
import sys
import tempfile
import time
import unittest
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ihm.assembly.cervical_inertia import combine_bodies  # noqa: E402
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402
from ihm.native import path_refitting as pf  # noqa: E402
from scripts import transfer_neck_muscle_paths as neck  # noqa: E402
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
MUSCLED_MODEL = ROOT / neck.MUSCLED_MODEL
MUSCLED_REGISTRATION = neck.MUSCLED_REGISTRATION
FOOT_REGISTRATION = neck.PATHS_REGISTRATION
FOOT_PATHS = ROOT / neck.PATHS
SHIPPED_PATHS = ROOT / ('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/'
                        'example3DWalking/subject_walk_scaled_FunctionBasedPathSet.xml')
MASI = ROOT / neck.MASI
SAMPLER_ARM_FLOOR_M = 1e-4

#: The 22 muscles whose OWN GeometryPath crosses subtalar (11 per side) -- the 4
#: that also cross mtp are among them.  Measured by M1, not assumed by it: M1
#: re-derives this set from the geometry and must reproduce it.
FOOT_CROSSING = sorted('%s_%s' % (m, side) for side in ('r', 'l') for m in (
    'edl', 'ehl', 'fdl', 'fhl', 'gaslat', 'gasmed', 'perbrev', 'perlong',
    'soleus', 'tibant', 'tibpost'))
MTP_CROSSING = sorted('%s_%s' % (m, side) for side in ('r', 'l')
                      for m in ('edl', 'ehl', 'fdl', 'fhl'))

#: `docs/research/LUMBAR_SHOULDER_MUSCLE_COVERAGE.md`'s independently computed
#: lumbar_extension arms, which docs/UPPER_BODY_ACTUATION.md already reproduced.
LUMBAR_TABLE_M = {'gait2392_ercspn_r': 0.04269, 'gait2392_intobl_r': -0.05282,
                  'gait2392_extobl_r': -0.06279}


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


def joint_paths(model_path):
    """{coordinate name: /jointset/<joint>/<coordinate>} in model order."""
    out = {}
    for joint in ET.parse(model_path).getroot().find('.//JointSet/objects'):
        for element in joint.iter('Coordinate'):
            out[element.get('name')] = '/jointset/%s/%s' % (joint.get('name'), element.get('name'))
    return out


def probe_table(model_path, destination, names, offset=0.2):
    """One row at rest, then +-offset rad on each named coordinate in turn."""
    paths = joint_paths(model_path)
    labels = ['time'] + [paths[n] + '/value' for n in names]
    rows = [[0.0] + [0.0] * len(names)]
    for index in range(len(names)):
        for sign in (-1.0, 1.0):
            row = [float(len(rows))] + [0.0] * len(names)
            row[1 + index] = sign * offset
            rows.append(row)
    pf.write_sto(destination, labels, rows)
    return ','.join(paths[n] for n in names)


def peak_arms(csv_path):
    """{(muscle, coordinate): peak |moment arm|} and {muscle: [lengths]}."""
    arms, lengths = {}, {}
    with open(csv_path) as handle:
        next(handle)
        for line in handle:
            _row, _t, actuator, quantity, column, value = line.rstrip('\n').split(',')
            muscle = actuator.rsplit('/', 1)[-1]
            if quantity == 'length':
                lengths.setdefault(muscle, []).append(float(value))
            else:
                key = (muscle, column.rsplit('/', 1)[-1])
                arms[key] = max(arms.get(key, 0.0), abs(float(value)))
    return arms, lengths


def rest_values(csv_path):
    """{(muscle, quantity, coordinate): value} for row 0 (the rest pose)."""
    out = {}
    with open(csv_path) as handle:
        next(handle)
        for line in handle:
            row, _t, actuator, quantity, column, value = line.rstrip('\n').split(',')
            if row == '0':
                out[(actuator.rsplit('/', 1)[-1], quantity, column.rsplit('/', 1)[-1])] = float(value)
    return out


def crossing(arms, coordinate):
    return sorted(m for (m, c), v in arms.items() if c == coordinate and v > SAMPLER_ARM_FLOOR_M)


class MusclePaths(unittest.TestCase):
    """M1-M3 -- why the arms were zero, and what the installed foot paths change.

    Measured with `data/runtime/opensim/native_polynomial_path_fit sample`, which
    evaluates the model's OWN GeometryPaths (wrap objects included) or, given a
    path set, the fitted functions through OpenSim's own evaluator.  Nothing is
    reimplemented in Python."""

    PROBED = tuple(sorted(NEW_RANGES)) + ('mtp_angle_r', 'mtp_angle_l', 'ankle_angle_r')

    @classmethod
    def setUpClass(cls):
        cls.work = Path(tempfile.mkdtemp(prefix='verify-muscle-paths-',
                                         dir=ROOT / 'data/derived'))
        table = cls.work / 'probe.sto'
        cls.columns = probe_table(MODEL, table, cls.PROBED)
        run = lambda name, pathset=None: pf.sample(  # noqa: E731
            MODEL, table, cls.work / (name + '.csv'), pathset=pathset,
            moment_arm_coordinates=cls.columns.split(','), log=cls.work / (name + '.log'))
        run('truth')
        run('shipped', SHIPPED_PATHS)
        cls.truth, _ = peak_arms(cls.work / 'truth.csv')
        cls.shipped, _ = peak_arms(cls.work / 'shipped.csv')
        cls.installed = None
        if FOOT_PATHS.exists():
            run('installed', FOOT_PATHS)
            cls.installed, _ = peak_arms(cls.work / 'installed.csv')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def test_the_zero_splits_into_representation_and_topology(self):
        """M1.  The two classes need opposite fixes, so the split is measured,
        per coordinate, by two routes that share no code: the arm OpenSim reports
        on the model's own GeometryPaths, and whether a path has points on both
        sides of the joint, read off the XML."""
        topology = tuple(n for n in NEW_RANGES if not n.startswith('subtalar'))
        for name in topology:
            self.assertEqual(crossing(self.truth, name), [], '%s: geometry has an arm' % name)
            self.assertEqual(crossing(self.shipped, name), [], name)
        for side in ('r', 'l'):
            real = crossing(self.truth, 'subtalar_angle_' + side)
            self.assertEqual(real, [m for m in FOOT_CROSSING if m.endswith('_' + side)])
            self.assertEqual(crossing(self.shipped, 'subtalar_angle_' + side), [],
                             'the shipped fit cannot express subtalar')
            real = crossing(self.truth, 'mtp_angle_' + side)
            self.assertEqual(real, [m for m in MTP_CROSSING if m.endswith('_' + side)])
            # mtp is a BASE-plant coordinate: this zero predates the variant.
            self.assertEqual(crossing(self.shipped, 'mtp_angle_' + side), [])
        # Control: a coordinate the fit DID carry reads the same set both ways.
        self.assertEqual(crossing(self.truth, 'ankle_angle_r'),
                         crossing(self.shipped, 'ankle_angle_r'))
        self.assertEqual(len(crossing(self.truth, 'ankle_angle_r')), 11)
        # The XML route: every muscle with a measured subtalar arm has path points
        # distal (calcn/toes) AND proximal to the joint, and no other muscle does.
        root = ET.parse(MODEL).getroot()
        spanning = []
        for force in root.find('.//ForceSet/objects'):
            bodies = {p.findtext('socket_parent_frame').rsplit('/', 1)[-1]
                      for p in force.iter('PathPoint')}
            distal = bodies & {'calcn_r', 'toes_r'}
            if distal and bodies - distal:
                spanning.append(force.get('name'))
        self.assertEqual(sorted(spanning), crossing(self.truth, 'subtalar_angle_r'))

    def test_the_installed_foot_paths_keep_the_shipped_coefficients(self):
        """M2.  The installed set is NOT a refit (two refits were
        pre-registered and both FAILED; scripts/refit_muscle_paths_articulated_spine.py).
        It is the shipped set minus the 22 muscles that cross subtalar or mtp.
        Every path it does name must be the shipped path, coefficient string for
        coefficient string."""
        self.assertTrue(FOOT_PATHS.exists(), 'muscle_paths.xml is not installed')
        read = lambda path: {  # noqa: E731
            e.get('name').rsplit('/', 1)[-1]:
            (e.findtext('coordinate_paths').split(),
             e.findtext('length_function/MultivariatePolynomialFunction/coefficients').split())
            for e in ET.parse(path).getroot().iter('FunctionBasedPath')}
        installed, shipped = read(FOOT_PATHS), read(SHIPPED_PATHS)
        self.assertEqual(sorted(set(shipped) - set(installed)), FOOT_CROSSING)
        self.assertEqual(len(installed), 58)
        for name, value in installed.items():
            self.assertEqual(value, shipped[name], name)

    def test_the_installed_set_gives_the_foot_its_true_arms(self):
        """M3.  With the installed set, the 22 run on their own GeometryPaths,
        so their subtalar and mtp arms must equal the geometry's exactly, and no
        coordinate of the topology class may acquire one.  The first half passes
        by construction and is here to catch a set that names a muscle it should
        not; the second half can fail."""
        self.assertIsNotNone(self.installed, 'muscle_paths.xml is not installed')
        for side in ('r', 'l'):
            for coordinate in ('subtalar_angle_' + side, 'mtp_angle_' + side):
                for muscle in crossing(self.truth, coordinate):
                    self.assertAlmostEqual(self.installed[(muscle, coordinate)],
                                           self.truth[(muscle, coordinate)], places=12,
                                           msg='%s about %s' % (muscle, coordinate))
        for name in NEW_RANGES:
            if not name.startswith('subtalar'):
                self.assertEqual(crossing(self.installed, name), [], name)


class NeckTransfer(unittest.TestCase):
    """N1-N5 -- the donor's own neck muscles, moved by the donor's own recipe."""

    #: Our lumped joints against the donor's two independent cervical joints.
    PAIRS = {'neck_extension': 'pitch2', 'neck_bending': 'roll2', 'neck_rotation': 'yaw2',
             'head_extension': 'pitch1', 'head_bending': 'roll1', 'head_rotation': 'yaw1'}
    DONOR_JOINT = {'pitch2': 'auxt1jnt', 'roll2': 'auxt1jnt', 'yaw2': 'auxt1jnt',
                   'pitch1': 'aux2jnt', 'roll1': 'aux2jnt', 'yaw1': 'aux2jnt'}

    @classmethod
    def setUpClass(cls):
        cls.recipe = json.loads((ROOT / neck.RECIPE).read_text())
        cls.registration = json.loads((ROOT / REGISTRATION).read_text())
        cls.transferred, cls.excluded = neck.plan(ROOT)
        # Both sides evaluated by OpenSim at the reference pose.  The donor file
        # does not load in OpenSim 4 as shipped ('/' in component names), so a
        # copy with ONLY those names changed is evaluated, and that no other byte
        # moved is asserted in N3.
        cls.work = Path(tempfile.mkdtemp(prefix='verify-neck-', dir=ROOT / 'data/derived'))
        text = MASI.read_text()
        cls.renamed = re.sub(r'name="([^"]*)"',
                             lambda m: 'name="%s"' % m.group(1).replace('/', '_'), text)
        cls.original = text
        (cls.work / 'masi.osim').write_text(cls.renamed)
        labels, values = ['time'], [0.0]
        for joint in ET.parse(cls.work / 'masi.osim').getroot().iter():
            if joint.tag.endswith('Joint') and re.fullmatch(r'aux(t1|\d)jnt', joint.get('name') or ''):
                for element in joint.iter('Coordinate'):
                    labels.append('/jointset/%s/%s/value' % (joint.get('name'), element.get('name')))
                    values.append(0.0)
        cls.donor_coordinates = len(labels) - 1
        pf.write_sto(cls.work / 'donor.sto', labels, [values])
        probe_table(MUSCLED_MODEL, cls.work / 'ours.sto', tuple(cls.PAIRS))
        ours = joint_paths(MUSCLED_MODEL)
        pf.sample(cls.work / 'masi.osim', cls.work / 'donor.sto', cls.work / 'donor.csv',
                  moment_arm_coordinates=['/jointset/%s/%s' % (cls.DONOR_JOINT[c], c)
                                          for c in cls.PAIRS.values()],
                  log=cls.work / 'donor.log')
        pf.sample(MUSCLED_MODEL, cls.work / 'ours.sto', cls.work / 'ours.csv',
                  moment_arm_coordinates=[ours[c] for c in cls.PAIRS]
                  + [ours[c] for c in NEW_RANGES if not c.startswith(('neck', 'head'))],
                  log=cls.work / 'ours.log')
        cls.donor_rest = rest_values(cls.work / 'donor.csv')
        cls.ours_rest = rest_values(cls.work / 'ours.csv')

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def test_the_recipe_frames_are_the_variants_joint_centres(self):
        """N1.  The transfer re-expresses donor points relative to the variant's
        joint centres, so the two must be the same points -- not close."""
        centres = self.registration['joint_centres_torso_frame_m']
        for donor, joint in (('cerv7', 'neck'), ('skull', 'atlantooccipital')):
            origin = np.asarray(self.recipe['bodies'][donor]
                                ['body_to_target_torso_reference'])[:3, 3]
            self.assertLess(np.abs(origin - np.asarray(centres[joint])).max(), 1e-12, joint)

    def test_every_donor_muscle_is_accounted_for(self):
        """N2.  78 in the donor; each is transferred or excluded with a reason."""
        donor = [m.get('name') for m in ET.parse(MASI).getroot().iter('Thelen2003Muscle')]
        self.assertEqual(len(donor), 78)
        accounted = [m['donor_muscle'] for m in self.transferred] + \
            [e['donor_muscle'] for e in self.excluded]
        self.assertEqual(sorted(accounted), sorted(donor))
        reasons = {}
        for e in self.excluded:
            reasons.setdefault(e['reason'], []).append(e['donor_muscle'])
        self.assertEqual(len(self.transferred), 50)
        self.assertEqual(len(reasons['girdle']), 10)
        self.assertEqual(len(reasons['internal_to_one_lumped_body']), 18)
        # The minimal set the neck needs, by name: sternocleidomastoid (sternal
        # head), splenius and semispinalis, both sides.
        names = {m['donor_muscle'] for m in self.transferred}
        for stem in ('stern_mast', 'splen_cap_sklthx', 'splen_cap_sklc6',
                     'splen_cerv_c3thx', 'semi_cap_sklthx', 'semi_cap_sklc5',
                     'semi_cerv_c3thx'):
            self.assertIn(stem, names)
            self.assertIn(stem + '_l', names)

    def test_path_lengths_equal_the_donors_at_the_reference_pose(self):
        """N3, the known answer.  One rigid transform per donor body, all from
        one rigid fit, so at the reference pose every transferred path is the
        donor's own path moved rigidly and its length cannot change.  A wrong
        transform direction, a wrong body origin or a wrong body map each moves
        lengths by centimetres."""
        self.assertEqual(re.sub(r'name="[^"]*"', '', self.original),
                         re.sub(r'name="[^"]*"', '', self.renamed))
        self.assertEqual(self.donor_coordinates, 24)
        worst = 0.0
        for record in self.transferred:
            theirs = self.donor_rest[(record['donor_muscle'].replace('/', '_'), 'length', '')]
            ours = self.ours_rest[(record['id'], 'length', '')]
            worst = max(worst, abs(theirs - ours))
        print('\n  worst |transferred - donor| path length at the reference pose: %.3g m' % worst)
        self.assertLess(worst, 1e-9)

    def test_moment_arm_signs_against_the_donor_are_the_measured_ones(self):
        """N4 -- MEASURED, NOT A GATE.  The donor's arm about pitch2 is a
        GENERALIZED arm summed over seven coupled levels; ours about
        neck_extension is one joint at C7/T1.  They are different quantities, so
        no bar on their agreement is derivable, and none is set.  What is pinned
        is the measurement itself, 2026-09-18, at the reference pose, counting
        (muscle, coordinate) pairs where BOTH arms exceed 0.1 mm:

            160 agree in sign, 8 disagree, 28 are the donor's only (an arm about a
            level our lumped cervical body has absorbed), 0 are ours only.

        All 8 disagreements are left/right pairs, i.e. structural, and all sit
        where the donor's summed arm nearly cancels (|donor| <= 1.7 mm) except
        longissimus capitis about head_rotation (4.4 vs 4.8 mm), whose donor
        joint is C2/C1 and ours the occiput.  If this list changes, something
        moved; read it before accepting the change."""
        agree, disagree, donor_only, ours_only = 0, [], 0, 0
        for record in self.transferred:
            name = record['donor_muscle'].replace('/', '_')
            for ours, theirs in self.PAIRS.items():
                a = self.ours_rest[(record['id'], 'moment_arm', ours)]
                b = self.donor_rest[(name, 'moment_arm', theirs)]
                if abs(a) > SAMPLER_ARM_FLOOR_M and abs(b) > SAMPLER_ARM_FLOOR_M:
                    if (a > 0) == (b > 0):
                        agree += 1
                    else:
                        disagree.append((record['id'], ours))
                elif abs(a) > SAMPLER_ARM_FLOOR_M:
                    ours_only += 1
                elif abs(b) > SAMPLER_ARM_FLOOR_M:
                    donor_only += 1
        self.assertEqual((agree, donor_only, ours_only), (160, 28, 0))
        self.assertEqual(sorted(disagree), sorted(
            [('masi_%s_%s' % (m, s), c) for s in ('r', 'l') for m, c in (
                ('scalenus_med', 'neck_extension'), ('scalenus_post', 'neck_extension'),
                ('long_col_c5thx', 'neck_rotation'), ('longissi_cap_sklc6', 'head_rotation'))]))

    def test_no_neck_muscle_reaches_the_thoracic_joint(self):
        """The most proximal body a transferred muscle touches is `thorax`, so
        none may have an arm about thoracic_* -- that stays topology-class."""
        for record in self.transferred:
            for name in ('thoracic_extension', 'thoracic_bending', 'thoracic_rotation'):
                self.assertLess(abs(self.ours_rest[(record['id'], 'moment_arm', name)]),
                                1e-9, '%s about %s' % (record['id'], name))

    def test_the_transfer_is_a_function(self):
        """N5.  Run it twice; identical bytes."""
        files = [ROOT / p for p in (neck.MUSCLED_MODEL, neck.MUSCLED_CATALOG,
                                    neck.MUSCLED_REGISTRATION, neck.PATHS_REGISTRATION)
                 if (ROOT / p).exists()]
        before = [f.read_bytes() for f in files]
        neck.build(ROOT)
        self.assertEqual([f.read_bytes() for f in files], before)


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

    def timed_tonic_excursions(self, registration, model):
        """The same protocol as `tonic_excursions`, also returning wall clock per
        10 ms advance.  A separate method so F2's own code path is untouched."""
        declared = {k: v for k, v in coordinates(model).items() if v is not None}
        stream = self.open(environment='supine', registration=registration)
        excitation = {m: 0.02 for m in stream.snapshot()['muscles']}
        worst = {}
        started = time.time()
        for _ in range(200):
            state = stream.advance(0.01, actuation=excitation)
            for name, bounds in declared.items():
                value = state['coordinates'][name]['value']
                worst[name] = max(worst.get(name, -1e9),
                                  bounds[0] - value, value - bounds[1])
        seconds = (time.time() - started) / 200
        self.drop()
        return worst, seconds

    def test_the_muscled_plant_has_arms_about_the_neck_and_the_foot(self):
        """D.  The engine's verdict on the muscled plant, with the controls
        docs/ARTICULATED_SPINE.md used: soleus_r about the ankle, the lumbar arms
        against their independently computed table, and a repeat call that must
        be bit-equal.  Then what changed: every transferred neck muscle has an
        arm about the neck or head, the 22 foot muscles have one about subtalar,
        and thoracic_* and both wrists stay at zero for all 148 muscles."""
        stream = self.open(registration=MUSCLED_REGISTRATION)
        state = stream.snapshot()
        self.assertEqual(len(state['muscles']), 148)
        self.assertEqual(len(state['bodies']), 25)
        self.assertEqual(len(state['coordinates']), 48)
        self.assertLess(abs(state['mass_kg'] - TARGET_MASS_KG), 1e-9)
        muscles = sorted(state['muscles'])
        query = sorted(NEW_RANGES) + ['ankle_angle_r', 'lumbar_extension', 'mtp_angle_r']
        first = stream.moment_arms(muscles=muscles, coordinates=query)['moment_arms_m']
        second = stream.moment_arms(muscles=muscles, coordinates=query)['moment_arms_m']
        self.assertEqual(first, second, 'moment_arms is not a function of its input')
        soleus = first['soleus_r']['ankle_angle_r']
        print('\n  soleus_r about ankle_angle_r: %+.5f m' % soleus)
        self.assertLess(abs(soleus - (-0.0497)), 5e-4)
        for name, expected in LUMBAR_TABLE_M.items():
            self.assertAlmostEqual(first[name]['lumbar_extension'], expected, delta=1e-5, msg=name)
        neck_coordinates = [c for c in NEW_RANGES if c.startswith(('neck', 'head'))]
        added = [m for m in muscles if m.startswith('masi_')]
        self.assertEqual(len(added), 50)
        for muscle in added:
            self.assertGreater(max(abs(first[muscle][c]) for c in neck_coordinates),
                               SAMPLER_ARM_FLOOR_M, muscle)
        for side in ('r', 'l'):
            moved = sorted(m for m in muscles
                           if abs(first[m].get('subtalar_angle_' + side, 0.0)) > SAMPLER_ARM_FLOOR_M)
            self.assertEqual(moved, [m for m in FOOT_CROSSING if m.endswith('_' + side)])
        moved = sorted(m for m in muscles if abs(first[m]['mtp_angle_r']) > SAMPLER_ARM_FLOOR_M)
        self.assertEqual(moved, [m for m in MTP_CROSSING if m.endswith('_r')])
        for muscle in muscles:
            for name in NEW_RANGES:
                if name.startswith(('thoracic', 'wrist')):
                    self.assertLess(abs(first[muscle][name]), 1e-9, '%s about %s' % (muscle, name))
        # Extensors extend and the sternocleidomastoid flexes the lower neck: the
        # anatomical sign, read from the engine rather than from the donor.
        self.assertGreater(first['masi_splen_cap_sklthx_r']['neck_extension'], 0.01)
        self.assertGreater(first['masi_semi_cap_sklthx_r']['head_extension'], 0.01)
        self.assertLess(first['masi_stern_mast_r']['neck_extension'], -0.01)

    def test_gate_the_foot_paths_repair_the_ankle(self):
        """G-S, PRE-REGISTERED 2026-09-18, written before any plant carrying the
        installed foot paths had been integrated.

        What it answers: docs/ARTICULATED_SPINE.md measured that un-welding the
        subtalar made the ankle leave its range by 1.36 rad instead of 0.12,
        because the plantarflexors loaded a hinge no muscle could control.  With
        the 22 foot muscles on their own GeometryPaths, that hinge has 11
        muscles per side.  Does the regression go away?

        Protocol: docs/NATIVE_JOINT_LIMITS.md's, unchanged -- supine, 0.02 tonic
        on every muscle, 2 s, worst excursion past each declared range.
        Bar: F2's bar, the BASE model's own worst excursion under the identical
        protocol, measured in the same run.  No new constant is introduced.
        Plant under test: registration_foot_paths.json -- model.osim UNCHANGED
        plus only the path-set override, so a pass or a fail is attributable to
        the foot paths and nothing else.
        Coordinates gated: ankle_angle_{l,r} and subtalar_angle_{l,r}, the four
        the intervention acts on.  Everything else is printed, not gated,
        including thoracic_extension, which no muscle in any of these plants
        crosses and on which F2 already failed.
        The muscled plant (neck muscles added too) is measured and printed in the
        same run, not gated.
        If it fails: recorded FAILED; the threshold does not move and nothing is
        tuned.

        Measured 2026-09-18, same run as F2:
            bar      0.2202 rad (base knee_angle_r)
            ankle_angle_r    base 0.1229   kinematic variant 1.3563   foot paths 1.1399
            ankle_angle_l    base 0.1207   kinematic variant 1.3522   foot paths 0.8709
            subtalar_angle_r                                            foot paths 0.1865
            subtalar_angle_l                                            foot paths 0.1889
            wall clock per advance: base 0.291 s, foot paths 0.405 s, muscled 0.412 s
        VERDICT: FAIL on both ankles, 5.2x and 4.0x the bar.  The subtalars stay
        inside it.  The arms are real (M3, D) and giving them to the muscles
        removes 16% and 36% of the ankle excursion, not the regression: the
        cause of the ankle collapse is NOT only that subtalar had no arm, and
        what else it is was not measured here.  The separating control not yet
        run: the kinematic variant with ONLY the subtalar re-welded.

        RUN 2026-09-18 (docs/FOOT_JOINTS.md, pre-registered in e00844d): it does
        NOT recover (1.3847 / 1.2665), subtalar freedom alone on the base trunk
        does not collapse it (0.1129 / 0.1080), and the variant with all seven
        new joints welded still collapses (1.4014 / 1.2876).  The cause is the
        torso partition's mass/contact redistribution, not the subtalar; the
        supine plane 48.5 mm lower is the leading candidate, not separated.  This gate's verdict is unchanged; its premise is withdrawn."""
        base, base_s = self.timed_tonic_excursions(BASE + '/registration.json',
                                                   ROOT / BASE / 'model.osim')
        foot, foot_s = self.timed_tonic_excursions(FOOT_REGISTRATION, MODEL)
        muscled, muscled_s = self.timed_tonic_excursions(MUSCLED_REGISTRATION, MUSCLED_MODEL)
        bar = max(base.values())
        gated = ('ankle_angle_r', 'ankle_angle_l', 'subtalar_angle_r', 'subtalar_angle_l')
        print('\n  bar (base worst): %.4f rad (%s)' % (bar, max(base, key=base.get)))
        print('  wall clock per 10 ms advance: base %.3f s, foot paths %.3f s, muscled %.3f s'
              % (base_s, foot_s, muscled_s))
        for name in gated + ('thoracic_extension', 'neck_extension', 'head_extension'):
            print('  %-20s base %s  foot %.4f  muscled %.4f' % (
                name, ('%.4f' % base[name]) if name in base else '   -  ',
                foot[name], muscled[name]))
        worse = {k: round(foot[k] - base[k], 4) for k in base if foot[k] - base[k] > 0.05}
        print('  EXISTING coordinates the foot-path plant made worse by >0.05 rad: %r' % worse)
        over = {k: round(foot[k], 4) for k in gated if foot[k] > bar}
        print('  GATED coordinates over the bar: %r' % over)
        self.assertEqual(over, {})


if __name__ == '__main__':
    unittest.main(verbosity=2)
