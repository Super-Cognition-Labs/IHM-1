"""Can the body press its own trunk up? The bound, and what the bound is worth.

    OPENBLAS_NUM_THREADS=1 nice -n 10 \\
        .venv/bin/python -m scripts.measure_shoulder_girdle_capacity

docs/UPPER_BODY_ACTUATION.md section 3 ends: "Pressing the trunk up off the floor
is a demand of a different order and NOTHING HERE BOUNDS IT."  This measures the
bound the girdle makes measurable, and it is worth being precise about why the
girdle is what makes the question answerable at all.

WITHOUT A GIRDLE THE QUESTION HAS NO ANSWER, and not because the plant is strong.
In `articulated_spine_v1` the humerus is jointed straight to `torso`, so at the
top of a press-up -- trunk horizontal, arms straight, hands under the shoulders --
the ground reaction runs hand -> radius -> ulna -> humerus -> glenohumeral joint
-> torso entirely through JOINT REACTIONS.  No muscle is in the load path, the
skeleton carries any force whatever, and "how much can it press" is unbounded for
a reason that says nothing about strength.

WITH the girdle it is bounded, because the scapula has no bony attachment to the
axial skeleton except the clavicle.  Every newton the hand pushes with reaches
the trunk through the thoracoscapular muscles or through the clavicular strut.

WHAT IS COMPUTED, and every term is measured:

  capacity    Sum(Fmax * |moment arm|) about each girdle coordinate, at the
              model's own rest pose, read from the engine.  Full activation at
              optimal fibre length with no force-velocity and no tendon state, so
              an UPPER BOUND on torque; and a per-axis sum is not a
              tension-feasible torque cone, which is the caution
              docs/UPPER_BODY_ACTUATION.md section 2 already attaches to exactly
              this quantity.
  lever       the distance from the glenohumeral centre to the scapulothoracic
              joint centre, in the torso frame, read from the model's own frames.
              The true lever of the hand load about the scapula is the
              PERPENDICULAR distance to the line of action, which is at most this,
              so using the full distance makes the force bound CONSERVATIVE.
  force       capacity / lever, per side, doubled for two arms.
  demand      the plant's own mass times g from the model's own <gravity>.  A
              press-up pivoting on the toes puts (toe-to-COM / toe-to-hand) of
              body weight on the hands, which is at most ONE body weight, so a
              capacity of one body weight suffices for any pivot geometry.

WHAT THIS IS NOT.  It is not a demonstration that the body presses itself up.
Nothing here integrates, nothing here is a controller, and a capacity computed
from Fmax at one pose is the ceiling of what the muscles could ever deliver, not
what they do deliver.  docs/SHOULDER_GIRDLE.md says so beside the number.
"""
from __future__ import annotations

import json
import shutil
import sys
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402

TARGET_MASS_KG = 77.6122029
GIRDLE = 'data/models/shoulder_girdle_v1'
BASE = 'data/models/articulated_spine_v1'
OUT = GIRDLE + '/press_up_capacity.json'

GIRDLE_COORDINATES = ('clav_prot', 'clav_elev', 'scapula_abduction',
                      'scapula_elevation', 'scapula_upward_rot', 'scapula_winging')
ARM_COORDINATES = ('arm_flex', 'arm_add', 'arm_rot', 'elbow_flex')


def model_gravity(path):
    return np.asarray([float(v) for v in
                       ET.parse(path).getroot().find('.//gravity').text.split()], float)


def max_forces(path):
    out = {}
    for element in ET.parse(path).getroot().iter():
        if element.tag.endswith('Muscle'):
            force = element.findtext('max_isometric_force')
            if force is not None:
                out[element.get('name')] = float(force)
    return out


def joint_frame_translation(path, joint_name, frame_name):
    for joint in ET.parse(path).getroot().iter():
        if joint.get('name') == joint_name and joint.find('frames') is not None:
            for frame in joint.find('frames'):
                if frame.get('name') == frame_name:
                    return np.asarray([float(v) for v in
                                       frame.findtext('translation').split()], float)
    raise KeyError('%s/%s' % (joint_name, frame_name))


def directional_bounds(stream, forces, coordinates):
    """{coordinate: (+bound, -bound, contributing muscles)} in N*m.

    +bound is Sum(Fmax * r) over the muscles with a positive arm, -bound the same
    over the negative ones, exactly as docs/UPPER_BODY_ACTUATION.md section 2
    defines them.
    """
    muscles = sorted(stream.snapshot()['muscles'])
    arms = stream.moment_arms(muscles=muscles, coordinates=list(coordinates))
    repeat = stream.moment_arms(muscles=muscles, coordinates=list(coordinates))
    if arms['moment_arms_m'] != repeat['moment_arms_m']:
        raise RuntimeError('moment_arms is not a function of its input')
    out = {}
    for name in coordinates:
        positive, negative, contributors = 0.0, 0.0, []
        for muscle in muscles:
            arm = arms['moment_arms_m'][muscle][name]
            if abs(arm) < 1e-9:
                continue
            contributors.append({'muscle': muscle, 'moment_arm_m': arm,
                                 'max_isometric_force_n': forces[muscle],
                                 'moment_nm': forces[muscle] * arm})
            if arm > 0:
                positive += forces[muscle] * arm
            else:
                negative += forces[muscle] * -arm
        contributors.sort(key=lambda row: -abs(row['moment_nm']))
        out[name] = {'plus_bound_nm': positive, 'minus_bound_nm': negative,
                     'muscles': len(contributors), 'top': contributors[:6]}
    return out, arms['moment_arms_m']


def session(registration, output):
    return NativeMechanicalStream(ROOT, ROOT / output, environment='free',
                                  target_mass_kg=TARGET_MASS_KG,
                                  augmented_registration=registration)


def measure(root: Path):
    gravity = model_gravity(root / GIRDLE / 'model.osim')
    weight_n = TARGET_MASS_KG * float(np.linalg.norm(gravity))

    # ---- the girdle plant --------------------------------------------
    forces = max_forces(root / GIRDLE / 'model.osim')
    coordinates = (['%s_r' % c for c in GIRDLE_COORDINATES]
                   + ['%s_l' % c for c in GIRDLE_COORDINATES]
                   + ['%s_r' % c for c in ARM_COORDINATES]
                   + ['%s_l' % c for c in ARM_COORDINATES])
    output = 'data/derived/girdle-capacity-' + uuid.uuid4().hex[:10]
    stream = session(GIRDLE + '/registration.json', output)
    try:
        girdle_bounds, _ = directional_bounds(stream, forces, coordinates)
        soleus = stream.moment_arms(muscles=['soleus_r'],
                                    coordinates=['ankle_angle_r'])['moment_arms_m']
        # Segment levers, from the engine's own body frames at the rest pose.
        origin = {b: np.asarray(stream.body_point(body=b, station_m=[0, 0, 0])
                                ['point_source_m'], float)
                  for b in ('humerus_r', 'ulna_r', 'hand_r')}
        lever_shoulder = float(np.linalg.norm(origin['hand_r'] - origin['humerus_r']))
        lever_elbow = float(np.linalg.norm(origin['hand_r'] - origin['ulna_r']))
    finally:
        stream.close()
        shutil.rmtree(root / output, ignore_errors=True)

    # ---- the same plant before the girdle -----------------------------
    base_forces = max_forces(root / BASE / 'model.osim')
    base_coordinates = (['%s_r' % c for c in ARM_COORDINATES]
                        + ['%s_l' % c for c in ARM_COORDINATES])
    output = 'data/derived/girdle-capacity-base-' + uuid.uuid4().hex[:10]
    stream = session(BASE + '/registration.json', output)
    try:
        base_bounds, _ = directional_bounds(stream, base_forces, base_coordinates)
        base_soleus = stream.moment_arms(muscles=['soleus_r'],
                                         coordinates=['ankle_angle_r'])['moment_arms_m']
    finally:
        stream.close()
        shutil.rmtree(root / output, ignore_errors=True)

    # ---- the lever, from the model's own frames -----------------------
    model = root / GIRDLE / 'model.osim'
    record = json.loads((root / GIRDLE / 'registration.json').read_text())
    glenohumeral = np.asarray(record['map']['our_glenohumeral_torso_m'], float)
    scapulothoracic = joint_frame_translation(model, 'scapulothoracic_r',
                                              'torso_st_r_offset')
    sternoclavicular = joint_frame_translation(model, 'sternoclavicular_r',
                                               'torso_sc_r_offset')
    lever_st = float(np.linalg.norm(glenohumeral - scapulothoracic))
    lever_sc = float(np.linalg.norm(glenohumeral - sternoclavicular))

    # ---- the press-up bound -------------------------------------------
    # The binding girdle axis is the one whose capacity, divided by the lever the
    # hand load applies about it, gives the smallest force.  Both directions are
    # reported; the smaller of the two is what a single arm can hold whichever
    # way the load happens to press.
    per_axis = {}
    for name in ['%s_r' % c for c in GIRDLE_COORDINATES]:
        bound = girdle_bounds[name]
        weaker = min(bound['plus_bound_nm'], bound['minus_bound_nm'])
        per_axis[name] = {
            'plus_bound_nm': bound['plus_bound_nm'],
            'minus_bound_nm': bound['minus_bound_nm'],
            'weaker_direction_nm': weaker,
            'force_at_st_lever_n': weaker / lever_st,
        }
    binding = min(per_axis, key=lambda k: per_axis[k]['force_at_st_lever_n'])
    girdle_per_arm_n = per_axis[binding]['force_at_st_lever_n']

    # The rest of the chain, on the same conservative rule: the lever is the full
    # straight-line distance from the joint to the hand, which is what a hand load
    # applied perpendicular to the segment would see. Any other direction sees
    # less, so these are lower bounds on the force each joint can hold.
    chain = {
        'girdle (scapulothoracic, binding axis %s)' % binding: {
            'capacity_nm': per_axis[binding]['weaker_direction_nm'],
            'lever_m': lever_st, 'force_n': girdle_per_arm_n},
        'shoulder (arm_flex_r, weaker direction)': {
            'capacity_nm': min(girdle_bounds['arm_flex_r']['plus_bound_nm'],
                               girdle_bounds['arm_flex_r']['minus_bound_nm']),
            'lever_m': lever_shoulder,
            'force_n': min(girdle_bounds['arm_flex_r']['plus_bound_nm'],
                           girdle_bounds['arm_flex_r']['minus_bound_nm']) / lever_shoulder},
        'elbow (elbow_flex_r, extension)': {
            'capacity_nm': girdle_bounds['elbow_flex_r']['minus_bound_nm'],
            'lever_m': lever_elbow,
            'force_n': girdle_bounds['elbow_flex_r']['minus_bound_nm'] / lever_elbow},
    }
    weakest = min(chain, key=lambda k: chain[k]['force_n'])
    per_arm_n = chain[weakest]['force_n']
    both_arms_n = 2 * per_arm_n

    # TWO readings, and they bracket the answer.
    #
    # Arms straight, hands under the shoulders -- the top of a press-up, and the
    # position the trunk has to be held in. The load line runs through the elbow
    # and the glenohumeral centre, so their levers go to zero and NO arm muscle is
    # in the path. The only muscular demand left is holding the scapula against
    # the thorax, which is the girdle bound and nothing else.
    #
    # Worst-case lever -- the hand force perpendicular to each segment at the
    # segment's full length. That is the hardest a bent arm can ever be loaded,
    # so it is a floor, not the bottom of a press-up specifically.
    readings = {
        'arms_straight_girdle_only': {
            'both_arms_n': 2 * girdle_per_arm_n,
            'body_weights': 2 * girdle_per_arm_n / weight_n,
            'limited_by': 'girdle, axis ' + binding,
            'geometry': 'load line through the elbow and glenohumeral centres, so '
                        'only the scapula needs holding'},
        'worst_case_lever': {
            'both_arms_n': both_arms_n,
            'body_weights': both_arms_n / weight_n,
            'limited_by': weakest,
            'geometry': 'hand force perpendicular to each segment at full segment '
                        'length; a floor on what a bent arm can hold'},
    }

    return {
        'schema': 'ihm.shoulder-girdle-press-up-capacity.v1',
        'measured_by': 'scripts/measure_shoulder_girdle_capacity.py',
        'plant': GIRDLE + '/registration.json',
        'baseline_plant': BASE + '/registration.json',
        'pose': "each model's own rest pose; moment arms are pose dependent and no "
                "press-up pose was evaluated",
        'controls': {
            'soleus_r_about_ankle_angle_r_m': soleus['soleus_r']['ankle_angle_r'],
            'baseline_soleus_r_about_ankle_angle_r_m':
                base_soleus['soleus_r']['ankle_angle_r'],
            'moment_arms_repeat_call_identical': True,
        },
        'gravity_m_s2': gravity.tolist(),
        'body_weight_n': weight_n,
        'levers': {
            'glenohumeral_torso_m': glenohumeral.tolist(),
            'scapulothoracic_torso_m': scapulothoracic.tolist(),
            'sternoclavicular_torso_m': sternoclavicular.tolist(),
            'glenohumeral_to_scapulothoracic_m': lever_st,
            'glenohumeral_to_sternoclavicular_m': lever_sc,
            'glenohumeral_to_hand_m': lever_shoulder,
            'elbow_to_hand_m': lever_elbow,
            'basis': 'straight-line distance, which is at least the perpendicular '
                     'distance from the joint to the hand load line of action, so '
                     'dividing by it makes the force bound conservative',
        },
        'girdle_bounds_nm': girdle_bounds,
        'baseline_arm_bounds_nm': base_bounds,
        'press_up': {
            'per_axis': per_axis,
            'chain': chain,
            'weakest_link': weakest,
            'girdle_per_arm_n': girdle_per_arm_n,
            'readings': readings,
            'binding_axis': binding,
            'per_arm_n': per_arm_n,
            'both_arms_n': both_arms_n,
            'body_weights': both_arms_n / weight_n,
            'baseline_bound_n': None,
            'baseline_reason':
                'In articulated_spine_v1 the humerus is jointed directly to torso, so '
                'at the top of a press-up the hand load reaches the trunk through '
                'joint reactions with no muscle in the path. The quantity is not '
                'unbounded because the plant is strong; it is undefined because '
                'nothing carries it.',
        },
        'not_claimed': [
            'This is capacity, not behaviour. Nothing here integrates and nothing '
            'here is a controller; the body has not been shown to press itself up.',
            'Sum(Fmax * |r|) is full activation at optimal fibre length with no '
            'force-velocity, no tendon state and no activation dynamics: a ceiling.',
            'A per-axis sum is not a tension-feasible torque cone.',
            "Moment arms are pose dependent and these are the rest pose's.",
            'The lever is a straight-line distance, not the perpendicular distance '
            'to a line of action in a press-up pose.',
        ],
    }


def main():
    record = measure(ROOT)
    (ROOT / OUT).write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({'body_weight_n': record['body_weight_n'],
                      'lever_m': record['levers']['glenohumeral_to_scapulothoracic_m'],
                      'press_up': {k: v for k, v in record['press_up'].items()
                                   if k != 'per_axis'},
                      'controls': record['controls']}, indent=2))
    for name, row in record['press_up']['readings'].items():
        print('%-28s %8.1f N = %.3f body weights   limited by %s'
              % (name, row['both_arms_n'], row['body_weights'], row['limited_by']))
    print()
    for name, row in record['press_up']['chain'].items():
        print('%-48s %7.1f N.m / %.4f m = %7.0f N'
              % (name, row['capacity_nm'], row['lever_m'], row['force_n']))
    print()
    for name, row in record['press_up']['per_axis'].items():
        print('%-24s +%.1f / -%.1f N.m   force %.0f N'
              % (name, row['plus_bound_nm'], row['minus_bound_nm'],
                 row['force_at_st_lever_n']))
    print()
    for name in sorted(record['baseline_arm_bounds_nm']):
        a = record['baseline_arm_bounds_nm'][name]
        b = record['girdle_bounds_nm'][name]
        print('%-14s base +%7.1f / -%7.1f   girdle +%7.1f / -%7.1f  N.m'
              % (name, a['plus_bound_nm'], a['minus_bound_nm'],
                 b['plus_bound_nm'], b['minus_bound_nm']))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
