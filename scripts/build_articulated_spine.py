"""Derive data/models/articulated_spine_v1 from engineering_stance_v1.

WHAT THIS ADDS, and where every number comes from.  Nothing here is typed by
hand; each quantity below names the artefact it was read from.

  torso -> thorax -> cervical -> head      three new bodies, nine new rotational
                                           coordinates
  radius_hand_{l,r}  WeldJoint -> UniversalJoint   wrist_flex/wrist_dev
  subtalar_{l,r}     WeldJoint -> PinJoint         subtalar_angle

MASS.  The torso of `engineering_stance_v1` is ONE rigid body from the sacrum to
the crown, 30.3832390808 kg unscaled.  It is repartitioned -- never added to --
into four bodies using the exclusive partition already measured in this repo:

  data/research/cervical_inertia/v2/manifest.json
      convex-envelope inertial priors for C1-C7, skull and jaw, expressed in the
      CURRENT TORSO FRAME, debited from the mass-scaled torso.  Method and its
      limits: docs/research/CERVICAL_INERTIAL_PRIORS.md.
  data/research/thoracic_mechanism/native_composition_v1/plan.json
      the same torso's 48 measured thoracic material shares (ribs 1-12, costal
      cartilages 1-7, manubrium/body/xiphoid, intercostals, diaphragm) and the
      17.937704 kg residual core left when they and the nine cervical bodies are
      taken out.

Those records are in MASS-SCALED units (the engine multiplies every body mass by
target_mass_kg/85.26984854173146 at load).  They are divided back out here so the
XML carries unscaled masses and the engine's own scaling reproduces them.

The four masses sum to the old torso's mass EXACTLY -- that is the check the
verifier runs, not a claim made here.

JOINT CENTRES.
  thoracic   midpoint of the twelfth thoracic and first lumbar vertebral
             centroids (data/derived/canonical/anatomy.json, FJ3156 and FJ3157),
             carried into the torso frame by the SAME canonical->torso transform
             the cervical inertia manifest used.
  neck       the cerv7 body-frame origin of the composition plan, i.e. the
             donor's C7/T1 articulation.
  head       the skull body-frame origin of the composition plan, i.e. the
             donor's atlanto-occipital level.

All three new body frames are AXIS-ALIGNED with the torso (identity rotation).
The registration's own 10.65-degree tilt is therefore NOT carried into the joint
axes; positions come from the registration and axes come from the torso, and the
`back` joint's own axis convention (rot1 z = extension, rot2 x = bending,
rot3 y = axial rotation) is reused unchanged.

RANGES.
  neck_*      MASI auxt1jnt pitch2/roll2/yaw2, the donor's own declared free
              lower-cervical ranges.
  head_*      MASI aux2jnt pitch1/roll1/yaw1, the donor's own declared free
              upper-cervical ranges.  MASI places them at C1/C2; this model lumps
              C1-C7 into one body, so they are applied at the head joint.
              (data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim)
  wrist_*     RajagopalLaiUhlrich2023.osim's own radius_hand UniversalJoint.
  subtalar_*  RajagopalLaiUhlrich2023.osim's own subtalar PinJoint.
  thoracic_*  NOT MEASURED.  An explicit engineering bound, see THORACIC_RANGE.

STOPS.  docs/NATIVE_JOINT_LIMITS.md is why: this model declares ranges, holds no
CoordinateLimitForce and OpenSim does not clamp during forward dynamics, so a
declared range without a stop is a new way for a search to buy travel.  Every new
rotational coordinate gets a CoordinateLimitForce IN THE MODEL at its declared
range, so the stop cannot be forgotten by a caller, plus the viscous damping term
the source model already applies to the same limb.

Run:
  cd <repo> && OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_articulated_spine
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ihm.assembly.cervical_inertia import combine_bodies  # noqa: E402
from ihm.body_parameters import (  # noqa: E402
    MECHANICAL_SOURCE_MASS_KG,
    MECHANICAL_TARGET_MASS_KG,
)

BASE = 'data/models/engineering_stance_v1'
OUT = 'data/models/articulated_spine_v1'
CERVICAL_INERTIA = 'data/research/cervical_inertia/v2/manifest.json'
COMPOSITION = 'data/research/thoracic_mechanism/native_composition_v1/plan.json'
ANATOMY = 'data/derived/canonical/anatomy.json'
MASI = 'data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim'
RAJAGOPAL = ('data/raw/anatomy/opensim-models/source/Models/Rajagopal/'
             'RajagopalLaiUhlrich2023.osim')

#: The engine divides target mass by the sum of the model's body masses.  The
#: research records are in the scaled world; the XML is not.
MASS_SCALE = MECHANICAL_TARGET_MASS_KG / MECHANICAL_SOURCE_MASS_KG

#: T12 and L1 vertebral bodies in data/derived/canonical/anatomy.json.
T12 = 'body-bp3d-FJ3156'
L1 = 'body-bp3d-FJ3157'

#: The joint stop constants scripts/crawl.py already uses, and the only sweep
#: this repo has measured (docs/NATIVE_JOINT_LIMITS.md: 30 N*m/rad is half the
#: wall clock of no stop at all and holds the plant 6.6x closer to its range).
#: They are ENGINEERING CONSTANTS stated by the caller, not ligament properties.
STOP_STIFFNESS_NM_PER_RAD = 30.0
STOP_DAMPING = 1.5
STOP_TRANSITION_RAD = 0.35

#: Viscous damping, transferred from the source model's OWN
#: subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml by limb: every
#: upper-body coordinate there carries -1.0*qdot (shoulder, elbow) and every
#: lower-limb one carries -0.1*qdot (hip, knee, ankle).  Not a measured neck or
#: thoracic damping; the source model has neither.
DAMPING_UPPER = 1.0
DAMPING_LOWER = 0.1

#: THE ONE NUMBER IN THIS MODEL WITH NO MEASURED PROVENANCE.  No donor in this
#: repository declares a thoracic range: MASI welds its rib cage (`ribcagejnt`),
#: the Rajagopal family has no thoracic coordinate, and
#: data/research/thoracic_mechanism/v2/manifest.json states its own rib-angle
#: bound (0.05 rad) as a "declared kinematic exploration bound, not physiological
#: ROM".  These are deliberately NARROWER than any published whole-thoracic range
#: because narrow is the conservative side for a stop whose only job is to keep a
#: search out of absurd configurations.  What would replace it: a per-level
#: thoracic ROM source, or a thoracic coordinate measured from recorded motion --
#: neither of which exists here, because the pose corpus is in gait2392 space and
#: gait2392 has no thoracic degree of freedom either.
THORACIC_RANGE = {'thoracic_extension': (-0.2617993877991494, 0.2617993877991494),
                  'thoracic_bending': (-0.2617993877991494, 0.2617993877991494),
                  'thoracic_rotation': (-0.3490658503988659, 0.3490658503988659)}

#: MASI's own declared ranges, read from its two FREE cervical joints.
NECK_RANGE = {'neck_extension': (-0.575958653158, 0.837758040957),
              'neck_bending': (-0.575958653158, 0.575958653158),
              'neck_rotation': (-0.471238898038, 0.471238898038)}
HEAD_RANGE = {'head_extension': (-0.279252680319, 0.418879020479),
              'head_bending': (-0.10471975512, 0.10471975512),
              'head_rotation': (-0.663225115758, 0.663225115758)}

#: RajagopalLaiUhlrich2023's own declared ranges for the joints this model welds.
WRIST_RANGE = {'wrist_flex': (-1.22173048, 1.22173048),
               'wrist_dev': (-0.43633231, 0.61086523999999998)}
SUBTALAR_RANGE = (-0.61086523819801497, 0.61086523819801497)

#: rotation1/2/3 axes, copied from the `back` joint of this very model.
AXES = ((0, 0, 1), (1, 0, 0), (0, 1, 0))

NEW_BODIES = ('thorax', 'cervical', 'head')
NEW_JOINTS = ('thoracic', 'neck', 'atlantooccipital')


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sub(parent, tag, text=None, **attrib):
    element = ET.SubElement(parent, tag, attrib)
    if text is not None:
        element.text = text
    return element


def triple(values) -> str:
    return ' '.join(repr(float(v)) for v in values)


def from_moments(mass, first, second):
    """(m, integral x dm, integral x x^T dm) about an origin -> body record."""
    mass = float(mass)
    first = np.asarray(first, float)
    second = np.asarray(second, float)
    center = first / mass
    central = second - mass * np.outer(center, center)
    inertia = np.trace(central) * np.eye(3) - central
    return {'mass_kg': mass, 'center_m': center.tolist(),
            'inertia_kg_m2': inertia.tolist()}


def read_partition(root: Path):
    """The four-way exclusive partition of the single torso body."""
    inertia = json.loads((root / CERVICAL_INERTIA).read_text())
    plan = json.loads((root / COMPOSITION).read_text())
    priors = inertia['bodies']
    cervical = combine_bodies([priors['cerv%d' % i]['current_torso_frame_prior']
                               for i in range(1, 8)])
    head = combine_bodies([priors['skull']['current_torso_frame_prior'],
                           priors['jaw']['current_torso_frame_prior']])
    shares = plan['thoracic_material_partitions']
    thorax = from_moments(
        sum(v['mass_kg'] for v in shares.values()),
        sum(np.asarray(v['H_first_moment_kg_m'], float) for v in shares.values()),
        sum(np.asarray(v['Q_second_moment_kg_m2'], float) for v in shares.values()))
    core = plan['thoracic_residual_core']
    torso = from_moments(core['mass_kg'], core['H_first_moment_kg_m'],
                         core['Q_second_moment_kg_m2'])
    return inertia, plan, {'head': head, 'cervical': cervical, 'thorax': thorax,
                           'torso': torso}


def joint_origins(root: Path, inertia, plan):
    """Joint centres in the torso frame, every one read from a measured record."""
    transform = np.asarray(inertia['canonical_to_current_torso'], float)
    entities = {e['id']: e for e in json.loads((root / ANATOMY).read_text())['entities']}
    midpoint = (np.asarray(entities[T12]['centroid_m'], float)
                + np.asarray(entities[L1]['centroid_m'], float)) / 2
    thoracic = transform[:3, :3] @ midpoint + transform[:3, 3]
    neck = np.asarray(plan['cervical_bodies']['cerv7']
                      ['body_to_target_torso_reference'], float)[:3, 3]
    head = np.asarray(plan['cervical_bodies']['skull']
                      ['body_to_target_torso_reference'], float)[:3, 3]
    return {'thoracic': thoracic, 'neck': neck, 'atlantooccipital': head}


def frame_geometry(parent):
    geometry = sub(parent, 'FrameGeometry', name='frame_geometry')
    sub(geometry, 'socket_frame', '..')
    sub(geometry, 'scale_factors', '0.20000000000000001 0.20000000000000001 '
                                   '0.20000000000000001')


def offset_frame(parent, name, socket_parent, translation):
    frame = sub(parent, 'PhysicalOffsetFrame', name=name)
    frame_geometry(frame)
    sub(frame, 'socket_parent', socket_parent)
    sub(frame, 'translation', triple(translation))
    sub(frame, 'orientation', '0 0 0')
    return frame


def make_body(name, record, origin, meshes):
    """A Body whose frame is axis-aligned with torso and offset by `origin`."""
    body = ET.Element('Body', name=name)
    frame_geometry(body)
    if meshes:
        # The .vtp vertices are expressed in the TORSO frame, so a mesh moved to
        # a body whose origin is elsewhere needs a frame that puts it back.  Move
        # the mesh without this and the skull is drawn at the joint centre.
        components = sub(body, 'components')
        offset_frame(components, name + '_display_origin', '..', -np.asarray(origin, float))
        attached = sub(body, 'attached_geometry')
        for mesh in meshes:
            mesh.find('socket_frame').text = '../components/%s_display_origin' % name
            attached.append(mesh)
    else:
        sub(body, 'attached_geometry')
    wrap = sub(body, 'WrapObjectSet', name='wrapobjectset')
    sub(wrap, 'objects')
    sub(wrap, 'groups')
    center = np.asarray(record['center_m'], float) - np.asarray(origin, float)
    inertia = np.asarray(record['inertia_kg_m2'], float) / MASS_SCALE
    sub(body, 'mass', repr(record['mass_kg'] / MASS_SCALE))
    sub(body, 'mass_center', triple(center))
    sub(body, 'inertia', ' '.join(repr(float(v)) for v in (
        inertia[0, 0], inertia[1, 1], inertia[2, 2],
        inertia[0, 1], inertia[0, 2], inertia[1, 2])))
    return body


def coordinate(parent, name, bounds):
    element = sub(parent, 'Coordinate', name=name)
    sub(element, 'default_value', '0')
    sub(element, 'default_speed_value', '0')
    sub(element, 'range', '%s %s' % (repr(float(bounds[0])), repr(float(bounds[1]))))
    sub(element, 'clamped', 'true')
    sub(element, 'locked', 'false')
    sub(element, 'prescribed_function')
    sub(element, 'prescribed', 'false')


def spatial_transform(parent, names):
    transform = sub(parent, 'SpatialTransform')
    for index, (name, axis) in enumerate(zip(names, AXES), start=1):
        rotation = sub(transform, 'TransformAxis', name='rotation%d' % index)
        sub(rotation, 'coordinates', name)
        sub(rotation, 'axis', ' '.join(str(v) for v in axis))
        function = sub(rotation, 'LinearFunction', name='function')
        sub(function, 'coefficients', ' 1 0')
    for index, axis in enumerate(((1, 0, 0), (0, 1, 0), (0, 0, 1)), start=1):
        translation = sub(transform, 'TransformAxis', name='translation%d' % index)
        sub(translation, 'coordinates')
        sub(translation, 'axis', ' '.join(str(v) for v in axis))
        constant = sub(translation, 'Constant', name='function')
        sub(constant, 'value', '0')


def make_joint(name, parent_body, child_body, origin, ranges):
    joint = ET.Element('CustomJoint', name=name)
    sub(joint, 'socket_parent_frame', '%s_offset' % parent_body)
    sub(joint, 'socket_child_frame', '%s_offset' % child_body)
    coordinates = sub(joint, 'coordinates')
    for coordinate_name, bounds in ranges.items():
        coordinate(coordinates, coordinate_name, bounds)
    frames = sub(joint, 'frames')
    offset_frame(frames, '%s_offset' % parent_body, '/bodyset/' + parent_body, origin)
    offset_frame(frames, '%s_offset' % child_body, '/bodyset/' + child_body, (0, 0, 0))
    spatial_transform(joint, list(ranges))
    return joint


def limit_force(name, bounds, stiffness=STOP_STIFFNESS_NM_PER_RAD):
    """A stop AT THE DECLARED RANGE, carried by the model rather than the caller.

    CoordinateLimitForce takes rotational limits and transition in DEGREES and
    its stiffness per degree; everything above is stated in radians.  The force
    is guaranteed zero strictly inside the limits, so a resting body feels it
    only if it is already outside its own declared range.
    """
    degrees = 180.0 / np.pi
    element = ET.Element('CoordinateLimitForce', name='articulated_stop_' + name)
    sub(element, 'coordinate', name)
    sub(element, 'upper_stiffness', repr(stiffness / degrees))
    sub(element, 'upper_limit', repr(float(bounds[1]) * degrees))
    sub(element, 'lower_stiffness', repr(stiffness / degrees))
    sub(element, 'lower_limit', repr(float(bounds[0]) * degrees))
    sub(element, 'damping', repr(STOP_DAMPING))
    sub(element, 'transition', repr(STOP_TRANSITION_RAD * degrees))
    sub(element, 'compute_dissipation_energy', 'false')
    return element


def damping_force(name, coefficient):
    element = ET.Element('ExpressionBasedCoordinateForce',
                         name='articulated_damping_' + name)
    sub(element, 'coordinate', name)
    sub(element, 'expression', '-%s*qdot' % repr(coefficient))
    return element


def find_joint(jointset, name):
    for index, joint in enumerate(jointset):
        if joint.get('name') == name:
            return index, joint
    raise KeyError(name)


def unweld(jointset, name, kind, coordinates):
    """Weld -> real joint, keeping this subject's own scaled offset frames."""
    index, weld = find_joint(jointset, name)
    joint = ET.Element(kind, name=name)
    sub(joint, 'socket_parent_frame', weld.findtext('socket_parent_frame'))
    sub(joint, 'socket_child_frame', weld.findtext('socket_child_frame'))
    element = sub(joint, 'coordinates')
    for coordinate_name, bounds in coordinates:
        coordinate(element, coordinate_name, bounds)
    joint.append(weld.find('frames'))
    jointset[index] = joint
    return joint


def build(root: Path):
    base = root / BASE
    out = root / OUT
    out.mkdir(parents=True, exist_ok=True)

    inertia, plan, partition = read_partition(root)
    origins = joint_origins(root, inertia, plan)

    tree = ET.parse(base / 'model.osim')
    model = tree.getroot().find('Model')
    bodyset = model.find('BodySet/objects')
    jointset = model.find('JointSet/objects')
    forceset = model.find('ForceSet/objects')

    torso = next(b for b in bodyset if b.get('name') == 'torso')
    attached = torso.find('attached_geometry')
    moved = {'head': [], 'thorax': []}
    for mesh in list(attached):
        source = mesh.findtext('mesh_file')
        if source in ('hat_skull.vtp', 'hat_jaw.vtp'):
            attached.remove(mesh)
            moved['head'].append(mesh)
        elif source == 'hat_ribs_scap.vtp':
            attached.remove(mesh)
            moved['thorax'].append(mesh)

    record = partition['torso']
    torso.find('mass').text = repr(record['mass_kg'] / MASS_SCALE)
    torso.find('mass_center').text = triple(record['center_m'])
    tensor = np.asarray(record['inertia_kg_m2'], float) / MASS_SCALE
    torso.find('inertia').text = ' '.join(repr(float(v)) for v in (
        tensor[0, 0], tensor[1, 1], tensor[2, 2],
        tensor[0, 1], tensor[0, 2], tensor[1, 2]))

    body_origins = {'thorax': origins['thoracic'], 'cervical': origins['neck'],
                    'head': origins['atlantooccipital']}
    for name in NEW_BODIES:
        bodyset.append(make_body(name, partition[name], body_origins[name],
                                 moved.get(name, [])))

    jointset.append(make_joint('thoracic', 'torso', 'thorax',
                               origins['thoracic'], THORACIC_RANGE))
    jointset.append(make_joint('neck', 'thorax', 'cervical',
                               origins['neck'] - origins['thoracic'], NECK_RANGE))
    jointset.append(make_joint('atlantooccipital', 'cervical', 'head',
                               origins['atlantooccipital'] - origins['neck'],
                               HEAD_RANGE))

    stops = []
    for ranges, damping in ((THORACIC_RANGE, DAMPING_UPPER),
                            (NECK_RANGE, DAMPING_UPPER),
                            (HEAD_RANGE, DAMPING_UPPER)):
        for name, bounds in ranges.items():
            stops.append((name, bounds, damping))
    for side in ('r', 'l'):
        unweld(jointset, 'radius_hand_%s' % side, 'UniversalJoint',
               [('wrist_flex_%s' % side, WRIST_RANGE['wrist_flex']),
                ('wrist_dev_%s' % side, WRIST_RANGE['wrist_dev'])])
        stops.append(('wrist_flex_%s' % side, WRIST_RANGE['wrist_flex'], DAMPING_UPPER))
        stops.append(('wrist_dev_%s' % side, WRIST_RANGE['wrist_dev'], DAMPING_UPPER))
        unweld(jointset, 'subtalar_%s' % side, 'PinJoint',
               [('subtalar_angle_%s' % side, SUBTALAR_RANGE)])
        stops.append(('subtalar_angle_%s' % side, SUBTALAR_RANGE, DAMPING_LOWER))

    for name, bounds, damping in stops:
        forceset.append(limit_force(name, bounds))
        forceset.append(damping_force(name, damping))

    ET.indent(tree, space='\t')
    model_path = out / 'model.osim'
    tree.write(model_path, encoding='utf-8', xml_declaration=True)
    shutil.copyfile(base / 'catalog.json', out / 'catalog.json')
    shutil.copyfile(base / 'initial_pose.json', out / 'initial_pose.json')
    shutil.copyfile(base / 'equilibrium_excitations.json',
                    out / 'equilibrium_excitations.json')

    sources = {}
    for relative in (BASE + '/model.osim', BASE + '/catalog.json', CERVICAL_INERTIA,
                     COMPOSITION, ANATOMY, MASI, RAJAGOPAL,
                     'scripts/build_articulated_spine.py',
                     'ihm/assembly/cervical_inertia.py'):
        sources[relative] = sha(root / relative)
    sources[OUT + '/model.osim'] = sha(model_path)
    sources[OUT + '/catalog.json'] = sha(out / 'catalog.json')

    registration = {
        'schema': 'ihm.articulated-spine-variant.v1',
        'model_path': OUT + '/model.osim',
        'model_sha256': sha(model_path),
        'catalog_path': OUT + '/catalog.json',
        'catalog_sha256': sha(out / 'catalog.json'),
        'base_model_path': BASE + '/model.osim',
        'base_model_sha256': sha(base / 'model.osim'),
        'sources': sources,
        'muscle_count': 98,
        'body_count': 25,
        'default_enabled': False,
        'native_acceptance_complete': False,
        'scope': ('Adds head/cervical/thoracic bodies, wrists and subtalar joints to '
                  'engineering_stance_v1 by repartitioning its single torso body. No '
                  'identified linearization, no stance acceptance, no muscle crosses '
                  'any new joint: this is a kinematic and inertial variant only.'),
        'not_claimed': [
            'No linearization.npz: the identified stance plant of engineering_stance_v1 '
            'was solved on 22 bodies and 33 coordinates and does not transfer.',
            'The nine spine/neck coordinates carry no actuator and no muscle. They are '
            'free coordinates with a passive stop and viscous damping only.',
            'thoracic_* ranges are an explicit engineering bound, not a measurement.',
            'Scapulae and clavicles are not in the thoracic material partition, so the '
            'shoulder girdle and both arms still ride torso, not thorax.',
        ],
        'mass_scale_used': MASS_SCALE,
        'joint_centres_torso_frame_m': {k: v.tolist() for k, v in origins.items()},
        'partition_kg_mass_scaled': {k: v['mass_kg'] for k, v in partition.items()},
    }
    (out / 'registration.json').write_text(json.dumps(registration, indent=2) + '\n')
    return registration


if __name__ == '__main__':
    result = build(ROOT)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('sources', 'not_claimed')}, indent=2))
