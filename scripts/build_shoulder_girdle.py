"""Derive data/models/shoulder_girdle_v1 from data/models/articulated_spine_v1.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_shoulder_girdle

WHAT THIS ADDS.  A shoulder girdle: `clavicle_{r,l}` and `scapula_{r,l}`, a
sternoclavicular and a scapulothoracic joint per side, the acromioclavicular
point constraint that closes each loop, and the 30 donor muscles per side that
live on those bodies -- deltoid, rotator cuff, pectoralis major and minor,
latissimus dorsi, trapezius, serratus anterior, rhomboids, levator scapulae,
teres major and minor, coracobrachialis.  `acromial_{r,l}` is re-parented from
`torso` onto `scapula_{r,l}`, so the arm now hangs off the girdle and the load
path from the hand to the trunk runs through the thoracoscapular muscles, which
is the whole point.

THE DONOR, and why it and not the other one.

  data/raw/mechanics/opensim-core/OpenSim/Tests/shared/ThoracoscapularShoulderModel.osim
  Seth, Dong, Matias and Delp 2019, Front. Neurorobot. 13:90.  CC BY 4.0 from the
  authors on SimTK (simtk.org/projects/thoracoscapular), Apache-2.0 as a file of
  the opensim-core repository it is read from here.  Both permit commercial use
  with attribution.  MoBL-ARMS 4.1, the other girdle donor on disk, is BSD-3
  RESTRICTED TO NON-COMMERCIAL USE by its SimTK notice
  (docs/UPPER_BODY_ACTUATION.md section 10) and is not used here, directly or
  indirectly.  Nothing in this file, in the model it writes or in the numbers it
  carries comes from MoBL-ARMS.

  Provenance, completeness and what could NOT be established:
  docs/SHOULDER_GIRDLE.md section 1.

THE MAP between the donor's frames and ours.  Two similarity factors, each
measured from a correspondence the two models both declare, and one translation:

  s_lat   this plant's glenohumeral half-width divided by the donor's, read from
          `acromial_r`'s own parent offset frame on `torso` and from the donor's
          glenohumeral centre in its thorax frame.  It scales the GIRDLE: the
          clavicle and the scapula span the midline to the glenoid, which is the
          distance this factor matches by construction.  Because the midline maps
          to the midline (the translation's z is exactly zero), every midline
          donor attachment -- trapezius and rhomboid origins on the spinous
          processes -- stays on the midline.
  s_long  this plant's glenohumeral-to-elbow distance divided by the donor's,
          both read from the models' own joint frames.  It scales the points that
          land on the HUMERUS, which is this plant's own body and is not being
          replaced.
  delta   chosen so the donor's glenohumeral centre lands exactly on this plant's.
          The arm therefore does not move: the verifier checks that against the
          base model through the engine, rather than trusting this sentence.

  Donor body poses in the donor thorax frame are not readable from its XML -- the
  scapulothoracic joint is a Simbody ellipsoid mobilizer and the clavicle/scapula
  loop is closed by a PointConstraint -- so they are MEASURED with the engine by
  scripts/measure_shoulder_donor_frames.py and read from its record here.

MASS.  The girdle is DEBITED from `torso`, never added, through
`ihm.assembly.cervical_inertia.partition_body`, which is the same instrument the
cervical/thoracic partition used and which composes with it exactly.
docs/ARTICULATED_SPINE.md says why this is the second half of one partition: its
thoracic share is rib cage only and has no mass for a scapula, so the girdle
comes out of the 17.938 kg residual core that is still inside `torso`.

  The donor's clavicle (0.1898 kg) and scapula (0.5016 kg) masses are transferred
  VERBATIM and are NOT scaled to this body, because the donor model states no body
  mass to scale from.  That is the same stance this repository already takes for
  donor muscle forces.  Their inertia tensors carry the s_lat^2 the geometry
  carries, so the tensor matches the size the body is modelled at.

RANGES.  The donor declares +-pi for five of the six girdle coordinates and
+-pi/2 for the sixth; neither is a range, both are the placeholder this plant's
own +-10 rad shoulder coordinates are.  No donor on disk declares a scapular or
clavicular range of motion, and the donor's shipped motion file is synthetic
(0.01, 0.02, 0.03 ... per column).  So the ranges here are an EXPLICIT
ENGINEERING BOUND with no measured provenance, chosen narrow because narrow is
the conservative side for a stop whose only job is to keep a search out of absurd
configurations -- exactly the position docs/ARTICULATED_SPINE.md takes for its
thoracic range, and labelled the same way.  Every one carries a
CoordinateLimitForce at that range IN THE MODEL FILE, because
docs/NATIVE_JOINT_LIMITS.md says nothing else enforces a declared range.

THE REST POSE is the donor's own girdle configuration: each new coordinate's
default is the donor's default, not zero, because the acromioclavicular
constraint is satisfied exactly there and at no value this file could invent.

THE LEFT SIDE IS A MIRROR PRIOR.  The donor has one arm.  Everything on the left
is the right side reflected in z.  docs/UPPER_BODY_ACTUATION.md section 8.5 warns
that this needs explicit validation of vectors, joint axes and wrap quadrants and
is not licensed by a name suffix; scripts/verify_shoulder_girdle.py measures it
and docs/SHOULDER_GIRDLE.md records what the measurement found.
"""
from __future__ import annotations

import copy
import hashlib
import json
import shutil
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ihm.assembly.cervical_inertia import (  # noqa: E402
    convex_geometry_prior,
    homogeneous_prior,
    partition_body,
)
from ihm.body_parameters import (  # noqa: E402
    MECHANICAL_SOURCE_MASS_KG,
    MECHANICAL_TARGET_MASS_KG,
)

BASE = 'data/models/articulated_spine_v1'
BASE_MODEL = BASE + '/model.osim'
BASE_CATALOG = BASE + '/catalog.json'
BASE_REGISTRATION = BASE + '/registration.json'
OUT = 'data/models/shoulder_girdle_v1'
DONOR = ('data/raw/mechanics/opensim-core/OpenSim/Tests/shared/'
         'ThoracoscapularShoulderModel.osim')
DONOR_FRAMES = OUT + '/donor_reference_frames.json'
CERVICAL_INERTIA = 'data/research/cervical_inertia/v2/manifest.json'
GEOMETRY = 'data/derived/canonical/geometry/%s.json.gz'

#: THE DONOR'S OWN CLAVICLE AND SCAPULA INERTIA TENSORS ARE NOT PHYSICALLY
#: REALISABLE and are therefore not used.  Measured here from the donor file:
#: the second-moment matrix tr(I)/2*eye - I has a NEGATIVE eigenvalue for the
#: clavicle (-3.07e-05), the scapula (-9.79e-05) and the radius (-1.40e-04), so
#: the principal moments violate the triangle inequality; the radius even has a
#: negative principal moment (-4.54e-06 kg m^2).  No rigid mass distribution has
#: those tensors.  `ihm.assembly.cervical_inertia._physical` refuses them, which
#: is how this was found, and the repository's standing decision for exactly this
#: case is on record: 'Replace donor tensor with geometry-conditioned prior; no
#: eigenvalue clipping and no Body.cpp fallback'
#: (data/research/cervical_inertia/v2/manifest.json, every cervical body).
#: So each new body's inertia is the convex-envelope prior of THIS body's own
#: anatomical scapula and clavicle, at the donor's stated mass, carried into the
#: new body's axes by the same canonical->torso rotation the cervical partition
#: used.  A rotation, not a translation: a central second moment does not depend
#: on where the mesh sits, so that registration's 62 mm positional residual does
#: not enter.  Each side uses its OWN side's mesh -- no mirror prior here.
INERTIA_GEOMETRY = {'clavicle_r': 'body-bp3d-FJ3362', 'scapula_r': 'body-bp3d-FJ3384',
                    'clavicle_l': 'body-bp3d-FJ3237', 'scapula_l': 'body-bp3d-FJ3279'}

#: The engine divides target mass by the sum of the model's body masses, so the
#: XML carries unscaled masses (build_articulated_spine.py, same constant).
MASS_SCALE = MECHANICAL_TARGET_MASS_KG / MECHANICAL_SOURCE_MASS_KG

#: Donor bodies that exist in this plant after the build, and what they become.
BODY_MAP = {'thorax': 'torso', 'clavicle': 'clavicle', 'scapula': 'scapula',
            'humerus': 'humerus'}
#: Donor bodies this plant does not take from this donor. A muscle touching one
#: is not transferred. ulna/radius carry the donor's own elbow muscles, which
#: would duplicate the arm26 TRIlong/BIClong/BICshort already in this plant
#: (docs/UPPER_BODY_ACTUATION.md section 8.4: adding girdle muscle does not
#: license duplicating or silently replacing them).
DONOR_EXCLUDED_BODIES = ('ulna', 'radius', 'hand')

#: The four new bodies, and which donor body each comes from.
NEW_BODIES = (('clavicle_r', 'clavicle', +1), ('scapula_r', 'scapula', +1),
              ('clavicle_l', 'clavicle', -1), ('scapula_l', 'scapula', -1))

#: Half-width of every new coordinate's declared range about the donor's own
#: default value. NOT MEASURED -- see the module docstring.
GIRDLE_RANGE_HALFWIDTH_RAD = 0.3490658503988659

#: The joint stop constants scripts/crawl.py uses and build_articulated_spine.py
#: wrote into the spine variant; the only stiffness this repository has swept
#: (docs/NATIVE_JOINT_LIMITS.md). Engineering constants, not ligament properties.
STOP_STIFFNESS_NM_PER_RAD = 30.0
STOP_DAMPING = 1.5
STOP_TRANSITION_RAD = 0.35
#: Viscous damping: the source model's own upper-body value (-1.0*qdot), the one
#: build_articulated_spine.py transferred to every upper-body coordinate.
DAMPING_UPPER = 1.0


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sub(parent, tag, text=None, **attrib):
    element = ET.SubElement(parent, tag, attrib)
    if text is not None:
        element.text = text
    return element


def triple(values):
    return ' '.join(repr(float(v)) for v in values)


def text_of(element, tag):
    found = element.find(tag)
    return None if found is None else (found.text or '').strip()


def numbers(element, tag):
    return np.asarray([float(v) for v in text_of(element, tag).split()], float)


def inertia_matrix(values):
    xx, yy, zz, xy, xz, yz = values
    return np.asarray([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], float)


def inertia_text(matrix):
    return ' '.join(repr(float(v)) for v in (matrix[0, 0], matrix[1, 1], matrix[2, 2],
                                             matrix[0, 1], matrix[0, 2], matrix[1, 2]))


def euler_xyz(rotation):
    """Body-fixed XYZ Euler angles, the convention OpenSim's `orientation` uses."""
    rotation = np.asarray(rotation, float)
    ry = np.arcsin(np.clip(rotation[0, 2], -1.0, 1.0))
    rx = np.arctan2(-rotation[1, 2], rotation[2, 2])
    rz = np.arctan2(-rotation[0, 1], rotation[0, 0])
    check = rotation_xyz((rx, ry, rz))
    if np.abs(check - rotation).max() > 1e-12:
        raise ValueError('Euler decomposition does not reproduce the rotation')
    return np.asarray([rx, ry, rz], float)


def rotation_xyz(angles):
    rx, ry, rz = (float(a) for a in angles)
    cx, sx = np.cos(rx), np.sin(rx)
    cy, sy = np.cos(ry), np.sin(ry)
    cz, sz = np.cos(rz), np.sin(rz)
    x = np.asarray([[1, 0, 0], [0, cx, -sx], [0, sx, cx]], float)
    y = np.asarray([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], float)
    z = np.asarray([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]], float)
    return x @ y @ z


MIRROR = np.diag([1.0, 1.0, -1.0])


def mirrored(rotation):
    """The reflection-conjugate of a rotation: proper, and det stays +1."""
    out = MIRROR @ np.asarray(rotation, float) @ MIRROR
    if abs(np.linalg.det(out) - 1) > 1e-12:
        raise ValueError('Mirrored rotation is not proper')
    return out


class Registration:
    """The donor -> this plant coordinate map, every factor measured."""

    def __init__(self, root: Path):
        self.root = root
        self.donor = ET.parse(root / DONOR).getroot().find('Model')
        self.base = ET.parse(root / BASE_MODEL)
        frames = json.loads((root / DONOR_FRAMES).read_text())
        if frames['donor_sha256'] != sha(root / DONOR):
            raise ValueError('donor_reference_frames.json was measured on other bytes')
        self.frames = {name: np.asarray(row['body_to_thorax_reference'], float)
                       for name, row in frames['frames'].items()}
        self.frame_record = frames

        # --- the two correspondences, read from the two models' own joint frames
        gh_in_scapula = self._donor_frame_translation('GlenoHumeral', 'scapula_offset')
        self.donor_gh_thorax = (self.frames['scapula'][:3, :3] @ gh_in_scapula
                                + self.frames['scapula'][:3, 3])
        self.donor_gh_in_scapula = gh_in_scapula
        donor_elbow = self._donor_frame_translation('elbow', 'humerus_offset')
        self.donor_humerus_length = float(np.linalg.norm(donor_elbow))

        acromial = self._base_joint('acromial_r')
        self.our_gh_torso = self._frame_translation(acromial, 'torso_offset')
        elbow = self._base_joint('elbow_r')
        self.our_humerus_length = float(np.linalg.norm(
            self._frame_translation(elbow, 'humerus_r_offset')))

        self.s_lat = float(self.our_gh_torso[2] / self.donor_gh_thorax[2])
        self.s_long = float(self.our_humerus_length / self.donor_humerus_length)
        self.delta = np.asarray([
            self.our_gh_torso[0] - self.s_lat * self.donor_gh_thorax[0],
            self.our_gh_torso[1] - self.s_lat * self.donor_gh_thorax[1],
            0.0], float)
        if abs(self.s_lat * self.donor_gh_thorax[2] + self.delta[2]
               - self.our_gh_torso[2]) > 1e-15:
            raise ValueError('Lateral scale does not reproduce the glenohumeral width')
        #: The scapula's reference orientation in the torso frame. The acromial
        #: joint's new parent offset frame carries its inverse, so the humerus
        #: frame at zero arm coordinates is exactly the base plant's.
        self.scapula_rotation = self.frames['scapula'][:3, :3]
        self.clavicle_rotation = self.frames['clavicle'][:3, :3]

    # -- readers -----------------------------------------------------------
    def _donor_joint(self, name):
        for joint in self.donor.find('JointSet/objects'):
            if joint.get('name') == name:
                return joint
        raise KeyError(name)

    def _base_joint(self, name):
        for joint in self.base.getroot().find('.//JointSet/objects'):
            if joint.get('name') == name:
                return joint
        raise KeyError(name)

    @staticmethod
    def _frame_translation(joint, frame_name):
        for frame in joint.find('frames'):
            if frame.get('name') == frame_name:
                return numbers(frame, 'translation')
        raise KeyError(frame_name)

    def _donor_frame_translation(self, joint, frame_name):
        return self._frame_translation(self._donor_joint(joint), frame_name)

    # -- the map -----------------------------------------------------------
    def thorax_to_torso(self, point, side=+1):
        """A donor THORAX-frame point in this plant's torso frame."""
        out = self.s_lat * np.asarray(point, float) + self.delta
        return out if side > 0 else MIRROR @ out

    def girdle_local(self, point, side=+1):
        """A donor clavicle/scapula-frame point in the same body's frame here.

        The body frames keep the donor's orientation and the map is a similarity,
        so a local station only scales -- and mirrors, on the left.
        """
        out = self.s_lat * np.asarray(point, float)
        return out if side > 0 else MIRROR @ out

    def humerus_local(self, point, side=+1):
        """A donor HUMERUS-frame point in this plant's humerus frame.

        The donor's humerus frame at zero glenohumeral coordinates is aligned
        with its scapula; this plant's is aligned with `torso`. The difference is
        the scapula's reference rotation, and the humerus is this plant's own
        body, so the length factor is s_long, not s_lat.
        """
        out = self.s_long * (self.scapula_rotation @ np.asarray(point, float))
        return out if side > 0 else MIRROR @ out

    def body_origin(self, donor_body, side=+1):
        transform = self.frames[donor_body]
        return self.thorax_to_torso(transform[:3, 3], side)

    def body_rotation(self, donor_body, side=+1):
        rotation = self.frames[donor_body][:3, :3]
        return rotation if side > 0 else mirrored(rotation)


def anatomical_inertia(root: Path, body_name, mass_kg, body_rotation, canonical_rotation):
    """Convex-envelope inertia of THIS body's own scapula/clavicle, in body axes."""
    import gzip

    entity = INERTIA_GEOMETRY[body_name]
    with gzip.open(root / (GEOMETRY % entity)) as handle:
        mesh = json.load(handle)
    vertices = np.asarray(mesh['positions'], float).reshape(-1, 3)
    faces = np.asarray(mesh['indices'], int).reshape(-1, 3)
    geometry = convex_geometry_prior(vertices, faces)
    prior = homogeneous_prior(geometry, mass_kg)
    canonical = np.asarray(prior['inertia_kg_m2'], float)
    # canonical axes -> torso axes -> this body's axes.
    in_torso = canonical_rotation @ canonical @ canonical_rotation.T
    return body_rotation.T @ in_torso @ body_rotation, geometry, prior, entity


def girdle_body_record(root, registration, body_name, donor_body, side, canonical_rotation):
    """Mass, COM and inertia of one new body, in its own frame and in the torso frame."""
    donor = None
    for body in registration.donor.find('BodySet/objects'):
        if body.get('name') == donor_body:
            donor = body
    if donor is None:
        raise KeyError(donor_body)
    physical_mass = float(text_of(donor, 'mass'))
    mass = physical_mass / MASS_SCALE
    com_local = registration.girdle_local(numbers(donor, 'mass_center'), side)
    rotation = registration.body_rotation(donor_body, side)
    origin = registration.body_origin(donor_body, side)
    tensor, geometry, prior, entity = anatomical_inertia(
        root, body_name, physical_mass, rotation, canonical_rotation)
    inertia_local = tensor / MASS_SCALE
    local = {'mass_kg': mass, 'center_m': com_local.tolist(),
             'inertia_kg_m2': inertia_local.tolist()}
    torso_frame = {
        'mass_kg': mass,
        'center_m': (rotation @ com_local + origin).tolist(),
        'inertia_kg_m2': (rotation @ inertia_local @ rotation.T).tolist()}
    evidence = {
        'donor_mass_kg': physical_mass,
        'donor_mass_scaled_to_this_body': False,
        'donor_inertia_kg_m2': numbers(donor, 'inertia').tolist(),
        'donor_inertia_principal_moments': sorted(
            np.linalg.eigvalsh(inertia_matrix(numbers(donor, 'inertia'))).tolist()),
        'donor_inertia_second_moment_min_eigenvalue': float(np.linalg.eigvalsh(
            np.trace(inertia_matrix(numbers(donor, 'inertia'))) / 2 * np.eye(3)
            - inertia_matrix(numbers(donor, 'inertia'))).min()),
        'donor_inertia_physically_realisable': False,
        'replacement_geometry_entity': entity,
        'replacement_basis': prior['basis'],
        'replacement_envelope_volume_m3': geometry['volume_m3'],
        'replacement_effective_density_kg_m3': prior['effective_density_kg_m3'],
        'replacement_principal_moments': sorted(
            np.linalg.eigvalsh(np.asarray(prior['inertia_kg_m2'], float)).tolist()),
        'hull_to_aabb_volume_ratio':
            geometry['geometry_error']['hull_to_aabb_volume_ratio'],
    }
    return local, torso_frame, origin, rotation, evidence


# --------------------------------------------------------------------------
# XML emitters
# --------------------------------------------------------------------------
def frame_geometry(parent):
    geometry = sub(parent, 'FrameGeometry', name='frame_geometry')
    sub(geometry, 'socket_frame', '..')
    sub(geometry, 'scale_factors', '0.20000000000000001 0.20000000000000001 '
                                   '0.20000000000000001')


def offset_frame(parent, name, socket_parent, translation, orientation=(0, 0, 0)):
    frame = sub(parent, 'PhysicalOffsetFrame', name=name)
    frame_geometry(frame)
    sub(frame, 'socket_parent', socket_parent)
    sub(frame, 'translation', triple(translation))
    sub(frame, 'orientation', triple(orientation))
    return frame


def make_body(name, local, wraps):
    body = ET.Element('Body', name=name)
    frame_geometry(body)
    sub(body, 'attached_geometry')
    wrap_set = sub(body, 'WrapObjectSet', name='wrapobjectset')
    objects = sub(wrap_set, 'objects')
    for wrap in wraps:
        objects.append(wrap)
    sub(wrap_set, 'groups')
    sub(body, 'mass', repr(float(local['mass_kg'])))
    sub(body, 'mass_center', triple(local['center_m']))
    sub(body, 'inertia', inertia_text(np.asarray(local['inertia_kg_m2'], float)))
    return body


def coordinate(parent, name, bounds, default):
    element = sub(parent, 'Coordinate', name=name)
    sub(element, 'default_value', repr(float(default)))
    sub(element, 'default_speed_value', '0')
    sub(element, 'range', '%s %s' % (repr(float(bounds[0])), repr(float(bounds[1]))))
    sub(element, 'clamped', 'true')
    sub(element, 'locked', 'false')
    sub(element, 'prescribed_function')
    sub(element, 'prescribed', 'false')


def limit_force(name, bounds):
    """A stop AT the declared range, in the model file rather than the caller.

    CoordinateLimitForce takes rotational limits, stiffness and transition in
    OpenSim's degree convention; the constants above are stated in radians. This
    is the identical conversion build_articulated_spine.py writes, including the
    damping units that docs/ARTICULATED_SPINE.md records as disagreeing with
    their declared value -- the two paths agree with each other, which is what
    matters for a plant meant to compare with that one.
    """
    degrees = 180.0 / np.pi
    element = ET.Element('CoordinateLimitForce', name='girdle_stop_' + name)
    sub(element, 'coordinate', name)
    sub(element, 'upper_stiffness', repr(STOP_STIFFNESS_NM_PER_RAD / degrees))
    sub(element, 'upper_limit', repr(float(bounds[1]) * degrees))
    sub(element, 'lower_stiffness', repr(STOP_STIFFNESS_NM_PER_RAD / degrees))
    sub(element, 'lower_limit', repr(float(bounds[0]) * degrees))
    sub(element, 'damping', repr(STOP_DAMPING))
    sub(element, 'transition', repr(STOP_TRANSITION_RAD * degrees))
    sub(element, 'compute_dissipation_energy', 'false')
    return element


def damping_force(name):
    element = ET.Element('ExpressionBasedCoordinateForce',
                         name='girdle_damping_' + name)
    sub(element, 'coordinate', name)
    sub(element, 'expression', '-%s*qdot' % repr(DAMPING_UPPER))
    return element


def scaled_wrap(registration, wrap, donor_body, side):
    """A donor wrap object moved onto the body that now carries it.

    torso    the donor's thorax ellipsoid: translated by delta and scaled by s_lat.
    scapula  a similarity in the body's own frame: translation and size scale.
    humerus  this plant's OWN body, so the donor's humerus-frame quantities are
             rotated by the scapula's reference rotation and scaled by s_long.
    Quadrant strings survive the z-mirror unchanged: the mirrored wrap frame's
    axes are the mirror images of the right side's, so the same half-space is
    named by the same string.
    """
    out = copy.deepcopy(wrap)
    translation = numbers(wrap, 'translation')
    rotation = rotation_xyz(numbers(wrap, 'xyz_body_rotation'))
    if donor_body == 'thorax':
        factor = registration.s_lat
        translation = registration.thorax_to_torso(translation, side)
        rotation = rotation if side > 0 else mirrored(rotation)
    elif donor_body == 'scapula':
        factor = registration.s_lat
        translation = registration.girdle_local(translation, side)
        rotation = rotation if side > 0 else mirrored(rotation)
    elif donor_body == 'humerus':
        factor = registration.s_long
        translation = registration.humerus_local(translation, side)
        rotation = registration.scapula_rotation @ rotation
        rotation = rotation if side > 0 else mirrored(rotation)
    else:
        raise ValueError('No wrap transfer rule for donor body %r' % donor_body)
    out.find('translation').text = triple(translation)
    out.find('xyz_body_rotation').text = triple(euler_xyz(rotation))
    for tag in ('radius', 'length'):
        node = out.find(tag)
        if node is not None:
            node.text = repr(factor * float(node.text))
    node = out.find('dimensions')
    if node is not None:
        node.text = triple(factor * numbers(wrap, 'dimensions'))
    out.set('name', wrap.get('name') + '_' + ('r' if side > 0 else 'l'))
    return out


def donor_wraps(registration):
    """{donor body: [wrap element]} for every wrap a transferred muscle uses."""
    used = set()
    for muscle in registration.donor.iter('Millard2012EquilibriumMuscle'):
        for wrap in muscle.iter('PathWrap'):
            used.add(text_of(wrap, 'wrap_object'))
    out = {}
    for body in registration.donor.find('BodySet/objects'):
        for wrap_set in body.iter('WrapObjectSet'):
            for wrap in wrap_set.find('objects'):
                if wrap.get('name') in used:
                    out.setdefault(body.get('name'), []).append(wrap)
    return out


def sternoclavicular_joint(registration, side):
    """The donor's own sternoclavicular CustomJoint, carried over.

    The axes are the donor's. On the left each rotation axis is -MIRROR*a: a
    positive coordinate then produces the MIRROR IMAGE of the right side's
    motion, which is the convention this plant's own `acromial_l` already uses
    (its axes are exactly -MIRROR of `acromial_r`'s).
    """
    suffix = 'r' if side > 0 else 'l'
    donor = registration._donor_joint('sternoclavicular')
    joint = ET.Element('CustomJoint', name='sternoclavicular_' + suffix)
    sub(joint, 'socket_parent_frame', 'torso_sc_%s_offset' % suffix)
    sub(joint, 'socket_child_frame', 'clavicle_%s_offset' % suffix)
    coordinates = sub(joint, 'coordinates')
    ranges = {}
    for element in donor.iter('Coordinate'):
        name = '%s_%s' % (element.get('name'), suffix)
        default = float(text_of(element, 'default_value'))
        bounds = (default - GIRDLE_RANGE_HALFWIDTH_RAD,
                  default + GIRDLE_RANGE_HALFWIDTH_RAD)
        coordinate(coordinates, name, bounds, default)
        ranges[name] = bounds
    frames = sub(joint, 'frames')
    origin = registration.body_origin('clavicle', side)
    offset_frame(frames, 'torso_sc_%s_offset' % suffix, '/bodyset/torso', origin)
    offset_frame(frames, 'clavicle_%s_offset' % suffix,
                 '/bodyset/clavicle_%s' % suffix, (0, 0, 0))
    transform = copy.deepcopy(donor.find('SpatialTransform'))
    for axis in transform.findall('TransformAxis'):
        names = text_of(axis, 'coordinates')
        if names:
            axis.find('coordinates').text = '%s_%s' % (names, suffix)
        if axis.get('name').startswith('rotation') and side < 0:
            axis.find('axis').text = triple(-(MIRROR @ numbers(axis, 'axis')))
    joint.append(transform)
    return joint, ranges


def scapulothoracic_joint(registration, side, left_signs):
    """The donor's ScapulothoracicJoint: a Simbody ellipsoid mobilizer.

    The ellipsoid is the donor's, scaled by s_lat, on `torso` instead of on the
    donor's `thorax`. On the left the parent frame is the mirror-conjugate; the
    coordinate SIGN convention that produces a mirror-image rest pose is not
    derivable from the joint's declaration, so it is MEASURED
    (scripts/verify_shoulder_girdle.py gates the mirror to 1e-9 m) and the
    measured signs are `left_signs`.
    """
    suffix = 'r' if side > 0 else 'l'
    donor = registration._donor_joint('scapulothoracic')
    joint = ET.Element('ScapulothoracicJoint', name='scapulothoracic_' + suffix)
    sub(joint, 'socket_parent_frame', 'torso_st_%s_offset' % suffix)
    sub(joint, 'socket_child_frame', 'scapula_%s_offset' % suffix)
    coordinates = sub(joint, 'coordinates')
    ranges = {}
    for element in donor.iter('Coordinate'):
        base_name = element.get('name')
        name = '%s_%s' % (base_name, suffix)
        sign = 1.0 if side > 0 else float(left_signs[base_name])
        default = sign * float(text_of(element, 'default_value'))
        bounds = (default - GIRDLE_RANGE_HALFWIDTH_RAD,
                  default + GIRDLE_RANGE_HALFWIDTH_RAD)
        coordinate(coordinates, name, bounds, default)
        ranges[name] = bounds
    sub(joint, 'thoracic_ellipsoid_radii_x_y_z',
        triple(registration.s_lat * numbers(donor, 'thoracic_ellipsoid_radii_x_y_z')))
    sub(joint, 'scapula_winging_axis_origin',
        ' '.join(repr(registration.s_lat * float(v))
                 for v in text_of(donor, 'scapula_winging_axis_origin').split()))
    sub(joint, 'scapula_winging_axis_direction',
        text_of(donor, 'scapula_winging_axis_direction'))
    frames = sub(joint, 'frames')
    parent = None
    child = None
    for frame in donor.find('frames'):
        if frame.get('name') == 'thorax_offset':
            parent = frame
        elif frame.get('name') == 'scapula_offset':
            child = frame
    parent_rotation = rotation_xyz(numbers(parent, 'xyz_body_rotation')
                                   if parent.find('xyz_body_rotation') is not None
                                   else numbers(parent, 'orientation'))
    child_rotation = rotation_xyz(numbers(child, 'orientation'))
    if side < 0:
        flip = np.diag(LEFT_SCAPULOTHORACIC_FRAME_FLIP)
        parent_rotation = MIRROR @ parent_rotation @ flip
        child_rotation = MIRROR @ child_rotation @ flip
        for name, rotation in (('parent', parent_rotation), ('child', child_rotation)):
            if abs(np.linalg.det(rotation) - 1) > 1e-12:
                raise ValueError('Mirrored %s frame is not a rotation' % name)
    offset_frame(frames, 'torso_st_%s_offset' % suffix, '/bodyset/torso',
                 registration.thorax_to_torso(numbers(parent, 'translation'), side),
                 euler_xyz(parent_rotation))
    offset_frame(frames, 'scapula_%s_offset' % suffix, '/bodyset/scapula_%s' % suffix,
                 registration.girdle_local(numbers(child, 'translation'), side),
                 euler_xyz(child_rotation))
    return joint, ranges


def acromioclavicular_constraint(registration, side):
    suffix = 'r' if side > 0 else 'l'
    donor = None
    for element in registration.donor.find('ConstraintSet/objects'):
        if element.tag == 'PointConstraint':
            donor = element
    if donor is None:
        raise KeyError('PointConstraint AC')
    out = ET.Element('PointConstraint', name='acromioclavicular_' + suffix)
    sub(out, 'isEnforced', 'true')
    sub(out, 'socket_body_1', '/bodyset/clavicle_' + suffix)
    sub(out, 'socket_body_2', '/bodyset/scapula_' + suffix)
    sub(out, 'location_body_1',
        triple(registration.girdle_local(numbers(donor, 'location_body_1'), side)))
    sub(out, 'location_body_2',
        triple(registration.girdle_local(numbers(donor, 'location_body_2'), side)))
    return out


def spline_value_at_zero(multiplier):
    """A MultiplierFunction over a two-knot SimmSpline, evaluated at q = 0.

    The five MovingPathPoints in this donor are the three latissimus dorsi
    via-points on the thorax, each a function of the donor's `shoulder_elv`.
    This plant has no `shoulder_elv`: its shoulder is a flexion/adduction/rotation
    decomposition, not an elevation plane. So they are FROZEN at the donor's own
    value at shoulder_elv = 0, which is the arm hanging -- this model's rest pose.
    The reduction that costs is at high elevation, where the donor moves them by
    at most the spread of the two knots; the builder records that spread.
    """
    scale = float(multiplier.findtext('scale'))
    spline = multiplier.find('function/SimmSpline')
    if spline is None:
        raise ValueError('Only MultiplierFunction over SimmSpline is handled')
    x = [float(v) for v in spline.findtext('x').split()]
    y = [float(v) for v in spline.findtext('y').split()]
    if len(x) != 2 or len(y) != 2 or x[0] != 0.0:
        raise ValueError('Expected a two-knot spline starting at zero')
    return scale * y[0], scale * (y[1] - y[0])


def transfer_muscle(registration, muscle, side):
    """One donor muscle on this plant's bodies. Every donor property is kept."""
    suffix = 'r' if side > 0 else 'l'
    donor_name = muscle.get('name')
    identity = 'seth_%s_%s' % (donor_name, suffix)
    out = copy.deepcopy(muscle)
    out.set('name', identity)
    frozen = []
    path = out.find('GeometryPath')
    # The `path` property of PathActuator is serialised as <GeometryPath name="path">:
    # the element tag is the concrete type and the name attribute is the PROPERTY
    # name. The donor file is Version 40000 and calls it "geometrypath", which
    # OpenSim renames on its version upgrade -- but this model is Version 40500, so
    # no upgrade runs and a wrongly named element is silently ignored, leaving an
    # empty default path. That failure reads as "a valid path must be connected to a
    # model by at least two PathPoints" on a muscle whose XML plainly has two.
    path.set('name', 'path')
    points = path.find('PathPointSet/objects')
    replacements = []
    for index, point in enumerate(list(points), 1):
        donor_body = text_of(point, 'socket_parent_frame').rsplit('/', 1)[-1]
        target = BODY_MAP[donor_body]
        body = 'torso' if target == 'torso' else '%s_%s' % (target, suffix)
        if point.tag == 'PathPoint':
            location = numbers(point, 'location')
        elif point.tag == 'MovingPathPoint':
            values, spreads = [], []
            for axis in ('x_location', 'y_location', 'z_location'):
                value, spread = spline_value_at_zero(point.find(axis + '/MultiplierFunction'))
                values.append(value)
                spreads.append(spread)
            location = np.asarray(values, float)
            frozen.append({'donor_point': point.get('name'), 'donor_body': donor_body,
                           'frozen_at': 'shoulder_elv = 0 (arm hanging)',
                           'donor_travel_to_shoulder_elv_pi_m': spreads})
        else:
            raise ValueError('%s carries an unhandled path point %s' % (donor_name, point.tag))
        if target == 'torso':
            placed = registration.thorax_to_torso(location, side)
        elif target == 'humerus':
            placed = registration.humerus_local(location, side)
        else:
            placed = registration.girdle_local(location, side)
        replacement = ET.Element('PathPoint', name='%s-P%d' % (identity, index))
        sub(replacement, 'socket_parent_frame', '/bodyset/' + body)
        sub(replacement, 'location', triple(placed))
        replacements.append((point, replacement))
    for old, new in replacements:
        points.remove(old)
    for _, new in replacements:
        points.append(new)
    for wrap in path.iter('PathWrap'):
        wrap.find('wrap_object').text = text_of(wrap, 'wrap_object') + '_' + suffix
        wrap.set('name', identity + '_' + wrap.get('name'))
    return identity, out, frozen


def transferable(registration):
    """Donor muscles whose every point lands on a body this plant will have."""
    keep, excluded = [], []
    for muscle in registration.donor.iter('Millard2012EquilibriumMuscle'):
        bodies = sorted({text_of(p, 'socket_parent_frame').rsplit('/', 1)[-1]
                         for p in muscle.find('GeometryPath/PathPointSet/objects')})
        blocked = sorted(set(bodies) & set(DONOR_EXCLUDED_BODIES))
        if blocked:
            excluded.append({'donor_muscle': muscle.get('name'),
                             'reason': 'attaches to a donor forearm body this plant '
                                       'does not take from this donor; its elbow '
                                       'muscles would duplicate the arm26 set already '
                                       'in the plant',
                             'donor_bodies': bodies})
            continue
        unknown = sorted(set(bodies) - set(BODY_MAP))
        if unknown:
            raise ValueError('%s attaches to unmapped donor bodies %r'
                             % (muscle.get('name'), unknown))
        keep.append(muscle)
    return keep, excluded


def reparent_acromial(registration, jointset, side):
    """acromial -> scapula, with the arm left exactly where it was.

    The new parent offset frame sits at the donor's glenohumeral centre in the
    scapula (scaled) and carries the INVERSE of the scapula's reference
    rotation, so at zero arm coordinates the humerus frame is identical to the
    base plant's -- position and orientation both. That is arithmetic, and the
    verifier measures it in the engine anyway.
    """
    suffix = 'r' if side > 0 else 'l'
    name = 'acromial_' + suffix
    index, joint = next((i, j) for i, j in enumerate(jointset)
                        if j.get('name') == name)
    rotation = registration.body_rotation('scapula', side)
    translation = registration.girdle_local(registration.donor_gh_in_scapula, side)
    frames = joint.find('frames')
    for frame in list(frames):
        if frame.get('name') == 'torso_offset':
            frames.remove(frame)
    offset_frame(frames, 'scapula_%s_gh_offset' % suffix,
                 '/bodyset/scapula_%s' % suffix, translation, euler_xyz(rotation.T))
    joint.find('socket_parent_frame').text = 'scapula_%s_gh_offset' % suffix
    jointset[index] = joint
    return translation, euler_xyz(rotation.T)


def catalog_row(identity, donor_name, muscle, bodies, side, donor_sha, crosses):
    hemi = 'lh' if side > 0 else 'rh'
    return {
        'id': identity, 'side': 'r' if side > 0 else 'l', 'body_group': 'shoulder_girdle',
        'attachment_bodies': sorted(bodies),
        'sensory_region': 'brain-%s-postcentral' % hemi,
        'motor_region': 'brain-%s-precentral' % hemi,
        'max_isometric_force_n': float(text_of(muscle, 'max_isometric_force')),
        'optimal_fiber_length_m': float(text_of(muscle, 'optimal_fiber_length')),
        'tendon_slack_length_m': float(text_of(muscle, 'tendon_slack_length')),
        'source_path': DONOR,
        'source_sha256': donor_sha,
        'source_muscle_name': donor_name,
        'source_muscle_law': 'Millard2012EquilibriumMuscle',
        'crosses_joints': crosses,
        'assignment_basis': ('Contralateral regional cortical engineering prior only, '
                             'the rule every other row of this catalog uses; NOT a '
                             'claim about shoulder motor control.'),
        'parameter_basis': ('Every Millard2012EquilibriumMuscle property is the Seth '
                            '2019 donor element copied verbatim, including the curves. '
                            'Forces are the donor\'s own and are NOT scaled to this '
                            'body: the donor model states no body mass to scale from.'),
        'geometry_basis': ('Donor path points carried by the girdle map: s_lat on the '
                           'clavicle, scapula and thorax attachments, s_long and the '
                           'scapula\'s reference rotation on the humerus attachments. '
                           'Left side is the right reflected in z.'),
        'native_control_ready': False,
        'default_excitation_assignment': None,
    }


#: MEASURED, not chosen, and the measurement is in docs/SHOULDER_GIRDLE.md.
#:
#: The scapulothoracic joint's coordinates are NOT user-supplied axes: they come
#: from a Simbody ellipsoid mobilizer whose conventions are fixed inside
#: ScapulothoracicJoint.cpp. Mirroring its parent frame therefore does NOT tell
#: you what the mirrored coordinates mean, and the naive choice -- conjugate the
#: rotation by diag(1,1,-1), keep the donor's defaults -- puts the left scapula
#: 162.5 mm away from the mirror of the right and makes the left
#: acromioclavicular constraint unsatisfiable (assembly error 83 mm, where the
#: right side assembles to 1e-10 with nothing moved).
#:
#: So the convention was measured. All 48 candidates -- three ways to restore a
#: right-handed mirrored frame (right-multiply by diag(-1,1,1), diag(1,-1,1) or
#: diag(1,1,-1)) times sixteen coordinate sign patterns -- were built and loaded,
#: and each was scored on the only thing that is not a matter of taste: does the
#: LEFT girdle's rest pose come out as the exact z-mirror of the right's?
#:
#:   diag(1,-1,1) with (+1, -1, -1, +1)   origin 1.4e-16 m, rotation 1.3e-15  <-- used
#:   diag(-1,1,1) with (-1, +1, -1, -1)   identical; the same configuration, named twice
#:   the next best of the 48                origin 6.4e-03 m, rotation 3.9e-02
#:
#: Two exact answers and a 6.4 mm gap to the third: this is a discrete fact about
#: the mobilizer, not a tuned threshold. ASSEMBLY SUCCESS WAS USELESS AS A TEST --
#: 16 of the 48 assembled, including wrong ones, because the acromioclavicular
#: constraint is 3 equations on a 6-coordinate girdle and the solver simply moves
#: the coordinates off their declared defaults until it closes. A control that
#: passes for a reason unrelated to what it tests.
#:
#: WHAT THIS MEANS FOR ANY CALLER, and it is the gait2392-knee trap in advance:
#: on the LEFT, `scapula_elevation_l` and `scapula_upward_rot_l` run OPPOSITE to
#: their right-side namesakes. `scapula_abduction_l` and `scapula_winging_l` run
#: the same way. `clav_prot_l` and `clav_elev_l` run the same way as the right
#: (their joint is a CustomJoint and its axes carry the -MIRROR*a rule this
#: plant's own `acromial_l` already uses).
LEFT_SCAPULOTHORACIC_FRAME_FLIP = (1.0, -1.0, 1.0)
LEFT_SCAPULOTHORACIC_SIGNS = {'scapula_abduction': 1.0, 'scapula_elevation': -1.0,
                              'scapula_upward_rot': -1.0, 'scapula_winging': 1.0}


def build(root: Path):
    root = Path(root)
    out = root / OUT
    out.mkdir(parents=True, exist_ok=True)
    registration = Registration(root)
    tree = registration.base
    model = tree.getroot().find('Model')
    bodyset = model.find('BodySet/objects')
    jointset = model.find('JointSet/objects')
    forceset = model.find('ForceSet/objects')
    constraintset = model.find('ConstraintSet/objects')
    bodies = {body.get('name'): body for body in bodyset}

    # ---- the one-time torso debit -------------------------------------
    torso = bodies['torso']
    parent_record = {'mass_kg': float(text_of(torso, 'mass')),
                     'center_m': numbers(torso, 'mass_center').tolist(),
                     'inertia_kg_m2': inertia_matrix(numbers(torso, 'inertia')).tolist()}
    canonical_rotation = np.asarray(
        json.loads((root / CERVICAL_INERTIA).read_text())['canonical_to_current_torso'],
        float)[:3, :3]
    locals_, torso_frames, origins, rotations, inertia_evidence = {}, {}, {}, {}, {}
    for name, donor_body, side in NEW_BODIES:
        local, in_torso, origin, rotation, evidence = girdle_body_record(
            root, registration, name, donor_body, side, canonical_rotation)
        locals_[name], torso_frames[name] = local, in_torso
        origins[name], rotations[name] = origin, rotation
        inertia_evidence[name] = evidence
    residual = partition_body(parent_record, [torso_frames[n] for n, _, _ in NEW_BODIES])
    torso.find('mass').text = repr(float(residual['mass_kg']))
    torso.find('mass_center').text = triple(residual['center_m'])
    torso.find('inertia').text = inertia_text(np.asarray(residual['inertia_kg_m2'], float))

    # ---- wrap objects --------------------------------------------------
    wraps = donor_wraps(registration)
    body_wraps = {name: [] for name, _, _ in NEW_BODIES}
    for side, suffix in ((+1, 'r'), (-1, 'l')):
        for wrap in wraps.get('thorax', []):
            bodies['torso'].find('WrapObjectSet/objects').append(
                scaled_wrap(registration, wrap, 'thorax', side))
        for wrap in wraps.get('humerus', []):
            bodies['humerus_' + suffix].find('WrapObjectSet/objects').append(
                scaled_wrap(registration, wrap, 'humerus', side))
        for wrap in wraps.get('scapula', []):
            body_wraps['scapula_' + suffix].append(
                scaled_wrap(registration, wrap, 'scapula', side))
        for wrap in wraps.get('clavicle', []):
            body_wraps['clavicle_' + suffix].append(
                scaled_wrap(registration, wrap, 'clavicle', side))

    for name, _, _ in NEW_BODIES:
        bodyset.append(make_body(name, locals_[name], body_wraps[name]))

    # ---- joints, constraints, and the arm's new parent -----------------
    ranges = {}
    for side in (+1, -1):
        joint, added = sternoclavicular_joint(registration, side)
        jointset.append(joint)
        ranges.update(added)
        joint, added = scapulothoracic_joint(registration, side,
                                             LEFT_SCAPULOTHORACIC_SIGNS)
        jointset.append(joint)
        ranges.update(added)
        constraintset.append(acromioclavicular_constraint(registration, side))
    acromial = {}
    for side in (+1, -1):
        translation, orientation = reparent_acromial(registration, jointset, side)
        acromial['r' if side > 0 else 'l'] = {'translation_m': translation.tolist(),
                                              'orientation_rad': orientation.tolist()}

    for name, bounds in ranges.items():
        forceset.append(limit_force(name, bounds))
        forceset.append(damping_force(name))

    # ---- muscles --------------------------------------------------------
    keep, excluded = transferable(registration)
    donor_sha = sha(root / DONOR)
    catalog = json.loads((root / BASE_CATALOG).read_text())
    existing = {f.get('name') for f in forceset}
    added, frozen_points = [], []
    for side in (+1, -1):
        for muscle in keep:
            identity, element, frozen = transfer_muscle(registration, muscle, side)
            if identity in existing:
                raise ValueError('%s already in the model' % identity)
            forceset.append(element)
            added.append(identity)
            frozen_points.extend(frozen)
            point_bodies = sorted({text_of(p, 'socket_parent_frame').rsplit('/', 1)[-1]
                                   for p in element.find('GeometryPath/PathPointSet/objects')})
            suffix = 'r' if side > 0 else 'l'
            crosses = []
            if any(b.startswith('scapula') for b in point_bodies) and 'torso' in point_bodies:
                crosses.append('scapulothoracic_' + suffix)
            if any(b.startswith('clavicle') for b in point_bodies) and 'torso' in point_bodies:
                crosses.append('sternoclavicular_' + suffix)
            if any(b.startswith('humerus') for b in point_bodies):
                crosses.append('acromial_' + suffix)
                if 'torso' in point_bodies or any(b.startswith('clavicle') for b in point_bodies):
                    crosses.append('scapulothoracic_' + suffix)
            catalog.append(catalog_row(identity, muscle.get('name'), muscle,
                                       point_bodies, side, donor_sha,
                                       sorted(set(crosses))))

    ET.indent(tree, space='\t')
    model_path = out / 'model.osim'
    tree.write(model_path, encoding='utf-8', xml_declaration=True)
    (out / 'catalog.json').write_text(json.dumps(catalog, indent=2) + '\n')
    for name in ('initial_pose.json', 'equilibrium_excitations.json'):
        shutil.copyfile(root / BASE / name, out / name)

    base_registration = json.loads((root / BASE_REGISTRATION).read_text())
    sources = dict(base_registration['sources'])
    for relative in (BASE_MODEL, BASE_CATALOG, BASE_REGISTRATION, DONOR, DONOR_FRAMES,
                     'scripts/build_shoulder_girdle.py',
                     'scripts/measure_shoulder_donor_frames.py',
                     CERVICAL_INERTIA,
                     'ihm/assembly/cervical_inertia.py'):
        sources[relative] = sha(root / relative)
    sources[OUT + '/model.osim'] = sha(model_path)
    sources[OUT + '/catalog.json'] = sha(out / 'catalog.json')

    record = {
        'schema': 'ihm.shoulder-girdle-variant.v1',
        'model_path': OUT + '/model.osim',
        'model_sha256': sha(model_path),
        'catalog_path': OUT + '/catalog.json',
        'catalog_sha256': sha(out / 'catalog.json'),
        'base_model_path': BASE_MODEL,
        'base_model_sha256': sha(root / BASE_MODEL),
        'sources': sources,
        'muscle_count': len(catalog),
        'body_count': len(list(bodyset)),
        'default_enabled': False,
        'native_acceptance_complete': False,
        'donor': {
            'path': DONOR, 'sha256': donor_sha,
            'model': 'Thoracoscapular Shoulder Model, Seth, Dong, Matias and Delp 2019, '
                     'Front. Neurorobot. 13:90',
            'licence': 'CC BY 4.0 from the authors on SimTK (simtk.org/projects/'
                       'thoracoscapular); Apache-2.0 as a file of opensim-core, the '
                       'repository these bytes were read from. Both permit commercial '
                       'use with attribution. See docs/SHOULDER_GIRDLE.md section 1.',
            'attribution_required': 'Seth A, Dong M, Matias R, Delp SL. Muscle '
                                    'contributions to upper-extremity movement and work '
                                    'from a musculoskeletal model of the human shoulder. '
                                    'Front. Neurorobot. 13:90 (2019).',
            'not_used': 'No number in this model comes from MoBL-ARMS 4.1, whose SimTK '
                        'notice restricts it to non-commercial use.',
        },
        'map': {
            's_lat': registration.s_lat,
            's_lat_basis': ('this plant\'s glenohumeral half-width (acromial_r parent '
                            'offset on torso, z) divided by the donor\'s glenohumeral '
                            'centre in its thorax frame (z). Applied to the clavicle, '
                            'the scapula and every donor thorax attachment.'),
            's_long': registration.s_long,
            's_long_basis': ('this plant\'s glenohumeral-to-elbow distance divided by '
                             'the donor\'s, both from the models\' own joint frames. '
                             'Applied to donor humerus attachments only, because the '
                             'humerus is this plant\'s own body.'),
            'delta_torso_frame_m': registration.delta.tolist(),
            'delta_basis': ('x and y put the donor\'s glenohumeral centre exactly on '
                            'this plant\'s; z is exactly zero, so the donor midline is '
                            'this plant\'s midline and midline attachments stay on it.'),
            'our_glenohumeral_torso_m': registration.our_gh_torso.tolist(),
            'donor_glenohumeral_thorax_m': registration.donor_gh_thorax.tolist(),
            'our_humerus_length_m': registration.our_humerus_length,
            'donor_humerus_length_m': registration.donor_humerus_length,
            'donor_thorax_origin_in_torso_m':
                registration.thorax_to_torso((0.0, 0.0, 0.0)).tolist(),
            'left_scapulothoracic_signs': LEFT_SCAPULOTHORACIC_SIGNS,
            'left_scapulothoracic_frame_flip': list(LEFT_SCAPULOTHORACIC_FRAME_FLIP),
            'left_scapulothoracic_convention_basis':
                ('MEASURED over all 48 candidates by exact-mirror residual; see '
                 'scripts/build_shoulder_girdle.py and docs/SHOULDER_GIRDLE.md. '
                 'scapula_elevation_l and scapula_upward_rot_l run OPPOSITE to their '
                 'right-side namesakes.'),
        },
        'partition': {
            'instrument': 'ihm.assembly.cervical_inertia.partition_body',
            'parent_torso_xml_kg': parent_record['mass_kg'],
            'residual_torso_xml_kg': residual['mass_kg'],
            'new_bodies_xml_kg': {n: locals_[n]['mass_kg'] for n, _, _ in NEW_BODIES},
            'mass_scale_used': MASS_SCALE,
            'reconstruction_residual': residual['reconstruction_residual'],
            'donor_masses_unscaled': True,
            'inertia': inertia_evidence,
        },
        'acromial_parent_offsets': acromial,
        'new_coordinates': {k: list(v) for k, v in ranges.items()},
        'range_basis': ('EXPLICIT ENGINEERING BOUND, NOT MEASURED: the donor\'s own '
                        'declaration is +-pi for five of the six and +-pi/2 for the '
                        'sixth, which is the same placeholder as this plant\'s +-10 rad '
                        'shoulder coordinates. Each range is the donor\'s default value '
                        '+- %r rad, chosen narrow because narrow is the conservative '
                        'side for a stop. Every one carries a CoordinateLimitForce at '
                        'exactly that range in the model file.'
                        % GIRDLE_RANGE_HALFWIDTH_RAD),
        'added_muscles': added,
        'excluded_donor_muscles': excluded,
        'frozen_moving_path_points': frozen_points,
        'scope': ('articulated_spine_v1 plus a Seth 2019 shoulder girdle: clavicle and '
                  'scapula per side, sternoclavicular and scapulothoracic joints, the '
                  'acromioclavicular point constraint, the arm re-parented onto the '
                  'scapula, and 30 donor muscles per side. No linearization, no stance '
                  'acceptance, no controller.'),
        'not_claimed': [
            'No linearization.npz. The identified stance plant was solved on 22 bodies '
            'and 33 coordinates; this one has 29 and 60.',
            'The girdle coordinate ranges are an engineering bound with no measured '
            'provenance. No donor on disk declares a scapular or clavicular range.',
            'The left side is the right reflected in z. The donor has one arm.',
            'Donor muscle forces are unscaled: the donor model states no body mass.',
            'The clavicle and scapula inertia tensors are NOT the donor\'s: the '
            'donor\'s are not physically realisable, so they are convex-envelope '
            'priors of this body\'s own anatomy at the donor\'s stated mass. The '
            'mass centres ARE the donor\'s, scaled.',
            'The three latissimus dorsi via-points are moving points in the donor, a '
            'function of its shoulder_elv; they are frozen at shoulder_elv = 0 here.',
            'Display geometry is unchanged: the scapulae and clavicles still ride the '
            'thorax mesh in data/derived/anatomy-segment-binding.',
            'Capacity measured from this model is Sum(Fmax * moment arm) at one pose. '
            'It is an upper bound on torque, and it is not demonstrated behaviour.',
        ],
    }
    (out / 'registration.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def main():
    record = build(ROOT)
    print(json.dumps({k: record[k] for k in ('model_path', 'model_sha256', 'muscle_count',
                                             'body_count')}, indent=2))
    print(json.dumps(record['map'], indent=2))
    print(json.dumps(record['partition'], indent=2))
    print('added muscles:', len(record['added_muscles']))
    print('excluded:', [e['donor_muscle'] for e in record['excluded_donor_muscles']])


if __name__ == '__main__':
    main()
