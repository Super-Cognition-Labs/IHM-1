"""Put the donor's own neck muscles on the articulated-spine variant's neck.

    .venv/bin/python -m scripts.transfer_neck_muscle_paths

`data/models/articulated_spine_v1` gave the body a thorax, a lumped cervical
segment and a head, and measured that no muscle crosses any of the joints between
them.  `scripts/refit_muscle_paths_articulated_spine.py` proves the zero is
TOPOLOGICAL there: the model's own GeometryPaths give zero too, because no muscle
in the plant has a path point above the torso.  A refit cannot help.  A muscle
has to be added.

THE DONOR IS ALREADY HERE AND ALREADY REGISTERED.  Nothing below is authored.

  data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim
      78 Thelen2003Muscle neck muscles, all plain PathPoints, no wrap objects.
  data/research/cervical_registration/v2/recipe.json
      the rigid donor-spine -> target-torso fit (rms 14.9 mm over eight
      bone-envelope centroids -- engineering proxies, not homologous landmarks;
      the recipe says so), and the resulting reference transform of EVERY donor
      body into this torso frame.  The variant's `neck` and `atlantooccipital`
      joint centres were placed from the same records and coincide with the
      recipe's `cerv7` and `skull` origins to every printed digit; the verifier
      checks that to 1e-12 m rather than trusting this sentence.

THE BODY MAP is a KINEMATIC correspondence, joint for joint, not a height rule:

  donor spine, torso (below auxt1jnt, T1-C7)   -> thorax    (below `neck`)
  donor cerv1 ... cerv7                        -> cervical  (between the two)
  donor skull, jaw (above aux1jnt, C1-skull)   -> head      (above `atlantooccipital`)

The variant's `neck` sits at the donor's auxt1jnt (C7/T1) and its
`atlantooccipital` at the donor's skull origin, so every donor joint is either one
of ours or internal to our lumped `cervical`.  Every donor `spine`/`torso` point a
transferred muscle uses is also checked to lie ABOVE the variant's `thoracic`
joint centre; one below it would be ambiguous between `thorax` and `torso`, and
the transfer refuses rather than guesses.

STATIONS.  p_torso = T_donor_body->torso . p_donor, from the recipe; then
p_ours = p_torso - origin_ours, because the three new bodies are axis-aligned
with `torso` and have their origins at their proximal joint centres
(`build_articulated_spine.make_joint`).  One rigid transform per donor body, all
derived from ONE rigid fit, so at the reference pose every path length is the
donor's own to floating point -- the known answer the verifier checks against the
donor model evaluated by OpenSim.

WHICH MUSCLES, by rule, not by choice:

  transferred   every donor muscle whose points map to at least two of our
                bodies (so it crosses at least one of our joints) and touch no
                excluded body.
  excluded, GIRDLE     any point on rclavicle/lclavicle/rscapula/lscapula:
                cleid_mast, cleid_occ, trap_cl, trap_acr, levator_scap.  This
                plant has no clavicle or scapula body; the girdle mass rides
                `torso` (docs/ARTICULATED_SPINE.md).  Putting these on `torso`
                would make the cleidomastoid cross the thoracic joint, which it
                does not.  Blocked on the same girdle as the shoulder.
  excluded, INTERNAL   every point on ONE of our bodies -- the short
                intersegmental muscles (multifidi C-C, obl_cap_inf,
                long_col_c1c5, deepmult-T2-T1).  In a lumped cervical body they
                cross no joint, their length is constant, and they would add
                muscle states with no mechanical effect.  This is a consequence
                of lumping C1-C7, which the variant already states.

PARAMETERS.  max_isometric_force, optimal_fiber_length, tendon_slack_length and
pennation_angle_at_optimal are the donor's own text, copied verbatim.  The donor
declares nothing else for its Thelen2003Muscles and has no <defaults> block, so
every other Thelen property is OpenSim's class default -- which is what the donor
itself runs with when this engine's OpenSim loads it.  That is stated, not
authored: no value here was chosen.  The donor is a 50th-percentile male model and
its forces are NOT scaled to this body, exactly as the gait2392 trunk muscles in
this plant were not (`catalog.json`, "Unscaled retained donor parameters").

OUTPUTS, all in data/models/articulated_spine_v1/, none of which touch the files
the recorded gates were measured on:

  model_muscled.osim, catalog_muscled.json, registration_muscled.json
      model.osim plus the transferred muscles; with `muscle_paths.xml` as the
      FunctionBasedPathSet source override when it is installed.
  registration_foot_paths.json
      model.osim UNCHANGED plus only `muscle_paths.xml` -- the plant that
      attributes any subtalar/ankle change to the foot paths alone.

`muscle_paths.xml` is NOT a refit.  It is the shipped path set with the 22
muscles whose geometry crosses subtalar or mtp removed from it, so the engine
runs those 22 on the model's own GeometryPaths.  Two refits were pre-registered
and both FAILED their no-regression gate; see
`scripts/refit_muscle_paths_articulated_spine.py` and docs/ARTICULATED_SPINE.md.

The transfer is a function of its inputs; the verifier runs it twice and compares
bytes.
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

OUT = 'data/models/articulated_spine_v1'
MODEL = OUT + '/model.osim'
CATALOG = OUT + '/catalog.json'
REGISTRATION = OUT + '/registration.json'
PATHS = OUT + '/muscle_paths.xml'
MUSCLED_MODEL = OUT + '/model_muscled.osim'
MUSCLED_CATALOG = OUT + '/catalog_muscled.json'
MUSCLED_REGISTRATION = OUT + '/registration_muscled.json'
PATHS_REGISTRATION = OUT + '/registration_foot_paths.json'
MASI = 'data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim'
RECIPE = 'data/research/cervical_registration/v2/recipe.json'
PATHSET_KEY = 'subject_walk_scaled_FunctionBasedPathSet.xml'

BODY_MAP = {'spine': 'thorax', 'torso': 'thorax', 'skull': 'head', 'jaw': 'head',
            **{'cerv%d' % level: 'cervical' for level in range(1, 8)}}
GIRDLE = ('rclavicle', 'lclavicle', 'rscapula', 'lscapula')
#: Which joint each of our new bodies' ORIGIN sits on (build_articulated_spine).
ORIGIN_JOINT = {'thorax': 'thoracic', 'cervical': 'neck', 'head': 'atlantooccipital'}
#: The only four properties the donor declares for its Thelen2003Muscles.
DONOR_PROPERTIES = ('max_isometric_force', 'optimal_fiber_length',
                    'tendon_slack_length', 'pennation_angle_at_optimal')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def donor_transforms(recipe):
    """{donor body: 4x4 donor-frame -> target-torso-frame at the reference pose}."""
    out = {name: np.asarray(record['body_to_target_torso_reference'], float)
           for name, record in recipe['bodies'].items()}
    for name in ('spine', 'torso'):
        out[name] = np.asarray(recipe['fixed_anchor_frames'][name]['source_to_target_torso'], float)
    for name, transform in out.items():
        rotation = transform[:3, :3]
        if (np.abs(rotation @ rotation.T - np.eye(3)).max() > 1e-9
                or abs(np.linalg.det(rotation) - 1) > 1e-9):
            raise ValueError('Recipe transform for %s is not a proper rotation' % name)
    return out


def our_name(donor):
    """masi_<donor name>_<side>.  The donor names its right muscles bare and its
    left ones with `_l`; characters OpenSim will not take in a component name
    ('/', '-') become '_'.  The donor name is kept in the catalog verbatim."""
    side = 'l' if donor.endswith('_l') else 'r'
    stem = donor[:-2] if side == 'l' else donor
    return 'masi_%s_%s' % (re.sub(r'[^A-Za-z0-9_]', '_', stem), side), side


#: The lumped chain, proximal to distal, and the joint between each neighbour.
CHAIN = ('thorax', 'cervical', 'head')
CHAIN_JOINTS = ('neck', 'atlantooccipital')


def crossed(bodies):
    """Joints of the lumped chain a path spanning these bodies crosses."""
    index = sorted(CHAIN.index(b) for b in set(bodies))
    return list(CHAIN_JOINTS[index[0]:index[-1]])


def plan(root):
    """Which donor muscles transfer, where each point lands, and why the rest do not."""
    recipe = json.loads((root / RECIPE).read_text())
    registration = json.loads((root / REGISTRATION).read_text())
    centres = {k: np.asarray(v, float)
               for k, v in registration['joint_centres_torso_frame_m'].items()}
    origin = {body: centres[joint] for body, joint in ORIGIN_JOINT.items()}
    transforms = donor_transforms(recipe)
    donor = ET.parse(root / MASI).getroot()

    transferred, excluded = [], []
    for muscle in donor.iter('Thelen2003Muscle'):
        name = muscle.get('name')
        points = [p for p in muscle.iter() if p.tag.endswith('PathPoint')]
        if any(p.tag != 'PathPoint' for p in points) or muscle.find('.//PathWrap') is not None:
            raise ValueError('%s: only plain PathPoints without wraps are handled' % name)
        bodies = [p.findtext('body').strip() for p in points]
        if set(bodies) & set(GIRDLE):
            excluded.append({'donor_muscle': name, 'reason': 'girdle',
                             'donor_bodies': sorted(set(bodies))})
            continue
        unknown = set(bodies) - set(BODY_MAP)
        if unknown:
            raise ValueError('%s attaches to unmapped donor bodies %r' % (name, sorted(unknown)))
        ours = [BODY_MAP[b] for b in bodies]
        if len(set(ours)) < 2:
            excluded.append({'donor_muscle': name, 'reason': 'internal_to_one_lumped_body',
                             'our_body': ours[0], 'donor_bodies': sorted(set(bodies))})
            continue
        stations = []
        for point, body, target in zip(points, bodies, ours):
            local = np.asarray([float(v) for v in point.findtext('location').split()])
            in_torso = (transforms[body] @ np.r_[local, 1.0])[:3]
            if target == 'thorax' and in_torso[1] <= centres['thoracic'][1]:
                raise ValueError('%s point %s lies below the thoracic joint; thorax/torso '
                                 'is ambiguous and the transfer refuses to guess'
                                 % (name, point.get('name')))
            stations.append({'donor_point': point.get('name'), 'donor_body': body,
                             'donor_station_m': local.tolist(), 'body': target,
                             'torso_frame_m': in_torso.tolist(),
                             'target_station_m': (in_torso - origin[target]).tolist()})
        identity, side = our_name(name)
        transferred.append({
            'id': identity, 'side': side, 'donor_muscle': name,
            'properties': {k: muscle.findtext(k).strip() for k in DONOR_PROPERTIES},
            'stations': stations,
            'crosses': crossed(ours)})
    if len({m['id'] for m in transferred}) != len(transferred):
        raise ValueError('Transferred muscle names collide')
    return transferred, excluded


def muscle_element(record):
    muscle = ET.Element('Thelen2003Muscle', name=record['id'])
    ET.SubElement(muscle, 'appliesForce').text = 'true'
    path = ET.SubElement(muscle, 'GeometryPath', name='path')
    point_set = ET.SubElement(path, 'PathPointSet')
    points = ET.SubElement(point_set, 'objects')
    for index, station in enumerate(record['stations'], 1):
        point = ET.SubElement(points, 'PathPoint', name='%s-P%d' % (record['id'], index))
        ET.SubElement(point, 'socket_parent_frame').text = '/bodyset/' + station['body']
        ET.SubElement(point, 'location').text = ' '.join(repr(float(v)) for v in
                                                         station['target_station_m'])
    ET.SubElement(point_set, 'groups')
    wraps = ET.SubElement(path, 'PathWrapSet')
    ET.SubElement(wraps, 'objects')
    ET.SubElement(wraps, 'groups')
    for key in DONOR_PROPERTIES:
        ET.SubElement(muscle, key).text = record['properties'][key]
    return muscle


def catalog_row(record, donor_sha):
    bodies = sorted({s['body'] for s in record['stations']})
    hemi = 'lh' if record['side'] == 'r' else 'rh'
    return {
        'id': record['id'], 'side': record['side'], 'body_group': 'neck',
        'attachment_bodies': bodies,
        'sensory_region': 'brain-%s-postcentral' % hemi,
        'motor_region': 'brain-%s-precentral' % hemi,
        'max_isometric_force_n': float(record['properties']['max_isometric_force']),
        'optimal_fiber_length_m': float(record['properties']['optimal_fiber_length']),
        'tendon_slack_length_m': float(record['properties']['tendon_slack_length']),
        'source_path': MASI, 'source_sha256': donor_sha,
        'source_muscle_name': record['donor_muscle'],
        'source_muscle_law': 'Thelen2003Muscle',
        'crosses_joints': record['crosses'],
        'path_points': [{k: s[k] for k in ('donor_point', 'donor_body', 'donor_station_m',
                                            'body', 'target_station_m')}
                        for s in record['stations']],
        'assignment_basis': ('Contralateral regional cortical engineering prior only, the '
                             'rule every other row of this catalog uses; neck motor '
                             'control is substantially bilateral and this is NOT a '
                             'claim about it.'),
        'parameter_basis': ('Donor text verbatim for max_isometric_force, '
                            'optimal_fiber_length, tendon_slack_length, '
                            'pennation_angle_at_optimal; every other Thelen2003Muscle '
                            'property is OpenSim\'s class default because the donor '
                            'declares none. Unscaled 50th-percentile-male donor forces.'),
        'geometry_basis': ('Donor PathPoints moved by the cervical registration recipe\'s '
                           'rigid reference transform (one fit, rms 14.9 mm over '
                           'bone-envelope centroid proxies) and re-expressed in the '
                           'lumped variant body. Exact at the reference pose; the lumped '
                           'kinematics differ from the donor\'s seven-level chain under '
                           'motion.'),
        'native_control_ready': False,
        'default_excitation_assignment': None,
    }


def sources_for(root, extra):
    base = json.loads((root / REGISTRATION).read_text())['sources']
    out = dict(base)
    for relative in extra:
        out[relative] = sha(root / relative)
    return out


def build(root):
    root = Path(root)
    transferred, excluded = plan(root)
    tree = ET.parse(root / MODEL)
    forces = tree.getroot().find('.//ForceSet/objects')
    existing = {f.get('name') for f in forces}
    for record in transferred:
        if record['id'] in existing:
            raise ValueError('%s already in the model' % record['id'])
        forces.append(muscle_element(record))
    ET.indent(tree, space='\t')
    tree.write(root / MUSCLED_MODEL, encoding='utf-8', xml_declaration=True)

    donor_sha = sha(root / MASI)
    catalog = json.loads((root / CATALOG).read_text())
    catalog += [catalog_row(r, donor_sha) for r in transferred]
    (root / MUSCLED_CATALOG).write_text(json.dumps(catalog, indent=2) + '\n')

    base = json.loads((root / REGISTRATION).read_text())
    paths = root / PATHS
    override = ({PATHSET_KEY: {'path': PATHS, 'sha256': sha(paths)}}
                if paths.exists() else {})
    inputs = [MASI, RECIPE, 'scripts/transfer_neck_muscle_paths.py', REGISTRATION,
              MODEL, CATALOG] + ([PATHS] if override else [])
    common = {k: base[k] for k in ('base_model_path', 'base_model_sha256', 'body_count',
                                   'mass_scale_used', 'joint_centres_torso_frame_m',
                                   'partition_kg_mass_scaled')}

    muscled = dict(
        schema='ihm.articulated-spine-variant.v1',
        stage='muscled',
        model_path=MUSCLED_MODEL, model_sha256=sha(root / MUSCLED_MODEL),
        catalog_path=MUSCLED_CATALOG, catalog_sha256=sha(root / MUSCLED_CATALOG),
        kinematic_registration=REGISTRATION,
        **common,
        muscle_count=len(catalog),
        added_muscles=[r['id'] for r in transferred],
        source_overrides=override,
        sources={**sources_for(root, inputs),
                 MUSCLED_MODEL: sha(root / MUSCLED_MODEL),
                 MUSCLED_CATALOG: sha(root / MUSCLED_CATALOG)},
        default_enabled=False, native_acceptance_complete=False,
        transfer={'donor': MASI, 'donor_sha256': donor_sha, 'recipe': RECIPE,
                  'recipe_sha256': sha(root / RECIPE), 'body_map': BODY_MAP,
                  'transferred': len(transferred),
                  'excluded': excluded},
        scope=('articulated_spine_v1 plus the donor MASI neck muscles that cross a '
               'joint of the lumped neck, and -- when installed -- the 22 foot muscles '
               'crossing subtalar or mtp on their own GeometryPaths. No '
               'linearization, no stance acceptance.'),
        not_claimed=[
            'No muscle crosses thoracic_* or either wrist: those remain free '
            'coordinates with a passive stop and viscous damping only.',
            'The neck muscles are placed exactly at the reference pose; under motion '
            'the lumped two-joint neck is not the donor\'s seven-level coupled chain, '
            'so arms differ from the donor\'s and are reported beside them.',
            'Clavicle- and scapula-anchored neck muscles are absent (girdle).',
            'Forces are the donor\'s unscaled 50th-percentile-male values.',
        ])
    (root / MUSCLED_REGISTRATION).write_text(json.dumps(muscled, indent=2) + '\n')

    written = [MUSCLED_REGISTRATION]
    if override:
        refit_only = dict(base)
        refit_only.update(
            stage='foot_geometry_paths_only',
            kinematic_registration=REGISTRATION,
            source_overrides=override,
            sources={**base['sources'], PATHS: sha(paths),
                     'scripts/transfer_neck_muscle_paths.py':
                         sha(root / 'scripts/transfer_neck_muscle_paths.py')},
            scope=('articulated_spine_v1 with its model file UNCHANGED and only the '
                   'FunctionBasedPathSet overridden: the 22 muscles whose geometry '
                   'crosses subtalar or mtp run on their own GeometryPaths, the other '
                   '58 on the shipped coefficients. The plant that attributes any '
                   'subtalar/ankle change to the foot paths alone.'))
        refit_only['not_claimed'] = [
            c for c in base['not_claimed']] + [
            'Subtalar and mtp now carry the GeometryPath moment arms of the 22 '
            'muscles that cross them; the spine, neck and wrist coordinates still '
            'carry none.']
        (root / PATHS_REGISTRATION).write_text(json.dumps(refit_only, indent=2) + '\n')
        written.append(PATHS_REGISTRATION)
    return {'transferred': [r['id'] for r in transferred], 'excluded': excluded,
            'written': written, 'override': bool(override)}


if __name__ == '__main__':
    result = build(ROOT)
    print(json.dumps({'transferred': len(result['transferred']),
                      'excluded': [(e['donor_muscle'], e['reason']) for e in result['excluded']],
                      'written': result['written'], 'override': result['override']},
                     indent=2))
