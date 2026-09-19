"""Can any muscle already in the plant be reassigned so that it drives `thoracic_*`?

    OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
        .venv/bin/python -m scripts.measure_thoracic_reassignability

WHY THIS EXISTS.  docs/ARTICULATED_SPINE.md records the thoracic joint as
BLOCKED, topology class, for two stated reasons:

    "The only candidates in the plant are the gait2392 trunk muscles.  Their
     `torso` insertion (y = 0.11 m) sits 32 mm above the thoracic joint centre,
     which is itself inherited from a registration with 62 mm RMS proxy
     residual, so reassigning them to `thorax` is not decidable from this data."

That verdict was written on `articulated_spine_v1`, whose only above-pelvis
muscles were the 12 arm26 and 6 gait2392 trunk elements.  `shoulder_girdle_v1`
has since added SIXTY Seth 2019 donor muscles, thirty per side, of which 32 take
an attachment on `torso`.  The verdict has to be re-asked against those, and
that is all this script does.  It does not build anything and it loads no engine.

WHAT IS AND IS NOT AT ISSUE.  The `thoracic` joint is `torso -> thorax`.  Its
parent frame sits on `torso` at a pure translation with ZERO orientation, and
its child frame is the `thorax` origin, so at `thoracic_* = 0` a point at `p` in
torso coordinates is at `p - t` in thorax coordinates and NOTHING ELSE CHANGES.
Reassigning an attachment from `torso` to `thorax` is therefore exact arithmetic
with no registration and no fitted quantity.  The only question is whether the
attachment BELONGS on the thorax, and that is an anatomy question, not an
arithmetic one.

THE TRAP THIS IS WRITTEN AGAINST.  IHM-1 CLAUDE.md: "A hand-placed coordinate is
not a measurement, even in the right frame ... Check each ENDPOINT against the
structure its own label names."  A muscle called "serratus anterior" is not
evidence that its transferred attachment point lands on a rib.  So the rule below
decides membership from the point's distance to the ACTUAL SURFACE of this body's
own anatomical meshes, and the muscle's name is never consulted.

======================================================================
PRE-REGISTRATION.  Written and committed to this file before the script was
first run, and not edited afterwards.
======================================================================

`thorax` is not "the upper trunk".  It is exactly the 48 structures that
`data/research/thoracic_mechanism/native_composition_v1/plan.json` debited from
the torso: ribs 1-12, costal cartilages 1-7, manubrium, body of sternum, xiphoid
process, the three intercostal layers and the diaphragm, both sides.  It holds NO
thoracic vertebra -- T1-T12 stayed in the 17.938 kg residual core inside `torso`
-- so a spinous-process origin does NOT belong on `thorax` however high it sits.

An attachment currently on `torso` is REASSIGNABLE to `thorax` iff BOTH hold.

  R1  ENDPOINT.  Its nearest anatomical surface, by point-to-triangle-vertex
      distance in the canonical frame over the whole skeletal candidate set
      (every thorax-partition member, every cervical/thoracic/lumbar vertebra,
      the sacrum, both hip bones, both clavicles and both scapulae), is a member
      of the thorax partition.  Nearest wins outright; there is no tolerance,
      because a tolerance here would be a number with no provenance.

  R2  JOINT-CENTRE MARGIN.  Its height above the thoracic joint centre in the
      torso frame exceeds 87.7 mm = sqrt(2) x 62 mm.

      62 mm is the declared RMS proxy residual of the registration that PLACED
      that joint centre (`data/research/cervical_inertia/v2/manifest.json`, and
      docs/SHOULDER_GIRDLE.md section 9 measures the same order for the girdle's
      own placement).  Both sides of this comparison are measured with an error
      of that order -- the joint centre and the attachment -- so the margin has
      to clear about sqrt(2) x sd and not a bare sd.  That is IBM-1 CLAUDE.md's
      standing rule, "A margin over a measured baseline must clear sampling error
      on BOTH sides", applied to a distance instead of a skill score.

A MUSCLE is reassignable iff EVERY ONE of its torso attachments is.  A muscle
with some points on ribs and some on vertebrae is NOT reassigned and is reported
as MIXED: splitting a path across two bodies at a via point is a different and
much larger claim than moving a whole attachment.

VERDICT RULE, fixed here in advance:
  * 0 reassignable muscles  -> the thoracic joint stays BLOCKED and the
    ARTICULATED_SPINE verdict stands, now with a measurement behind it.
  * >= 1 per side           -> the thoracic joint is UNBLOCKED for those muscles
    and `scripts/build_thoracic_drive.py` may transfer exactly that set.

Anything this script reports that the rule does not cover is reported and NOT
acted on.

OUTPUT: out/thoracic_reassignability.json
"""
from __future__ import annotations

import gzip
import json
import pathlib
import sys
import xml.etree.ElementTree as ET

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]

MODEL = ROOT / 'data/models/shoulder_girdle_v1/model.osim'
PLAN = ROOT / 'data/research/thoracic_mechanism/native_composition_v1/plan.json'
CERVICAL = ROOT / 'data/research/cervical_inertia/v2/manifest.json'
ANATOMY = ROOT / 'data/derived/canonical/anatomy.json'
GEOMETRY = ROOT / 'data/derived/canonical/geometry'
OUT = ROOT / 'out/thoracic_reassignability.json'

#: sqrt(2) x the 62 mm RMS proxy residual declared for the registration that
#: placed the thoracic joint centre.  See the pre-registration above.
JOINT_CENTRE_RESIDUAL_M = 0.062
MARGIN_BAR_M = float(np.sqrt(2.0) * JOINT_CENTRE_RESIDUAL_M)

#: Every vertebra, the sacrum, both hip bones, both clavicles and both scapulae.
#: These are the structures a trunk or girdle attachment could plausibly be ON
#: that are NOT part of the thorax.  Named by FMA-derived label, matched against
#: `anatomy.json` `name`; the ids are resolved at run time and printed, so a
#: label that matches nothing is visible rather than silently dropped.
NON_THORAX_LABEL_SUBSTRINGS = (
    'thoracic vertebra', 'cervical vertebra', 'lumbar vertebra',
    'atlas', 'axis', 'sacrum', 'coccyx',
    'hip bone', 'clavicle', 'scapula',
)


def joint_frames(root: ET.Element, joint: str):
    for element in root.iter():
        if element.tag.endswith('Joint') and element.get('name') == joint:
            frames = {}
            for frame in element.iter('PhysicalOffsetFrame'):
                frames[frame.findtext('socket_parent')] = (
                    np.array([float(x) for x in frame.findtext('translation').split()]),
                    np.array([float(x) for x in frame.findtext('orientation').split()]))
            return frames
    raise SystemExit('joint %s not found' % joint)


def torso_attachments(root: ET.Element):
    """Every muscle path point whose parent frame is the torso body itself."""
    rows = []
    for muscle in root.iter():
        name = muscle.get('name')
        if not name or not muscle.tag.endswith('Muscle'):
            continue
        points = []
        moving = 0
        for point in muscle.iter('PathPoint'):
            body = point.findtext('socket_parent_frame')
            if body and body.rstrip('/').endswith('torso'):
                points.append((point.get('name'),
                               np.array([float(x) for x in point.findtext('location').split()])))
        for point in muscle.iter('MovingPathPoint'):
            body = point.findtext('socket_parent_frame')
            if body and body.rstrip('/').endswith('torso'):
                moving += 1
        if points or moving:
            rows.append({'muscle': name, 'law': muscle.tag.split('}')[-1],
                         'points': points, 'moving_points_on_torso': moving})
    return rows


def surface(entity_id: str):
    path = GEOMETRY / ('%s.json.gz' % entity_id)
    if not path.exists():
        return None
    data = json.loads(gzip.open(path).read())
    return np.asarray(data['positions'], float).reshape(-1, 3)


def main() -> int:
    root = ET.parse(MODEL).getroot()

    frames = joint_frames(root, 'thoracic')
    torso_frame = [v for k, v in frames.items() if k.endswith('torso')][0]
    thorax_frame = [v for k, v in frames.items() if k.endswith('thorax')][0]
    centre, centre_rot = torso_frame
    if float(np.abs(centre_rot).max()) != 0.0 or float(np.abs(thorax_frame[0]).max()) != 0.0:
        raise SystemExit('the thoracic joint is not the pure translation this '
                         'script assumes; re-derive the offset before trusting it')

    plan = json.loads(PLAN.read_text())
    thorax_ids = set(plan['thoracic_material_partitions'])

    cervical = json.loads(CERVICAL.read_text())
    canonical_to_torso = np.asarray(cervical['canonical_to_current_torso'], float)
    torso_to_canonical = np.linalg.inv(canonical_to_torso)

    anatomy = json.loads(ANATOMY.read_text())
    by_id = {e['id']: e for e in anatomy['entities']}

    candidates = {}
    for entity_id in sorted(thorax_ids):
        if entity_id in by_id:
            candidates[entity_id] = ('thorax', by_id[entity_id]['name'])
    resolved_non_thorax = []
    for entity_id, entity in by_id.items():
        if entity_id in candidates:
            continue
        label = (entity.get('name') or '').lower()
        if any(s in label for s in NON_THORAX_LABEL_SUBSTRINGS):
            candidates[entity_id] = ('other', entity['name'])
            resolved_non_thorax.append(entity['name'])

    meshes = {}
    missing = []
    for entity_id in candidates:
        points = surface(entity_id)
        if points is None or not len(points):
            missing.append(entity_id)
            continue
        meshes[entity_id] = points

    rows = torso_attachments(root)
    report_muscles = []
    for row in rows:
        detail = []
        for point_name, location in row['points']:
            canonical = (torso_to_canonical @ np.append(location, 1.0))[:3]
            best_id, best_d = None, np.inf
            for entity_id, vertices in meshes.items():
                d = float(np.min(np.linalg.norm(vertices - canonical, axis=1)))
                if d < best_d:
                    best_id, best_d = entity_id, d
            group, label = candidates[best_id]
            margin = float(location[1] - centre[1])
            detail.append({
                'point': point_name,
                'location_torso_m': [float(x) for x in location],
                'height_above_thoracic_joint_centre_m': margin,
                'nearest_structure': label,
                'nearest_structure_id': best_id,
                'nearest_surface_distance_m': best_d,
                'R1_endpoint_is_a_thorax_structure': group == 'thorax',
                'R2_margin_clears_sqrt2_x_62mm': margin > MARGIN_BAR_M,
            })
        verdicts = [(d['R1_endpoint_is_a_thorax_structure'] and
                     d['R2_margin_clears_sqrt2_x_62mm']) for d in detail]
        if row['moving_points_on_torso']:
            verdict = 'EXCLUDED_MOVING_POINT'
        elif all(verdicts) and verdicts:
            verdict = 'REASSIGNABLE'
        elif any(verdicts):
            verdict = 'MIXED'
        else:
            verdict = 'STAYS_ON_TORSO'
        report_muscles.append({'muscle': row['muscle'], 'law': row['law'],
                               'verdict': verdict,
                               'moving_points_on_torso': row['moving_points_on_torso'],
                               'points': detail})

    reassignable = sorted(m['muscle'] for m in report_muscles
                          if m['verdict'] == 'REASSIGNABLE')
    per_side = {'r': len([m for m in reassignable if m.endswith('_r')]),
                'l': len([m for m in reassignable if m.endswith('_l')])}
    verdict = ('UNBLOCKED' if per_side['r'] >= 1 and per_side['l'] >= 1 else 'BLOCKED')

    report = {
        'schema': 'ihm.thoracic-reassignability.v1',
        'model': str(MODEL.relative_to(ROOT)),
        'thoracic_joint_centre_torso_m': [float(x) for x in centre],
        'margin_bar_m': MARGIN_BAR_M,
        'margin_bar_basis': 'sqrt(2) x the 62 mm RMS proxy residual declared for '
                            'the registration that placed the thoracic joint centre',
        'thorax_partition_structures': len(thorax_ids),
        'thorax_partition_with_surface': len([i for i in meshes if candidates[i][0] == 'thorax']),
        'non_thorax_candidates_resolved': sorted(resolved_non_thorax),
        'candidate_surfaces_missing': sorted(missing),
        'muscles_with_a_torso_attachment': len(report_muscles),
        'reassignable': reassignable,
        'reassignable_per_side': per_side,
        'verdict': verdict,
        'muscles': report_muscles,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + '\n')

    print('thoracic joint centre in the torso frame: %.5f %.5f %.5f m' % tuple(centre))
    print('margin bar: %.4f m (sqrt(2) x %.3f m)' % (MARGIN_BAR_M, JOINT_CENTRE_RESIDUAL_M))
    print('thorax partition structures with a surface on disk: %d of %d'
          % (report['thorax_partition_with_surface'], len(thorax_ids)))
    print('non-thorax candidate structures resolved: %d' % len(resolved_non_thorax))
    print()
    header = '%-34s %-22s %9s %9s %s'
    print(header % ('muscle', 'nearest structure', 'dist_mm', 'margin_mm', 'verdict'))
    for muscle in sorted(report_muscles, key=lambda m: m['muscle']):
        first = muscle['points'][0] if muscle['points'] else None
        if first is None:
            print(header % (muscle['muscle'], '(moving point only)', '', '', muscle['verdict']))
            continue
        worst = min(muscle['points'], key=lambda p: p['height_above_thoracic_joint_centre_m'])
        print(header % (muscle['muscle'], worst['nearest_structure'][:22],
                        '%.1f' % (1000 * worst['nearest_surface_distance_m']),
                        '%.1f' % (1000 * worst['height_above_thoracic_joint_centre_m']),
                        muscle['verdict']))
    print()
    print('REASSIGNABLE: %d (%d right, %d left)'
          % (len(reassignable), per_side['r'], per_side['l']))
    print('VERDICT: the thoracic joint is %s' % verdict)
    print('report: %s' % OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
