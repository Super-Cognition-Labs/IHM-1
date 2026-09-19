#!/usr/bin/env python3
"""Known answers for the live plant's display pose.

`ArticulatedBodyPlant(display_pose=...)` carries the anatomy's pose beside the
plant's own material embedding. The pose itself is already verified against
Simbody by `scripts/verify_anatomy_pose.py`; what is NOT covered there is the
thing this file adds -- that the compact 22-segment form a frame actually ships
reconstructs the per-entity pose EXACTLY, and that the field means what its basis
says it means.

Every check is a known answer: a case whose result is fixed before it runs.

    .venv/bin/python scripts/verify_display_pose.py
"""
import sys, uuid
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
REGISTRATION = 'data/derived/mechanics/whole_body_arm26_v2/registration.json'
MASS = 77.6122029


def check(name, ok, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {name}{"  " + detail if detail else ""}')
    return bool(ok)


def main() -> int:
    from ihm.assembly.articulated import ArticulatedBodyPlant
    out = ROOT / f'data/derived/scratch-verifydisp-{uuid.uuid4().hex[:8]}'
    plant = ArticulatedBodyPlant(ROOT, out, environment='upright', target_mass_kg=MASS,
                                 augmented_registration=REGISTRATION, display_pose='opensim')
    ok = True
    try:
        poser = plant.display_poser
        native = plant.native.snapshot()
        frame = plant.state
        dp = frame['display_pose']
        m = plant.display_pose_map

        print('reconstruction')
        # the whole reason the frame is small: applying the segment motion an entity is
        # mapped to must give back the per-entity pose, to the float, not approximately.
        pose = poser.pose_from_native(native)
        M = np.array(dp['segment_motion'])
        seg = np.array(m['segment_of'])
        rest = np.array(m['rest_centroid_m'])
        rebuilt = np.einsum('nij,nj->ni', M[seg][:, :3, :3], rest) + M[seg][:, :3, 3]
        direct = np.einsum('nij,nj->ni', pose.rotation, rest) + pose.translation
        gap = np.abs(rebuilt - direct).max()
        ok &= check('22 segment motions reconstruct 3,995 entity poses', gap == 0.0,
                    f'max |difference| = {gap:.3e} m')
        ok &= check('map covers every posed entity',
                    len(m['entity_ids']) == len(poser.entity_ids) == len(seg) == len(rest),
                    f'{len(m["entity_ids"])} entities')
        ok &= check('segment count is the model\'s own', len(dp['segments']) == 22,
                    f'{len(dp["segments"])} segments')

        print('idempotence -- the CLAUDE.md check: call it twice at the same input')
        a = plant._display_pose(native)
        b = plant._display_pose(native)
        ok &= check('same state twice gives bit-identical motions',
                    a['segment_motion'] == b['segment_motion'])

        print('rigidity')
        det = np.linalg.det(M[:, :3, :3])
        orth = np.abs(np.einsum('nij,nkj->nik', M[:, :3, :3], M[:, :3, :3]) - np.eye(3)).max()
        ok &= check('every segment motion is a proper rotation',
                    np.allclose(det, 1.0, atol=1e-9) and orth < 1e-9,
                    f'det in [{det.min():.12f}, {det.max():.12f}], max |RR^T - I| = {orth:.2e}')

        print('the display pose is NOT the force frame, and differs by a measured amount')
        ids = [i for i in m['entity_ids'] if i in frame['entities']]
        idx = {i: n for n, i in enumerate(m['entity_ids'])}
        shipped = np.array([frame['entities'][i]['centroid_m'] for i in ids])
        shown = np.array([direct[idx[i]] for i in ids])
        d = np.linalg.norm(shown - shipped, axis=1)
        # zero to machine precision, not exactly zero: the global fit is an SVD, so the
        # residual is float64 round-off. Measured at 1.28e-16 m, which is 0.13 femtometres.
        force_moved = np.abs(np.array([frame['entities'][i]['translation_m'] for i in ids])).max()
        ok &= check('force frame has not moved the anatomy at all', force_moved < 1e-12,
                    f'max translation {force_moved:.3e} m -- zero to machine precision')
        ok &= check('display pose carries the body where it actually is',
                    np.median(d) > 0.05,
                    f'median {np.median(d)*1e3:.2f} mm, max {d.max()*1e3:.2f} mm apart')

        print('nothing is silently left at rest')
        ok &= check('unbound entities are listed with a reason',
                    all(isinstance(v, str) and v for v in dp['unbound'].values()),
                    f'{len(dp["unbound"])} unbound: {sorted(dp["unbound"])}')

        print('payload')
        import json
        ok &= check('a frame carries the compact form', len(json.dumps(dp)) < 50_000,
                    f'{len(json.dumps(dp)):,} bytes/frame against 1,159,032 for the per-entity form')
    finally:
        plant.native.close()
    print('\nALL PASS' if ok else '\nFAILURES ABOVE')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
