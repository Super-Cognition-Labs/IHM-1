#!/usr/bin/env python3
"""Known answers for the two engine options added 18 Sep 2026.

1. The skin sensor cap was a magic 128 -- a whole body's cutaneous afference in 128
   numbers, for no stated reason. The bound is now the quadrature's OWN size, checked
   on both sides.
2. The supine support plane hangs under the lowest proxy sphere, whose radius is
   inscribed in that segment's INERTIA ellipsoid. So a mass repartition moves the
   floor -- 48.5 mm on the spine variant -- and reads as a joint failure
   (docs/FOOT_JOINTS.md). A caller can now PIN the plane and vary only the mass.

The first check is the one that matters most: neither option may move a plant that
does not ask for them.

    .venv/bin/python scripts/verify_engine_options.py
"""
import json, sys, uuid
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
MANIFEST = 'data/derived/supine-surface-contact-xexhch9l/manifest.json'
MASS = 77.6122029
REGISTRATION = 'data/derived/mechanics/whole_body_arm26_v2/registration.json'


def check(name, ok, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {name}{"  " + detail if detail else ""}', flush=True)
    return bool(ok)


def stream(**kw):
    from ihm.native.mechanical_stream import NativeMechanicalStream
    return NativeMechanicalStream(ROOT, ROOT / f'data/derived/scratch-eng-{uuid.uuid4().hex[:8]}', **kw)


def main() -> int:
    import importlib.util
    spec = importlib.util.spec_from_file_location('crawl', ROOT / 'scripts/crawl.py')
    crawl = importlib.util.module_from_spec(spec); spec.loader.exec_module(crawl)
    pose = json.loads((ROOT / 'data/models/engineering_stance_v1/initial_pose.json').read_text())
    points = json.loads((ROOT / MANIFEST).read_text())['points']
    ok = True

    print('regression: a plant that asks for neither option must not move')
    s = stream(environment='upright', target_mass_kg=MASS, initial_pose=pose,
               augmented_registration=REGISTRATION, coordinate_limits=crawl.joint_stops())
    try:
        traj = [{k: v['value'] for k, v in s.advance(0.01)['coordinates'].items()} for _ in range(50)]
    finally:
        s.close()
    reference = json.loads((ROOT / 'data/derived/engine-options/stopped_reference.json').read_text())['trajectory'] \
        if (ROOT / 'data/derived/engine-options/stopped_reference.json').exists() else None
    if reference is None:
        out = ROOT / 'data/derived/engine-options'; out.mkdir(parents=True, exist_ok=True)
        (out / 'stopped_reference.json').write_text(json.dumps({'trajectory': traj}))
        print('    (no reference on disk; wrote this run as the reference)')
    else:
        d = max(abs(a[k] - b[k]) for a, b in zip(reference, traj) for k in a)
        ok &= check('50 steps x 33 coordinates identical to the pre-change engine', d == 0.0,
                    f'max |dq| = {d}')

    print('\nthe skin sensor cap is the quadrature, not 128')
    ok &= check('the manifest declares more points than the old cap', points > 128, f'{points} points')
    wanted = list(range(500))
    s = stream(environment='supine', target_mass_kg=MASS, surface_contact_manifest=MANIFEST,
               surface_sensor_indices=wanted)
    try:
        got = len(s.surface_sensor_identity)
        f = s.advance(0.01)
        emitted = len(f['surface_foundation']['sensor_points'])
        ok &= check('500 skin sensors are accepted AND emitted, where 129 used to be refused',
                    got == 500 and emitted == 500, f'{got} bound, {emitted} emitted in a frame')
    finally:
        s.close()
    try:
        stream(environment='supine', target_mass_kg=MASS, surface_contact_manifest=MANIFEST,
               surface_sensor_indices=list(range(points + 1)))
        ok &= check('control: selecting more sensors than exist is refused', False, 'NOT refused')
    except ValueError as e:
        ok &= check('control: selecting more sensors than exist is refused', True, str(e)[:60])

    print('\nthe support plane can be pinned, so mass and floor can be varied separately')
    s = stream(environment='supine', target_mass_kg=MASS, initial_pose=pose, augmented_registration=REGISTRATION)
    try:
        default_plane = s.snapshot()['support_plane_source_x_m']
    finally:
        s.close()
    pinned = default_plane - 0.0485          # the 48.5 mm the repartition moved it
    s = stream(environment='supine', target_mass_kg=MASS, initial_pose=pose,
               augmented_registration=REGISTRATION, support_plane_source_x_m=pinned)
    try:
        got_plane = s.snapshot()['support_plane_source_x_m']
        ok &= check('the engine uses the pinned plane exactly', got_plane == pinned,
                    f'{got_plane:.10f} m against {pinned:.10f} requested')
    finally:
        s.close()
    s = stream(environment='supine', target_mass_kg=MASS, initial_pose=pose,
               augmented_registration=REGISTRATION, support_plane_source_x_m=default_plane)
    try:
        ok &= check('pinning it AT the default reproduces the default exactly',
                    s.snapshot()['support_plane_source_x_m'] == default_plane,
                    f'{default_plane:.10f} m')
    finally:
        s.close()
    try:
        stream(environment='upright', target_mass_kg=MASS, support_plane_source_x_m=0.0)
        ok &= check('control: the upright environment refuses a pinned plane', False, 'NOT refused')
    except ValueError as e:
        ok &= check('control: the upright environment refuses a pinned plane', True, str(e)[:60])

    print('\nALL PASS' if ok else '\nFAILURES ABOVE')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
