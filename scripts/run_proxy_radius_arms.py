#!/usr/bin/env python3
"""Q4: is the ankle collapse the proxy RADIUS or the repartitioned INERTIA?

Pre-registered in docs/FOOT_JOINTS.md before any arm here ran. Q3 excluded the
support plane; these two candidates are what survived, and both travel with the
torso repartition. Protocol, bar and scoring are F2's, reused from
run_foot_joint_arms.py rather than reimplemented.

    .venv/bin/python scripts/run_proxy_radius_arms.py
"""
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
sys.path.insert(0, str(ROOT))
import run_foot_joint_arms as F

BASE = 'data/models/engineering_stance_v1/registration.json'
TWELD = 'data/derived/foot-joint-arms/tweld/registration.json'
OUT = ROOT / 'data/models/articulated_spine_v1/proxy_radius_report.json'


def arm(label, registration, radius=None):
    # NOTE: F.proxy_spheres() RECOMPUTES the radius from inertia -- it reproduces the
    # engine's default rule and by construction cannot see a pinned radius. The first
    # run of this file reported that reconstruction and its ball column was wrong for
    # every pinned arm. Read the engine's own records instead: the contact element's
    # radius_m, and the override and basis out of the run's own execution.json.
    summary, trace = F.supine_tonic(registration, proxy_radius_m=radius)
    w = summary['worst_excursion_rad']
    radii = summary.get('engine_contact_radius_m') or {}
    effective = next((v for k, v in radii.items() if 'torso' in k), None)
    row = {'arm': label, 'registration': registration,
           'requested_torso_radius_m': None if radius is None else radius['torso'],
           'ankle_r': w['ankle_angle_r'], 'ankle_l': w['ankle_angle_l'],
           'plane_m': summary['support_plane_source_x_m'],
           'torso_radius_m': effective,
           'reconstructed_radius_m': (summary.get('proxy_spheres') or {}).get('torso', {}).get('radius_m'),
           'override_in_record': summary.get('proxy_radius_override_m'),
           'basis': summary.get('proxy_radius_basis')}
    print(f"  {label:12s} ankle {row['ankle_r']:7.4f}/{row['ankle_l']:7.4f}  "
          f"ball {row['torso_radius_m']}  plane {row['plane_m']:.6f}", flush=True)
    return row


def main() -> int:
    tweld = TWELD if (ROOT / TWELD).is_file() else None
    if tweld is None:
        cands = sorted((ROOT / 'data/models/articulated_spine_v1').glob('registration*tweld*.json'))
        if not cands:
            print('no tweld registration on disk; arms cannot be built', file=sys.stderr); return 2
        tweld = str(cands[0].relative_to(ROOT))
    print(f'tweld registration: {tweld}\narms:')
    rows = [arm('base', BASE), arm('tweld', tweld)]
    r0 = rows[0]['torso_radius_m']; rt = rows[1]['torso_radius_m']
    rows.append(arm('base@R0', BASE, {'torso': r0}))          # instrument: own value, must reproduce
    rows.append(arm('tweld@RT', tweld, {'torso': rt}))        # instrument: own value, must reproduce
    rows.append(arm('tweld@R0', tweld, {'torso': r0}))        # the separating arm
    rows.append(arm('base@RT', BASE, {'torso': rt}))          # the mirror
    by = {r['arm']: r for r in rows}

    checks = {
        'pin at own radius reproduces base bitwise':
            by['base@R0']['ankle_r'] == by['base']['ankle_r'] and by['base@R0']['ankle_l'] == by['base']['ankle_l'],
        'pin at own radius reproduces tweld bitwise':
            by['tweld@RT']['ankle_r'] == by['tweld']['ankle_r'] and by['tweld@RT']['ankle_l'] == by['tweld']['ankle_l'],
        'base reprints Q1/Q3': abs(by['base']['ankle_r'] - 0.1229) < 5e-4,
        'tweld reprints Q1/Q3': abs(by['tweld']['ankle_r'] - 1.4014) < 5e-4,
        'override read back from each run record':
            all(r['override_in_record'] == ({'torso': r['requested_torso_radius_m']} if r['requested_torso_radius_m'] else None) for r in rows),
        'pinned runs declare a pinned basis':
            all(('PINNED' in (r['basis'] or '')) == (r['requested_torso_radius_m'] is not None) for r in rows),
    }
    print('\ninstrument checks')
    for k, v in checks.items():
        print(f'  {"PASS" if v else "FAIL"}  {k}')

    bar = F.PRIOR['bar']
    recovered = lambda r: max(r['ankle_r'], r['ankle_l']) <= bar
    verdict = ('the RADIUS is the mechanism: a contact artefact' if recovered(by['tweld@R0']) and not recovered(by['base@RT'])
               else 'the INERTIA is the mechanism: not a contact artefact' if not recovered(by['tweld@R0']) and recovered(by['base@RT'])
               else 'NOT SEPARABLE by this arm: each necessary, neither sufficient' if recovered(by['tweld@R0']) and recovered(by['base@RT'])
               else 'BOTH independently sufficient')
    print(f'\nbar {bar}')
    print(f"  tweld@R0 (repartitioned, base ball)  {by['tweld@R0']['ankle_r']:.4f}/{by['tweld@R0']['ankle_l']:.4f}  "
          f"{'recovered' if recovered(by['tweld@R0']) else 'COLLAPSES'}")
    print(f"  base@RT  (base, repartitioned ball)  {by['base@RT']['ankle_r']:.4f}/{by['base@RT']['ankle_l']:.4f}  "
          f"{'fine' if recovered(by['base@RT']) else 'COLLAPSES'}")
    print(f'\nVERDICT: {verdict}')
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({'schema': 'ihm.proxy-radius-arms.v1', 'bar': bar,
                               'arms': rows, 'instrument_checks': checks,
                               'verdict': verdict}, indent=2) + '\n')
    print(f'written: {OUT.relative_to(ROOT)}')
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
