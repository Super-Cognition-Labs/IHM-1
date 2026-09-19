"""Q3: separate the supine support PLANE from the MASS repartition that moved it.

Pre-registration, arms, bar and reading rules: docs/FOOT_JOINTS.md, section
"PRE-REGISTRATION -- Q3", committed before any arm here ran.

docs/FOOT_JOINTS.md's Q1 showed the ankle collapse on articulated_spine_v1 is
carried by the trunk repartition (T +1.254 rad; W -0.001, S -0.019) and persists
with every new joint welded.  The leading candidate was the CONTACT PROXY: the
engine hangs the supine plane under the lowest sphere, each sphere inscribed in
that segment's inertia ellipsoid, so repartitioning the torso grew its ball
0.258 -> 0.309 m and dropped the plane 48.5 mm.  Q1 could not separate the plane
drop from the mass change, because the plane was computed inside the engine.

It now takes an option.  This runs the two arms that decide it between them --

    tweld@P0   repartitioned trunk, plane PINNED at the base plant's  -0.40350 m
    base@PT    base trunk,          plane PINNED at the variant's     -0.45201 m

-- plus the pin-at-its-own-default controls, the six derived-plane arms that
carry Q1's bar and its instrument checks, and a five-point ladder that lowers the
base plant's floor in equal steps.

Everything is Q1's protocol and Q1's scoring, through Q1's own `supine_tonic`,
so the numbers are comparable line for line.  Nothing here rescores F2 or G-S.

Resources: ONE native session at a time (each prlimit-capped at 4 GB by
NativeMechanicalStream); `open_stream` waits for 12 GB of MemAvailable first.
Every arm's summary is written the moment it finishes.

    cd <repo> && OPENBLAS_NUM_THREADS=1 nohup nice -n 10 \\
        .venv/bin/python -u -m scripts.run_plane_arms > logs/plane_arms.log 2>&1 &
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.build_foot_joint_arms import OUT  # noqa: E402
from scripts.run_foot_joint_arms import (  # noqa: E402
    BASE_REG, FOOT_REG, PRINT_TOLERANCE, PRIOR, RESULTS, SIDES, TRACE_BODIES,
    VARIANT_REG, ankle_timing, default, save, supine_tonic,
)

REPORT = ROOT / 'data/models/articulated_spine_v1/plane_arms_report.json'

#: 4-dp prints from a910288 (docs/FOOT_JOINTS.md Q1), on the PREVIOUS engine
#: build.  Reproduction tolerance is that print's rounding, 5e-5.
PRIOR_DERIVED = {'tweld': {'ankle_angle_r': 1.4014, 'ankle_angle_l': 1.2876},
                 't1w0s0': {'ankle_angle_r': 1.3851, 'ankle_angle_l': 1.2663}}

#: The two planes to 5 dp, as Q1 printed them.  The run uses the engine's own
#: full-precision values; these only have to reproduce, or the run is void.
PRINTED_PLANES = {'P0': -0.40350, 'PT': -0.45201}
PLANE_PRINT_TOLERANCE = 5e-6

#: Inscribed-ellipsoid radii Q1 measured for the torso, as a known answer on the
#: reconstruction of the engine's own rule.
PRINTED_TORSO_RADIUS = {'base': 0.2577, 'tweld': 0.3090}
RADIUS_PRINT_TOLERANCE = 5e-5

#: The ladder, in metres below P0.  0.0 is `base@P0` and the last is `base@PT`,
#: so the ladder is five points and two of them are already arms.
LADDER_MM = (0.0, 12.125, 24.25, 36.375, 48.5)

ARM_REG = OUT + '/%s/registration.json'
DERIVED_ARMS = (('base', BASE_REG), ('base_repeat', BASE_REG),
                ('t0w0s0', ARM_REG % 't0w0s0'), ('t1w1s1', VARIANT_REG),
                ('foot_paths', FOOT_REG), ('tweld', ARM_REG % 'tweld'),
                ('t1w0s0', ARM_REG % 't1w0s0'))


def git_sha():
    try:
        return subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, check=True,
                              capture_output=True, text=True).stdout.strip()
    except Exception:
        return 'git-unknown'


def run(arm, registration, plane=None):
    print('Q3 %-22s %-62s plane %s' % (
        arm, registration, 'derived' if plane is None else '%.10f' % plane), flush=True)
    try:
        summary, trace = supine_tonic(registration, support_plane_source_x_m=plane)
    except Exception as exc:
        save('q3_%s_FAILED.json' % arm, {'registration': registration, 'requested_plane': plane,
                                         'failure': type(exc).__name__ + ': ' + str(exc)[:500]})
        raise
    save('q3_%s.json' % arm, summary)
    save('q3_%s_trace.json' % arm, trace)
    w = summary['worst_excursion_rad']
    print('   ankle r %.4f l %.4f   plane %.10f   worst %.4f (%s)   %.3f s/adv' % (
        w['ankle_angle_r'], w['ankle_angle_l'], summary['support_plane_source_x_m'],
        max(w.values()), max(w, key=w.get), summary['wall_s_per_advance']), flush=True)
    return summary, trace


def lowest_sphere(summary):
    """The body the engine's own rule would hang the plane under, and its gap."""
    spheres = summary['proxy_spheres']
    name = min(spheres, key=lambda b: spheres[b]['lowest_point_x_m'])
    return name, spheres[name]


def contact_at(trace, index):
    return {b: f for b, f in trace[index]['contact_n'].items()}


def penetration_report(summary, trace):
    """The declared hazard, measured: how far inside (or above) the plane every
    proxy starts, and what the contact forces do about it."""
    name, low = lowest_sphere(summary)
    peaks = {}
    for r in trace:
        for b, f in r['contact_n'].items():
            peaks[b] = max(peaks.get(b, 0.0), f)
    return {'lowest_sphere_body': name,
            'lowest_sphere_gap_to_plane_m': low['gap_to_plane_m'],
            'initial_penetration_m': max(0.0, -low['gap_to_plane_m']),
            'torso_gap_to_plane_m': summary['proxy_spheres']['torso']['gap_to_plane_m'],
            'torso_radius_m': summary['proxy_spheres']['torso']['radius_m'],
            'foot_gaps_to_plane_m': {b: summary['proxy_spheres'][b]['gap_to_plane_m']
                                     for b in TRACE_BODIES if b in summary['proxy_spheres']},
            'contact_n_at_t0': contact_at(trace, 0),
            'peak_contact_n': dict(peaks),
            'peak_contact_n_all_max': max(peaks.values()) if peaks else 0.0,
            'peak_contact_n_all_max_body': max(peaks, key=peaks.get) if peaks else None}


def plane_rule_check(summary):
    """Check 10: the reconstruction IS the engine's rule.  For a derived plane the
    smallest reconstructed gap must be 0."""
    name, low = lowest_sphere(summary)
    return {'lowest_body': name, 'residual_m': low['gap_to_plane_m'],
            'ok': abs(low['gap_to_plane_m']) <= 1e-12}


def score(runs, traces):
    base, rep = runs['base'], runs['base_repeat']
    bar = max(base['worst_excursion_rad'].values())
    culprit = max(base['worst_excursion_rad'], key=base['worst_excursion_rad'].get)

    checks = {
        '1 base == base_repeat, exact':
            base['worst_excursion_rad'] == rep['worst_excursion_rad'],
        '2 round trip t0w0s0 == base, exact':
            runs['t0w0s0']['worst_excursion_rad'] == base['worst_excursion_rad'],
        '3 bar reproduces F2 0.2202': abs(bar - PRIOR['bar']) <= PRINT_TOLERANCE,
    }
    for arm in ('t1w1s1', 'foot_paths'):
        for k, v in PRIOR[arm].items():
            checks['%s %s reproduces %.4f' % (arm, k, v)] = \
                abs(runs[arm]['worst_excursion_rad'][k] - v) <= PRINT_TOLERANCE
    q1_checks = dict(checks)               # the seven; 8-10 are scored apart

    new_checks = {
        '8a base@P0 bit-equal to base':
            runs['base@P0']['worst_excursion_rad'] == base['worst_excursion_rad'],
        '8b tweld@PT bit-equal to tweld':
            runs['tweld@PT']['worst_excursion_rad'] == runs['tweld']['worst_excursion_rad'],
    }
    pinned = {a: r for a, r in runs.items() if r['support_plane_requested_x_m'] is not None}
    new_checks['9a every pinned arm reports the pinned plane exactly'] = all(
        r['support_plane_source_x_m'] == r['support_plane_requested_x_m'] for r in pinned.values())
    new_checks['9b every pinned arm records the override in its own execution.json'] = all(
        r['support_plane_override_x_m'] == r['support_plane_requested_x_m']
        and 'PINNED' in (r['support_plane_basis'] or '') for r in pinned.values())
    new_checks['9c no derived arm records an override'] = all(
        r['support_plane_override_x_m'] is None and 'PINNED' not in (r['support_plane_basis'] or '')
        for a, r in runs.items() if r['support_plane_requested_x_m'] is None)
    new_checks['9d the pin bites: tweld@P0 plane differs from tweld by 48.5 mm'] = \
        abs((runs['tweld@P0']['support_plane_source_x_m']
             - runs['tweld']['support_plane_source_x_m']) - 0.0485) < 5e-5
    rule = {a: plane_rule_check(r) for a, r in runs.items()
            if r['support_plane_requested_x_m'] is None}
    new_checks['10a derived plane == lowest reconstructed sphere, every arm'] = \
        all(v['ok'] for v in rule.values())
    new_checks['10b the lowest body is torso, every derived arm'] = \
        all(v['lowest_body'] == 'torso' for v in rule.values())
    new_checks['10c torso radius reproduces 0.2577 (base) and 0.3090 (tweld)'] = all(
        abs(runs[a]['proxy_spheres']['torso']['radius_m'] - v) <= RADIUS_PRINT_TOLERANCE
        for a, v in PRINTED_TORSO_RADIUS.items())

    planes = {'P0': runs['base']['support_plane_source_x_m'],
              'PT': runs['tweld']['support_plane_source_x_m']}
    plane_reproduction = {k: {'measured_m': planes[k], 'printed_m': PRINTED_PLANES[k],
                              'ok': abs(planes[k] - PRINTED_PLANES[k]) <= PLANE_PRINT_TOLERANCE}
                          for k in planes}
    derived_reproduction = {
        '%s %s reproduces %.4f' % (arm, k, v):
            abs(runs[arm]['worst_excursion_rad'][k] - v) <= PRINT_TOLERANCE
        for arm, d in PRIOR_DERIVED.items() for k, v in d.items()}

    table = {}
    for arm, r in runs.items():
        w = r['worst_excursion_rad']
        table[arm] = {
            'registration': r['registration'],
            'ankle_angle_r': w['ankle_angle_r'], 'ankle_angle_l': w['ankle_angle_l'],
            'A': max(w['ankle_angle_r'], w['ankle_angle_l']),
            'recovered': w['ankle_angle_r'] <= bar and w['ankle_angle_l'] <= bar,
            'worst_rad': max(w.values()), 'worst_coordinate': max(w, key=w.get),
            'plane_m': r['support_plane_source_x_m'],
            'plane_pinned': r['support_plane_requested_x_m'] is not None,
            'plane_vs_P0_mm': (r['support_plane_source_x_m'] - planes['P0']) * 1000.0,
            'hazard': penetration_report(r, traces[arm]),
            'wall_s_per_advance': r['wall_s_per_advance']}

    rec = {a: table[a]['recovered'] for a in table}
    quadrant = ('recovers' if rec['tweld@P0'] else 'collapses',
                'recovers' if rec['base@PT'] else 'collapses')
    meaning = {
        ('recovers', 'collapses'):
            'THE PLANE. The 48.5 mm drop is sufficient on an unchanged mass and necessary '
            'under the repartition: the collapse is the contact proxy\'s plane placement.',
        ('collapses', 'recovers'):
            'THE MASS. The plane drop is neither sufficient nor necessary; the leading '
            'candidate in docs/FOOT_JOINTS.md is WITHDRAWN.',
        ('recovers', 'recovers'):
            'NEITHER ALONE. Each term is necessary and neither is sufficient; the '
            'single-factor story is wrong and the candidate is downgraded, not confirmed.',
        ('collapses', 'collapses'):
            'OVER-DETERMINED, and this experiment cannot assign it: two independently '
            'sufficient routes, or a third cause common to both. Not a separation.',
    }[quadrant]

    ladder = []
    for mm in LADDER_MM:
        arm = 'base@P0' if mm == 0.0 else ('base@PT' if mm == LADDER_MM[-1]
                                           else 'base@P0-%gmm' % mm)
        ladder.append({'drop_mm': mm, 'arm': arm, 'A': table[arm]['A'],
                       'ankle_angle_r': table[arm]['ankle_angle_r'],
                       'ankle_angle_l': table[arm]['ankle_angle_l'],
                       'recovered': table[arm]['recovered'],
                       'plane_m': table[arm]['plane_m']})
    a_values = [p['A'] for p in ladder]
    steps = np.diff(a_values)
    ladder_shape = ('monotone increasing' if all(d > 0 for d in steps) else
                    'monotone decreasing' if all(d < 0 for d in steps) else 'NON-MONOTONE')
    first_collapse = next((p['drop_mm'] for p in ladder if not p['recovered']), None)

    voided = not all(new_checks.values())
    return {
        'bar_rad': bar, 'bar_coordinate': culprit,
        'planes_m': planes, 'plane_reproduction': plane_reproduction,
        'instrument_checks_q1_seven': q1_checks,
        'instrument_checks_q1_seven_ok': all(q1_checks.values()),
        'instrument_checks_new': new_checks,
        'instrument_checks_new_ok': not voided,
        'plane_rule_reconstruction': rule,
        'derived_arm_reproduction': derived_reproduction,
        'table': table,
        'verdict': {
            'tweld@P0 (repartitioned mass, floor held UP)': quadrant[0],
            'base@PT (base mass, floor lowered)': quadrant[1],
            'reading': 'VOID -- a new instrument check failed; no separation is claimed'
                       if voided else meaning,
            'secondary_t1w0s0@P0': 'recovers' if rec['t1w0s0@P0'] else 'collapses',
            'A_rad': {a: table[a]['A'] for a in
                      ('base', 'base@P0', 'base@PT', 'tweld', 'tweld@PT', 'tweld@P0',
                       't1w0s0', 't1w0s0@P0')},
        },
        'ladder': {'points': ladder, 'shape': ladder_shape,
                   'first_drop_mm_that_collapses': first_collapse,
                   'A_span_rad': max(a_values) - min(a_values)},
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true',
                    help='print the plan and exit without opening a native session')
    a = ap.parse_args()

    runs, traces = {}, {}
    if a.dry_run:
        for arm, reg in DERIVED_ARMS:
            print('derived  %-22s %s' % (arm, reg))
        for arm in ('base@P0', 'tweld@PT', 'tweld@P0', 'base@PT', 't1w0s0@P0'):
            print('pinned   %s' % arm)
        for mm in LADDER_MM[1:-1]:
            print('ladder   base@P0-%gmm' % mm)
        return 0

    started = time.time()
    for arm, reg in DERIVED_ARMS:
        runs[arm], traces[arm] = run(arm, reg)

    p0 = runs['base']['support_plane_source_x_m']
    pt = runs['tweld']['support_plane_source_x_m']
    print('\nP0 = %.12f m (base)   PT = %.12f m (tweld)   drop %.4f mm\n'
          % (p0, pt, (p0 - pt) * 1000.0), flush=True)

    for arm, reg, plane in (('base@P0', BASE_REG, p0),
                            ('tweld@PT', ARM_REG % 'tweld', pt),
                            ('tweld@P0', ARM_REG % 'tweld', p0),
                            ('base@PT', BASE_REG, pt),
                            ('t1w0s0@P0', ARM_REG % 't1w0s0', p0)):
        runs[arm], traces[arm] = run(arm, reg, plane)
    for mm in LADDER_MM[1:-1]:
        arm = 'base@P0-%gmm' % mm
        runs[arm], traces[arm] = run(arm, BASE_REG, p0 - mm / 1000.0)

    scored = score(runs, traces)
    scored['timing'] = {arm: ankle_timing(traces[arm], scored['bar_rad'], runs[arm]['declared'])
                        for arm in runs}
    save('q3_scored.json', scored)

    report = {
        'schema': 'ihm.plane-arms-report.v1',
        'preregistration': 'docs/FOOT_JOINTS.md, section "PRE-REGISTRATION -- Q3"',
        'question': 'Does the supine support PLANE\'s 48.5 mm drop, or the trunk mass '
                    'repartition that moved it, drive the ankle collapse of '
                    'articulated_spine_v1? The scaffold throughout; not a statement about '
                    'a human foot.',
        'git_sha': git_sha(),
        'engine_build': runs['base']['engine_build'],
        'wall_s': time.time() - started,
        'arms': {arm: {'registration': r['registration'], 'model': r['model'],
                       'support_plane_source_x_m': r['support_plane_source_x_m'],
                       'support_plane_override_x_m': r['support_plane_override_x_m'],
                       'support_plane_basis': r['support_plane_basis'],
                       'worst_excursion_rad': r['worst_excursion_rad'],
                       'proxy_spheres': r['proxy_spheres']}
                 for arm, r in runs.items()},
        'scored': scored,
    }
    REPORT.write_text(json.dumps(report, indent=1, default=default) + '\n')

    print(json.dumps({k: scored[k] for k in
                      ('bar_rad', 'planes_m', 'plane_reproduction',
                       'instrument_checks_q1_seven', 'instrument_checks_new',
                       'derived_arm_reproduction', 'verdict', 'ladder')},
                     indent=1, default=default), flush=True)
    print('\nreport: %s' % REPORT.relative_to(ROOT), flush=True)
    print('results: %s' % RESULTS.relative_to(ROOT), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
