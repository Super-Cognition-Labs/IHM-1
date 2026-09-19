"""Run the pre-registered arms of `docs/FOOT_GEOMETRY.md` and score them by the
rules written there BEFORE any of them was integrated.

Three questions, in this order:

A  **the sign pattern.**  Engine moment arms at the rest pose for the eight toe
   muscles about `mtp_angle_{l,r}`, on the 2x2 of {shipped, corrected} geometry
   against {shipped fitted paths, toe muscles on GeometryPaths}.  The gate is that
   the toe EXTENSORS and the toe FLEXORS come out with OPPOSITE signs on the
   corrected geometry.  An impossible pattern after the fix means the fix is wrong.

B  **what the toe hinge costs now.**  `docs/FOOT_JOINTS.md`'s Q2 protocols,
   unchanged and imported from its own runner, on the same four arms: supine tonic
   unstopped and stopped, and the 3 s crawl seed unstopped and stopped, against the
   recorded base numbers 0.1229 / 0.1178 / 1.4498 / 0.2196 rad.

C  **mtp welded in the articulated spine variant**, F2's protocol, with F2's bar
   measured on the base model in the same run.  **F2 itself is not rescored.** It
   FAILED on `articulated_spine_v1/registration.json` and stays FAILED; the welded
   model is a different plant and its numbers are a new measurement of it.

Every protocol here is IMPORTED from `scripts.run_foot_joint_arms` rather than
restated, so "same protocol" is a fact about the call graph and not a claim.  That
module is under concurrent edit by another worker; its sha256 at run time goes into
the report, so a later reader can tell which revision produced these numbers.

Resources: ONE native session at a time (each prlimit-capped at 4 GB by
NativeMechanicalStream); `open_stream` waits until MemAvailable is at least 12 GB.
Every arm's summary is written the moment it finishes.

    cd <repo> && OPENBLAS_NUM_THREADS=1 nohup nice -n 10 \\
        .venv/bin/python -u -m scripts.run_corrected_foot_arms --part all \\
        > logs/corrected_foot_arms.log 2>&1 &
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.run_foot_joint_arms import (  # noqa: E402
    close, crawl_seed, crawl_stops, default, model_of, open_stream, supine_tonic,
)

RESULTS = ROOT / 'data/derived/corrected-foot-arms'
REPORT = ROOT / 'data/models/corrected_foot_v1/corrected_foot_report.json'

BASE_REG = 'data/models/engineering_stance_v1/registration.json'
CORRECTED = 'data/models/corrected_foot_v1'
SPINE = 'data/models/articulated_spine_v1'

#: The 2x2.  G = geometry corrected, P = the toe muscles on their own GeometryPaths.
ARMS = (
    ('base', BASE_REG, {'G': 0, 'P': 0}),
    ('base_toe_paths', CORRECTED + '/registration_base_toe_paths.json', {'G': 0, 'P': 1}),
    ('corrected_fitted', CORRECTED + '/registration.json', {'G': 1, 'P': 0}),
    ('corrected_toe_paths', CORRECTED + '/registration_toe_paths.json', {'G': 1, 'P': 1}),
)

SPINE_ARMS = (
    ('base', BASE_REG),
    ('spine', SPINE + '/registration.json'),
    ('spine_mtp_welded', SPINE + '/registration_mtp_welded.json'),
)

SIDES = ('r', 'l')
TOE_MUSCLES = tuple('%s_%s' % (m, s) for s in SIDES for m in ('edl', 'ehl', 'fdl', 'fhl'))
EXTENSORS = ('edl', 'ehl')
FLEXORS = ('fdl', 'fhl')
PROBE_MUSCLES = list(TOE_MUSCLES) + ['soleus_r', 'tibant_r',
                                     'gait2392_ercspn_r', 'gait2392_intobl_r',
                                     'gait2392_extobl_r']
PROBE_COORDS = ['ankle_angle_r', 'ankle_angle_l', 'mtp_angle_r', 'mtp_angle_l',
                'lumbar_extension']

# ------------------------------------------------------------------ the bars
#: S1.  The gate.  Toe extensors and toe flexors must have OPPOSITE signs about
#: mtp_angle, on both sides, with no arm inside the floor -- a near-zero arm must
#: not be allowed to satisfy a sign test on rounding.
SIGN_FLOOR_M = 1.0e-3

#: I6.  The offline prediction, computed 18 Sep 2026 from the XML alone --
#: `arm = a . ((p - T) x u)` with `a = R(orientation) z`, `T` the corrected calcn
#: offset, `p` the last calcaneal point and `u` the unit vector to the first
#: toes-frame point mapped into calcn at q = 0 -- BEFORE any engine ran on the
#: corrected model.  The sign is the engine's own convention, calibrated on the
#: shipped geometry where the same routine reproduces the recorded engine values
#: (-4.7804 / +24.1416 / +8.3729 / +10.3578 mm) to 0.01 mm.  See docs/FOOT_GEOMETRY.md.
PREDICTED_MTP_ARM_M = {'edl': -0.006052, 'ehl': -0.007551,
                       'fdl': +0.006834, 'fhl': +0.007129}
PREDICTION_TOLERANCE_M = 5.0e-5          # 0.05 mm: 5x the agreement measured on the known answer

#: I3.  docs/models/articulated_spine_v1/foot_joints_report.json, 18 Sep 2026,
#: `q2.moment_arms_m_rest_pose.base_toe_geometry_paths` -- full precision, so the
#: bar is a reproduction bar and not a rounding one.
RECORDED_SHIPPED_MTP_ARM_M = {'edl': -0.0047804314685816855, 'ehl': 0.024141634365440916,
                              'fdl': 0.008372874673443862, 'fhl': 0.010357840979963293}
REPRODUCTION_TOLERANCE_M = 1.0e-6

#: I4.  docs/FOOT_JOINTS.md's own control.
SOLEUS_ANKLE_M = -0.0497
SOLEUS_TOLERANCE_M = 5.0e-5

#: I7.  docs/UPPER_BODY_ACTUATION.md:110, which does not record the pose it was
#: measured at.  Declared in advance: a failure here does NOT void this run; it is
#: recorded FAILED and diagnosed as a pose difference or as something else.
LUMBAR_ARM_MM = {'gait2392_ercspn_r': 42.69, 'gait2392_intobl_r': -52.82,
                 'gait2392_extobl_r': -62.79}
LUMBAR_TOLERANCE_MM = 0.05

#: I8.  The base plant's own recorded numbers, from the same report, full precision.
RECORDED_BASE_ANKLE_R = {'supine_unstopped': 0.12288047711162375,
                         'supine_stopped': 0.11780461495122352,
                         'crawl_unstopped': 1.4497937772102554,
                         'crawl_stopped': 0.21958669127417052}
#: I9.  The same report's delta_ankle_vs_base_rad for base_toe_geometry_paths.
RECORDED_TOE_PATH_DELTA = {'supine_unstopped': 0.07836913853629168,
                           'supine_stopped': 0.010618938121581212,
                           'crawl_unstopped': -0.05741417866623877,
                           'crawl_stopped': -0.013730683040444447}
RUN_TOLERANCE_RAD = 5.0e-5
DELTA_TOLERANCE_RAD = 5.0e-4

#: docs/FOOT_JOINTS.md's Q2 reporting line, carried over unchanged: 10% of F2's
#: bar, a reporting convention with no measured provenance.
MATERIAL_RAD = 0.022

#: docs/ARTICULATED_SPINE.md, gate F2, printed to 4 dp.  Reported for context only.
F2_BAR = 0.2202
F2_THORACIC_EXTENSION = 0.3940

PROTOCOLS = ('supine_unstopped', 'supine_stopped', 'crawl_unstopped', 'crawl_stopped')


def sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(name, payload):
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / name).write_text(json.dumps(payload, indent=1, default=default) + '\n')


def probe_moment_arms(registration):
    """Engine moment arms at the rest pose, called TWICE at the same input; the two
    must be equal, or the instrument is not a function and nothing below stands."""
    stream, out = open_stream(registration, 'free')
    try:
        a = stream.moment_arms(muscles=PROBE_MUSCLES, coordinates=PROBE_COORDS)['moment_arms_m']
        b = stream.moment_arms(muscles=PROBE_MUSCLES, coordinates=PROBE_COORDS)['moment_arms_m']
    finally:
        close(stream, out)
    assert a == b, 'moment_arms is not a function of its input on ' + registration
    return a


def sign_pattern(arms):
    """The S1 reading on one arm's moment-arm dictionary."""
    out = {}
    for side in SIDES:
        coordinate = 'mtp_angle_' + side
        got = {m: arms['%s_%s' % (m, side)][coordinate] for m in EXTENSORS + FLEXORS}
        smallest = min(abs(v) for v in got.values())
        ext = {np.sign(got[m]) for m in EXTENSORS}
        flx = {np.sign(got[m]) for m in FLEXORS}
        out[side] = {
            'arms_m': got,
            'extensors_agree': len(ext) == 1,
            'flexors_agree': len(flx) == 1,
            'extensors_oppose_flexors': len(ext) == 1 and len(flx) == 1 and ext != flx,
            'smallest_abs_arm_m': smallest,
            'above_floor': smallest >= SIGN_FLOOR_M,
            'possible': (len(ext) == 1 and len(flx) == 1 and ext != flx
                         and smallest >= SIGN_FLOOR_M),
        }
    out['verdict'] = 'POSSIBLE' if all(out[s]['possible'] for s in SIDES) else 'IMPOSSIBLE'
    return out


def score_arms(measured):
    checks, notes = {}, {}
    # I2 / I5: a fitted path set can never see mtp, whatever the geometry
    for arm in ('base', 'corrected_fitted'):
        checks['I%d %s reads exactly 0 about mtp (all 8)' % (2 if arm == 'base' else 5, arm)] = \
            all(measured[arm]['%s_%s' % (m, s)]['mtp_angle_' + s] == 0
                for s in SIDES for m in EXTENSORS + FLEXORS)
    # I3: the shipped-geometry geometry-path arm reproduces the recorded values
    checks['I3 base_toe_paths reproduces the recorded shipped-geometry arms'] = all(
        abs(measured['base_toe_paths']['%s_%s' % (m, s)]['mtp_angle_' + s]
            - RECORDED_SHIPPED_MTP_ARM_M[m]) <= REPRODUCTION_TOLERANCE_M
        for s in SIDES for m in EXTENSORS + FLEXORS)
    # I4: the ankle still reads as a muscle-driven ankle, on every arm
    checks['I4 soleus_r about ankle_angle_r = -0.0497 on every arm'] = all(
        abs(measured[a]['soleus_r']['ankle_angle_r'] - SOLEUS_ANKLE_M) <= SOLEUS_TOLERANCE_M
        for a in measured)
    # I6: the engine reproduces the offline prediction on the corrected geometry
    worst = max(abs(measured['corrected_toe_paths']['%s_%s' % (m, s)]['mtp_angle_' + s]
                    - PREDICTED_MTP_ARM_M[m])
                for s in SIDES for m in EXTENSORS + FLEXORS)
    checks['I6 corrected arms reproduce the offline prediction within 0.05 mm'] = \
        worst <= PREDICTION_TOLERANCE_M
    notes['I6_worst_disagreement_mm'] = 1000.0 * worst
    # I7: the lumbar table.  Declared non-voiding.
    lumbar = {m: 1000.0 * measured['base'][m]['lumbar_extension'] for m in LUMBAR_ARM_MM}
    checks['I7 lumbar arms reproduce UPPER_BODY_ACTUATION.md:110 (NON-VOIDING)'] = all(
        abs(lumbar[m] - LUMBAR_ARM_MM[m]) <= LUMBAR_TOLERANCE_MM for m in LUMBAR_ARM_MM)
    notes['I7_measured_mm'] = lumbar
    notes['I7_expected_mm'] = LUMBAR_ARM_MM
    voiding = {k: v for k, v in checks.items() if 'NON-VOIDING' not in k}
    return {'instrument_checks': checks, 'instrument_ok': all(voiding.values()),
            'notes': notes,
            'sign_pattern': {a: sign_pattern(measured[a]) for a in
                             ('base_toe_paths', 'corrected_toe_paths')},
            'mtp_arms_mm': {a: {'%s_%s' % (m, s): 1000.0 * measured[a]['%s_%s' % (m, s)]['mtp_angle_' + s]
                                for s in SIDES for m in EXTENSORS + FLEXORS}
                            for a in measured},
            'ankle_arms_mm': {a: {'%s_%s' % (m, s): 1000.0 * measured[a]['%s_%s' % (m, s)]['ankle_angle_' + s]
                                  for s in SIDES for m in EXTENSORS + FLEXORS}
                              for a in measured}}


def run_arms():
    measured = {}
    for arm, registration, _ in ARMS:
        print('A  moment arms %-22s %s' % (arm, registration), flush=True)
        measured[arm] = probe_moment_arms(registration)
        save('moment_arms.json', measured)
        for side in SIDES:
            print('     %s  %s' % (side, '  '.join(
                '%s %+7.3f mm' % (m, 1000 * measured[arm]['%s_%s' % (m, side)]['mtp_angle_' + side])
                for m in EXTENSORS + FLEXORS)), flush=True)
    scored = score_arms(measured)
    save('moment_arms_scored.json', scored)
    print(json.dumps({'instrument_checks': scored['instrument_checks'],
                      'notes': scored['notes'],
                      'verdict': {a: scored['sign_pattern'][a]['verdict']
                                  for a in scored['sign_pattern']}},
                     indent=1, default=default), flush=True)
    return measured, scored


def run_protocols():
    runs = {}
    for protocol in PROTOCOLS:
        for arm, registration, _ in ARMS:
            stops = crawl_stops(model_of(registration)) if protocol.endswith('_stopped') else None
            print('B  %-18s %-22s' % (protocol, arm), flush=True)
            fn = supine_tonic if protocol.startswith('supine') else crawl_seed
            try:
                summary, trace = fn(registration, stops)
            except Exception as exc:
                save('b_%s_%s_FAILED.json' % (protocol, arm),
                     {'registration': registration,
                      'failure': type(exc).__name__ + ': ' + str(exc)[:500]})
                raise
            runs.setdefault(protocol, {})[arm] = summary
            save('b_%s_%s.json' % (protocol, arm), summary)
            save('b_%s_%s_trace.json' % (protocol, arm), trace)
            w, e = summary['worst_excursion_rad'], summary['extrema_rad']
            print('     mtp r %s  ankle r %.4f l %.4f  worst %.4f (%s)  travel %s  %.3f s/adv' % (
                ('[%+.4f, %+.4f]' % tuple(e['mtp_angle_r'])) if 'mtp_angle_r' in e else 'welded',
                w['ankle_angle_r'], w['ankle_angle_l'], max(w.values()), max(w, key=w.get),
                summary.get('pelvis_forward_travel_m'), summary['wall_s_per_advance']), flush=True)
    scored = score_protocols(runs)
    save('protocols_scored.json', scored)
    print(json.dumps(scored['instrument_checks'], indent=1, default=default), flush=True)
    return runs, scored


def score_protocols(runs):
    checks = {}
    for protocol, expected in RECORDED_BASE_ANKLE_R.items():
        if protocol in runs:
            checks['I8 base %s reproduces %.6f' % (protocol, expected)] = \
                abs(runs[protocol]['base']['worst_excursion_rad']['ankle_angle_r']
                    - expected) <= RUN_TOLERANCE_RAD
    table = {}
    for protocol, arms in runs.items():
        base = arms['base']['worst_excursion_rad']
        base_ankle = max(base['ankle_angle_r'], base['ankle_angle_l'])
        rows = {}
        for arm, r in arms.items():
            w = r['worst_excursion_rad']
            other = {k: v for k, v in w.items() if not k.startswith('mtp_angle')}
            ankle = max(w['ankle_angle_r'], w['ankle_angle_l'])
            item = {'ankle_excursion_rad': {s: w['ankle_angle_' + s] for s in SIDES},
                    'delta_ankle_vs_base_rad': ankle - base_ankle,
                    'ankle_changes_materially': abs(ankle - base_ankle) > MATERIAL_RAD,
                    'mtp_extrema_rad': {s: r['extrema_rad'].get('mtp_angle_' + s) for s in SIDES},
                    'mtp_excursion_past_declared_rad': {s: w.get('mtp_angle_' + s) for s in SIDES},
                    'worst_excluding_mtp_rad': max(other.values()),
                    'worst_excluding_mtp_coordinate': max(other, key=other.get),
                    'support_plane_source_x_m': r.get('support_plane_source_x_m'),
                    'wall_s_per_advance': r['wall_s_per_advance']}
            if protocol.startswith('crawl'):
                item.update({k: r.get(k) for k in ('pelvis_forward_travel_m',
                                                   'first_crawl_divergence', 'stop_reason',
                                                   'failure', 'simulated_s')})
            rows[arm] = item
        table[protocol] = rows
        if protocol in RECORDED_TOE_PATH_DELTA:
            checks['I9 base_toe_paths %s reproduces delta %+.4f' % (
                protocol, RECORDED_TOE_PATH_DELTA[protocol])] = abs(
                rows['base_toe_paths']['delta_ankle_vs_base_rad']
                - RECORDED_TOE_PATH_DELTA[protocol]) <= DELTA_TOLERANCE_RAD
    # the 2x2: G is the geometry edit (which also moves the toes body and its two
    # contact spheres 27.8 mm distally), P is the toe muscles' own GeometryPaths
    effects = {}
    for protocol, rows in table.items():
        def A(arm):
            return max(rows[arm]['ankle_excursion_rad'].values())
        effects[protocol] = {
            'G_at_P0 (corrected_fitted - base)': A('corrected_fitted') - A('base'),
            'G_at_P1 (corrected_toe_paths - base_toe_paths)': A('corrected_toe_paths') - A('base_toe_paths'),
            'P_at_G0 (base_toe_paths - base)': A('base_toe_paths') - A('base'),
            'P_at_G1 (corrected_toe_paths - corrected_fitted)': A('corrected_toe_paths') - A('corrected_fitted'),
            'total (corrected_toe_paths - base)': A('corrected_toe_paths') - A('base')}
    return {'instrument_checks': checks, 'instrument_ok': all(checks.values()),
            'table': table, 'effects_on_A_rad': effects,
            'material_reporting_line_rad': MATERIAL_RAD}


def run_spine():
    """C: mtp welded in the articulated spine variant, under F2's protocol.
    F2 itself is NOT rescored -- it failed on registration.json and stays failed."""
    runs = {}
    for arm, registration in SPINE_ARMS:
        print('C  supine tonic %-20s %s' % (arm, registration), flush=True)
        summary, trace = supine_tonic(registration)
        runs[arm] = summary
        save('c_%s.json' % arm, summary)
        save('c_%s_trace.json' % arm, trace)
        w = summary['worst_excursion_rad']
        print('     worst %.4f (%s)  ankle r %.4f l %.4f  plane %s' % (
            max(w.values()), max(w, key=w.get), w['ankle_angle_r'], w['ankle_angle_l'],
            summary.get('support_plane_source_x_m')), flush=True)
    bar = max(runs['base']['worst_excursion_rad'].values())
    scored = {
        'bar_rad': bar,
        'bar_coordinate': max(runs['base']['worst_excursion_rad'],
                              key=runs['base']['worst_excursion_rad'].get),
        'instrument_checks': {
            'bar reproduces F2 0.2202': abs(bar - F2_BAR) <= RUN_TOLERANCE_RAD,
            'spine reproduces F2 thoracic_extension 0.3940':
                abs(runs['spine']['worst_excursion_rad']['thoracic_extension']
                    - F2_THORACIC_EXTENSION) <= RUN_TOLERANCE_RAD},
        'worst_excursions_rad': {a: runs[a]['worst_excursion_rad'] for a in runs},
        'F2': ('FAILED on data/models/articulated_spine_v1/registration.json, '
               'thoracic_extension 0.3940 against a 0.2202 bar. NOT RESCORED. The welded '
               'model is a different plant; the numbers beside it are a new measurement '
               'of that plant, never a new score for F2.'),
    }
    spine, welded = (runs['spine']['worst_excursion_rad'],
                     runs['spine_mtp_welded']['worst_excursion_rad'])
    shared = sorted(set(spine) & set(welded))
    scored['weld_deltas_rad'] = {k: welded[k] - spine[k] for k in shared}
    scored['weld_worst_delta'] = max(((abs(v), k) for k, v in scored['weld_deltas_rad'].items()),
                                     default=(0.0, None))
    scored['welded_over_bar'] = sorted(k for k in welded if welded[k] > bar)
    scored['spine_over_bar'] = sorted(k for k in spine if spine[k] > bar)
    save('c_scored.json', scored)
    print(json.dumps({k: scored[k] for k in ('bar_rad', 'instrument_checks',
                                             'weld_worst_delta', 'welded_over_bar',
                                             'spine_over_bar')},
                     indent=1, default=default), flush=True)
    return runs, scored


def write_report(parts):
    body = {'schema': 'ihm.corrected-foot-report.v1',
            'preregistration': 'docs/FOOT_GEOMETRY.md',
            'protocol_source': {'scripts/run_foot_joint_arms.py':
                                sha(ROOT / 'scripts/run_foot_joint_arms.py')},
            'builders': {p: sha(ROOT / p) for p in
                         ('scripts/build_corrected_foot.py',
                          'scripts/build_mtp_welded_spine.py',
                          'scripts/run_corrected_foot_arms.py')},
            'arms': {name: {'registration': reg, 'factors': f,
                            'model_sha256': json.loads((ROOT / reg).read_text())['model_sha256']}
                     for name, reg, f in ARMS},
            'spine_arms': {name: {'registration': reg,
                                  'model_sha256': json.loads((ROOT / reg).read_text())['model_sha256']}
                           for name, reg in SPINE_ARMS},
            **parts}
    REPORT.write_text(json.dumps(body, indent=1, sort_keys=True, default=default) + '\n')
    print('wrote', REPORT.relative_to(ROOT), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=('arms', 'protocols', 'spine', 'all'), default='all')
    a = ap.parse_args()
    parts = {}
    if REPORT.exists():                     # keep what an earlier part already wrote
        parts = {k: v for k, v in json.loads(REPORT.read_text()).items()
                 if k in ('A_moment_arms', 'B_protocols', 'C_spine_weld')}
    if a.part in ('arms', 'all'):
        measured, scored = run_arms()
        parts['A_moment_arms'] = {'measured_m': measured, 'scored': scored}
        write_report(parts)
    if a.part in ('protocols', 'all'):
        runs, scored = run_protocols()
        parts['B_protocols'] = {'runs': runs, 'scored': scored}
        write_report(parts)
    if a.part in ('spine', 'all'):
        runs, scored = run_spine()
        parts['C_spine_weld'] = {'runs': runs, 'scored': scored}
        write_report(parts)


if __name__ == '__main__':
    main()
