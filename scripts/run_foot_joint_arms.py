"""Run the pre-registered foot-joint arms (docs/FOOT_JOINTS.md) and score them
by the rules written there BEFORE any of them was integrated.

Q1  why the ankle collapses on articulated_spine_v1: a 2x2x2 factorial over
    T (trunk partition + spine joints), W (wrists), S (subtalar), plus `tweld`
    (T's bodies with its joints welded) and the G-S foot-paths plant, all under
    F2/G-S's protocol -- supine, 0.02 tonic on every muscle, 2 s, worst
    excursion past each declared range, bar = the BASE model's own worst
    measured in the same run.
Q2  what the base plant's undriven toe hinges cost: engineering_stance_v1 vs
    mtp welded vs the toe muscles on GeometryPaths, under the supine tonic
    protocol unstopped and stopped, and the crawl seed pattern (3 s, prone,
    upright environment) unstopped and stopped -- the protocols of
    docs/NATIVE_JOINT_LIMITS.md.

Resources: ONE native session at a time (each prlimit-capped at 4 GB by
NativeMechanicalStream); before each session this waits until MemAvailable is
at least 12 GB.  Every arm's summary is written the moment it finishes.

    cd <repo> && OPENBLAS_NUM_THREADS=1 nohup nice -n 10 \\
        .venv/bin/python -u -m scripts.run_foot_joint_arms --part all \\
        > logs/foot_joint_arms.log 2>&1 &
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import time
import uuid
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402
from scripts import crawl  # noqa: E402
from scripts.build_foot_joint_arms import OUT, TOE_MUSCLES  # noqa: E402
from scripts.verify_articulated_spine import TARGET_MASS_KG, coordinates  # noqa: E402

RESULTS = ROOT / OUT / 'results'
BASE_REG = 'data/models/engineering_stance_v1/registration.json'
VARIANT_REG = 'data/models/articulated_spine_v1/registration.json'
FOOT_REG = 'data/models/articulated_spine_v1/registration_foot_paths.json'

MIN_AVAILABLE_GB = 12.0
TONIC = 0.02
SUPINE_STEPS = 200
CRAWL_STEPS = 300           # 3 s: the horizon of the NATIVE_JOINT_LIMITS stop sweep
CRAWL_WALL_BUDGET_S = 1200.0

#: Printed by F2 / G-S on 2026-09-18 to 4 dp.  Reproduction tolerance is the
#: rounding of that print, 5e-5 -- not a physical tolerance.
PRIOR = {'bar': 0.2202,
         't1w1s1': {'ankle_angle_r': 1.3563, 'ankle_angle_l': 1.3522},
         'foot_paths': {'ankle_angle_r': 1.1399, 'ankle_angle_l': 0.8709}}
PRINT_TOLERANCE = 5e-5

#: Reporting threshold for Q2, NOT a measured quantity: 10% of F2's bar.
Q2_MATERIAL_RAD = 0.022

SIDES = ('r', 'l')
TRACE_COORDS = tuple('%s_%s' % (c, s) for s in SIDES for c in
                     ('ankle_angle', 'subtalar_angle', 'mtp_angle', 'knee_angle', 'hip_flexion'))
TRACE_BODIES = tuple('%s_%s' % (b, s) for s in SIDES for b in ('tibia', 'talus', 'calcn', 'toes'))
TRACE_MUSCLES = tuple('%s_%s' % (m, s) for s in SIDES for m in
                      ('soleus', 'gasmed', 'gaslat', 'tibpost', 'perlong', 'perbrev',
                       'fhl', 'fdl', 'tibant', 'edl', 'ehl'))


def default(value):
    """Never let a stray numpy value destroy a record (IBM-1 CLAUDE.md)."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return repr(value)


def save(name, payload):
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / name).write_text(json.dumps(payload, indent=1, default=default) + '\n')


def available_gb():
    for line in Path('/proc/meminfo').read_text().splitlines():
        if line.startswith('MemAvailable:'):
            return int(line.split()[1]) / 1024 / 1024
    return 0.0


def wait_for_memory():
    waited = 0
    while available_gb() < MIN_AVAILABLE_GB:
        if waited % 300 == 0:
            print('  waiting: %.1f GB available < %.0f' % (available_gb(), MIN_AVAILABLE_GB), flush=True)
        time.sleep(30)
        waited += 30


def model_of(registration):
    return ROOT / json.loads((ROOT / registration).read_text())['model_path']


def open_stream(registration, environment, stops=None, pose=None):
    wait_for_memory()
    out = 'data/derived/foot-joint-arms/sessions/' + uuid.uuid4().hex[:12]
    stream = NativeMechanicalStream(ROOT, ROOT / out, environment=environment,
                                    target_mass_kg=TARGET_MASS_KG,
                                    augmented_registration=registration,
                                    coordinate_limits=stops, initial_pose=pose)
    return stream, out


def close(stream, out):
    try:
        stream.close()
    finally:
        shutil.rmtree(ROOT / out, ignore_errors=True)


def crawl_stops(model_path):
    """crawl.py's stop set, restricted to coordinates the arm actually has."""
    ranges = crawl.declared_ranges(model_path)
    return [dict(coordinate=n, lower_rad=lo, upper_rad=hi, **crawl.JOINT_STOP)
            for n, (lo, hi) in sorted(ranges.items())]


def row(state, t):
    q = state['coordinates']
    # summed per body: in the upright environment calcn carries four source
    # contact spheres and toes two, so a dict keyed by body would keep only one
    contacts = {}
    for c in state['contacts']:
        f = math.sqrt(sum(v * v for v in c['force_n']))
        contacts[c['body_frame']] = contacts.get(c['body_frame'], 0.0) + f
    toes = {b: contacts.get(b, 0.0) for b in ('toes_r', 'toes_l')}
    return {'t': t,
            'q': {k: q[k]['value'] for k in TRACE_COORDS if k in q},
            'qd': {k: q[k]['speed'] for k in TRACE_COORDS if k in q},
            'contact_n': {b: contacts.get(b, 0.0) for b in TRACE_BODIES},
            'toe_contact_n': toes,
            'tendon_n': {m: state['muscles'][m]['tendon_force_n']
                         for m in TRACE_MUSCLES if m in state['muscles']}}


# ------------------------------------------------------------------ protocols
def supine_tonic(registration, stops=None):
    """F2's `tonic_excursions`, line for line, plus a trace."""
    model = model_of(registration)
    declared = {k: v for k, v in coordinates(model).items() if v is not None}
    stream, out = open_stream(registration, 'supine', stops=stops)
    try:
        first = stream.snapshot()
        plane = first['support_plane_source_x_m']
        excitation = {m: TONIC for m in first['muscles']}
        worst, extrema, trace = {}, {}, [row(first, 0.0)]
        started = time.time()
        for _ in range(SUPINE_STEPS):
            state = stream.advance(0.01, actuation=excitation)
            for name, bounds in declared.items():
                value = state['coordinates'][name]['value']
                worst[name] = max(worst.get(name, -1e9), bounds[0] - value, value - bounds[1])
                lo, hi = extrema.get(name, (value, value))
                extrema[name] = (min(lo, value), max(hi, value))
            trace.append(row(state, state['time_s']))
        seconds = (time.time() - started) / SUPINE_STEPS
    finally:
        close(stream, out)
    return {'registration': registration, 'model': str(model.relative_to(ROOT)),
            'stopped': stops is not None, 'n_stops': 0 if stops is None else len(stops),
            'support_plane_source_x_m': plane, 'worst_excursion_rad': worst,
            'extrema_rad': extrema, 'declared': declared,
            'wall_s_per_advance': seconds}, trace


def crawl_seed(registration, stops=None):
    """The crawl seed pattern, prone, 3 s -- the stop sweep's protocol.  Runs the
    full horizon unless the native call fails or the wall budget is spent;
    where crawl.diverged()'s rules would have stopped a search, the first such
    event is RECORDED, not acted on."""
    model = model_of(registration)
    ranges = crawl.declared_ranges(model)
    stream, out = open_stream(registration, 'upright', stops=stops, pose=crawl.PRONE_POSE)
    trace, worst, extrema = [], {}, {}
    first_divergence = None
    fibre_peak = (0.0, None)
    failure, reason = None, 'horizon'
    started = time.time()
    try:
        state = stream.snapshot()
        names = sorted(state['muscles'])
        pattern = crawl.CrawlPattern(crawl.SEED, names)
        tx0 = state['coordinates']['pelvis_tx']['value']
        for _ in range(CRAWL_STEPS):
            t = state['time_s']
            q = state['coordinates']
            for name, (lo, hi) in ranges.items():
                v = q[name]['value']
                worst[name] = max(worst.get(name, -1e9), lo - v, v - hi)
                a, b = extrema.get(name, (v, v))
                extrema[name] = (min(a, v), max(b, v))
            f, m = crawl.maximum_normalized_fibre_velocity(state)
            if f > fibre_peak[0]:
                fibre_peak = (f, m)
            if first_divergence is None:
                why = divergence(state, ranges)
                if why is not None:
                    first_divergence = {'time_s': t, 'reason': why}
            trace.append(row(state, t))
            if time.time() - started > CRAWL_WALL_BUDGET_S:
                reason = 'wall_budget'
                break
            u = pattern.excitation(t)
            state = stream.advance(crawl.DT, actuation=dict(zip(names, map(float, u))),
                                   coordinate_actuation=crawl.limb_pd(state, pattern.limb_targets(t)))
        trace.append(row(state, state['time_s']))
        travel = state['coordinates']['pelvis_tx']['value'] - tx0
        simulated = state['time_s']
    except Exception as exc:                       # recorded, never swallowed silently
        failure, reason = type(exc).__name__ + ': ' + str(exc)[:300], 'native_failure'
        travel, simulated = None, trace[-1]['t'] if trace else 0.0
    finally:
        close(stream, out)
    wall = time.time() - started
    return {'registration': registration, 'model': str(model.relative_to(ROOT)),
            'stopped': stops is not None, 'n_stops': 0 if stops is None else len(stops),
            'stop_reason': reason, 'failure': failure, 'simulated_s': simulated,
            'pelvis_forward_travel_m': travel, 'first_crawl_divergence': first_divergence,
            'peak_normalized_fibre_velocity': fibre_peak,
            'worst_excursion_rad': worst, 'extrema_rad': extrema,
            'wall_s_per_advance': wall / max(1, len(trace) - 1)}, trace


def divergence(state, ranges):
    """crawl.diverged()'s rules, with the ARM's declared ranges."""
    q = state['coordinates']
    if q['pelvis_ty']['value'] > 1.20:
        return 'pelvis_above_1.20m'
    for name, c in q.items():
        if abs(c['speed']) > 25.0:
            return 'joint_speed_beyond_25_' + name
    f, m = crawl.maximum_normalized_fibre_velocity(state)
    if f > crawl.FIBRE_VELOCITY_GUARD_OFL_S:
        return 'fibre_velocity_beyond_15_ofl_s_' + m
    past, culprit = 0.0, None
    for name, (lo, hi) in ranges.items():
        p = max(lo - q[name]['value'], q[name]['value'] - hi, 0.0)
        if p > past:
            past, culprit = p, name
    if past > crawl.DECLARED_RANGE_TOLERANCE_RAD:
        return 'beyond_declared_range_%s_by_%.2frad' % (culprit, past)
    return None


def moment_arms(registration):
    """Engine moment arms at the rest pose, called twice; the two must be equal."""
    stream, out = open_stream(registration, 'free')
    try:
        muscles = list(TOE_MUSCLES) + ['soleus_r', 'tibant_r']
        coords = ['ankle_angle_r', 'ankle_angle_l', 'mtp_angle_r', 'mtp_angle_l']
        a = stream.moment_arms(muscles=muscles, coordinates=coords)['moment_arms_m']
        b = stream.moment_arms(muscles=muscles, coordinates=coords)['moment_arms_m']
    finally:
        close(stream, out)
    assert a == b, 'moment_arms is not a function of its input'
    return a


# ------------------------------------------------------------------- analysis
def ankle_timing(trace, bar, declared):
    """Per side: when the ankle crosses the bar, when the subtalar reaches its
    own bound, and how the two co-vary over the whole 2 s."""
    out = {}
    for s in SIDES:
        a, st = 'ankle_angle_' + s, 'subtalar_angle_' + s
        lo, hi = declared[a]
        t = np.array([r['t'] for r in trace])
        qa = np.array([r['q'][a] for r in trace])
        ex = np.maximum(lo - qa, qa - hi)
        cross = np.nonzero(ex > bar)[0]
        side = {'ankle_min_rad': float(qa.min()), 'ankle_max_rad': float(qa.max()),
                'direction': 'plantarflexion (negative)' if -qa.min() > qa.max() else 'dorsiflexion (positive)',
                't_ankle_past_bar_s': float(t[cross[0]]) if len(cross) else None}
        if st in trace[0]['q']:
            slo, shi = declared[st]
            qs = np.array([r['q'][st] for r in trace])
            sx = np.maximum(slo - qs, qs - shi)
            reach = np.nonzero(sx >= 0)[0]
            side.update({
                't_subtalar_at_bound_s': float(t[reach[0]]) if len(reach) else None,
                'subtalar_at_ankle_crossing_rad': float(qs[cross[0]]) if len(cross) else None,
                'pearson_ankle_subtalar': float(np.corrcoef(qa, qs)[0, 1]),
                'subtalar_min_rad': float(qs.min()), 'subtalar_max_rad': float(qs.max())})
            ts, ta = side['t_subtalar_at_bound_s'], side['t_ankle_past_bar_s']
            if ta is None:
                side['reading'] = 'no collapse on this side'
            elif ts is not None and ts < ta:
                side['reading'] = 'subtalar at its bound BEFORE the ankle crossed the bar: consistent with a loaded hinge'
            else:
                side['reading'] = 'ankle crossed the bar BEFORE (or without) the subtalar reaching its bound: the hinge is not the onset driver'
        # contact on the foot at the moment of crossing, and plantarflexor load
        k = cross[0] if len(cross) else len(trace) - 1
        side['foot_contact_at_crossing_n'] = {b: trace[k]['contact_n'][b] for b in TRACE_BODIES if b.endswith('_' + s)}
        side['plantarflexor_tendon_n_at_crossing'] = sum(
            trace[k]['tendon_n'].get('%s_%s' % (m, s), 0.0)
            for m in ('soleus', 'gasmed', 'gaslat', 'tibpost', 'perlong', 'perbrev', 'fhl', 'fdl'))
        side['dorsiflexor_tendon_n_at_crossing'] = sum(
            trace[k]['tendon_n'].get('%s_%s' % (m, s), 0.0) for m in ('tibant', 'edl', 'ehl'))
        out[s] = side
    return out


def score_q1(runs):
    base, bar = runs['base'], max(runs['base']['worst_excursion_rad'].values())
    culprit = max(base['worst_excursion_rad'], key=base['worst_excursion_rad'].get)
    checks = {
        'base_is_a_function (base == base_repeat, exact)':
            runs['base']['worst_excursion_rad'] == runs['base_repeat']['worst_excursion_rad'],
        'round_trip t0w0s0 == base (exact)':
            runs['t0w0s0']['worst_excursion_rad'] == runs['base']['worst_excursion_rad'],
        'bar reproduces F2 0.2202':
            abs(bar - PRIOR['bar']) <= PRINT_TOLERANCE,
    }
    for arm in ('t1w1s1', 'foot_paths'):
        for k, v in PRIOR[arm].items():
            checks['%s %s reproduces %.4f' % (arm, k, v)] = \
                abs(runs[arm]['worst_excursion_rad'][k] - v) <= PRINT_TOLERANCE
    table = {}
    for arm, r in runs.items():
        w = r['worst_excursion_rad']
        a = max(w['ankle_angle_r'], w['ankle_angle_l'])
        table[arm] = {'ankle_angle_r': w['ankle_angle_r'], 'ankle_angle_l': w['ankle_angle_l'],
                      'A': a, 'recovered': w['ankle_angle_r'] <= bar and w['ankle_angle_l'] <= bar,
                      'subtalar_angle_r': w.get('subtalar_angle_r'),
                      'subtalar_angle_l': w.get('subtalar_angle_l'),
                      'support_plane_source_x_m': r['support_plane_source_x_m'],
                      'wall_s_per_advance': r['wall_s_per_advance']}
    A = {arm: table[arm]['A'] for arm in table}
    effects = {}
    for i, factor in enumerate('TWS'):
        pairs = []
        for t in (0, 1):
            for w in (0, 1):
                for s in (0, 1):
                    bits = [t, w, s]
                    if bits[i] == 1:
                        continue
                    lo = 't%dw%ds%d' % tuple(bits)
                    bits[i] = 1
                    hi = 't%dw%ds%d' % tuple(bits)
                    pairs.append({'without': lo, 'with': hi, 'delta_A': A[hi] - A[lo]})
        effects[factor] = {'mean_delta_A_rad': float(np.mean([p['delta_A'] for p in pairs])),
                           'pairs': pairs}
    rec = {arm: table[arm]['recovered'] for arm in table}
    verdict = {
        'separating_control_t1w1s0': 'RECOVERS' if rec['t1w1s0'] else 'DOES NOT RECOVER',
        'S_alone_t0w0s1': 'recovered' if rec['t0w0s1'] else 'COLLAPSES',
        'T_alone_t1w0s0': 'recovered' if rec['t1w0s0'] else 'COLLAPSES',
        'W_alone_t0w1s0': 'recovered' if rec['t0w1s0'] else 'COLLAPSES',
        'tweld': 'recovered' if rec['tweld'] else 'COLLAPSES',
    }
    collapsed = sorted(a for a in FACTORIAL_ARMS if not rec[a])
    needs = {f: all(('%s1' % f.lower()) in a for a in collapsed) for f in 'TWS'}
    verdict['collapsed_arms'] = collapsed
    verdict['factor_present_in_every_collapsed_arm'] = needs
    verdict['sufficient_alone'] = {f: not rec[{'T': 't1w0s0', 'W': 't0w1s0', 'S': 't0w0s1'}[f]]
                                   for f in 'TWS'}
    return {'bar_rad': bar, 'bar_coordinate': culprit, 'instrument_checks': checks,
            'instrument_ok': all(checks.values()), 'table': table, 'effects': effects,
            'verdict': verdict}


FACTORIAL_ARMS = tuple('t%dw%ds%d' % (t, w, s) for t in (0, 1) for w in (0, 1) for s in (0, 1))


def run_q1():
    arm_dir = OUT + '/%s/registration.json'
    plan = [('base', BASE_REG), ('base_repeat', BASE_REG)]
    plan += [(a, VARIANT_REG if a == 't1w1s1' else arm_dir % a) for a in FACTORIAL_ARMS]
    plan += [('tweld', arm_dir % 'tweld'), ('foot_paths', FOOT_REG)]
    runs, traces = {}, {}
    for arm, registration in plan:
        print('Q1 %-12s %s' % (arm, registration), flush=True)
        try:
            summary, trace = supine_tonic(registration)
        except Exception as exc:
            # recorded and fatal for scoring: every Q1 reading needs every arm
            save('q1_%s_FAILED.json' % arm, {'registration': registration,
                                              'failure': type(exc).__name__ + ': ' + str(exc)[:500]})
            raise
        runs[arm], traces[arm] = summary, trace
        save('q1_%s.json' % arm, summary)
        save('q1_%s_trace.json' % arm, trace)
        w = summary['worst_excursion_rad']
        print('   ankle r %.4f l %.4f  subtalar %s  plane %.4f  %.3f s/adv' % (
            w['ankle_angle_r'], w['ankle_angle_l'],
            ('%.4f/%.4f' % (w['subtalar_angle_r'], w['subtalar_angle_l'])) if 'subtalar_angle_r' in w else '-',
            summary['support_plane_source_x_m'], summary['wall_s_per_advance']), flush=True)
    scored = score_q1(runs)
    bar = scored['bar_rad']
    scored['timing'] = {arm: ankle_timing(traces[arm], bar, runs[arm]['declared'])
                        for arm in runs if arm in FACTORIAL_ARMS + ('tweld', 'foot_paths', 'base')}
    save('q1_scored.json', scored)
    print(json.dumps({k: scored[k] for k in ('bar_rad', 'instrument_checks', 'verdict')},
                     indent=1, default=default), flush=True)
    return scored


Q2_ARMS = (('base', BASE_REG),
           ('base_mtp_welded', OUT + '/base_mtp_welded/registration.json'),
           ('base_toe_geometry_paths', OUT + '/base_toe_geometry_paths/registration.json'))


def run_q2():
    results = {'moment_arms': {}, 'runs': {}}
    for arm, registration in (Q2_ARMS[0], Q2_ARMS[2]):
        print('Q2 moment arms %s' % arm, flush=True)
        results['moment_arms'][arm] = moment_arms(registration)
        save('q2_moment_arms.json', results['moment_arms'])
    for protocol in ('supine_unstopped', 'supine_stopped', 'crawl_unstopped', 'crawl_stopped'):
        for arm, registration in Q2_ARMS:
            stops = crawl_stops(model_of(registration)) if protocol.endswith('_stopped') else None
            print('Q2 %-17s %-24s' % (protocol, arm), flush=True)
            fn = supine_tonic if protocol.startswith('supine') else crawl_seed
            try:
                summary, trace = fn(registration, stops)
            except Exception as exc:
                save('q2_%s_%s_FAILED.json' % (protocol, arm), {
                    'registration': registration,
                    'failure': type(exc).__name__ + ': ' + str(exc)[:500]})
                raise
            results['runs'].setdefault(protocol, {})[arm] = summary
            save('q2_%s_%s.json' % (protocol, arm), summary)
            save('q2_%s_%s_trace.json' % (protocol, arm), trace)
            w, e = summary['worst_excursion_rad'], summary['extrema_rad']
            print('   mtp r %s  ankle r %.4f l %.4f  worst %.4f (%s)  travel %s  %.3f s/adv' % (
                ('[%.4f, %.4f]' % tuple(e['mtp_angle_r'])) if 'mtp_angle_r' in e else 'welded',
                w['ankle_angle_r'], w['ankle_angle_l'], max(w.values()), max(w, key=w.get),
                summary.get('pelvis_forward_travel_m'), summary['wall_s_per_advance']), flush=True)
    results['scored'] = score_q2(results)
    save('q2_scored.json', results['scored'])
    print(json.dumps(results['scored'], indent=1, default=default), flush=True)
    return results


def score_q2(results):
    out = {}
    base_bar = max(results['runs']['supine_unstopped']['base']['worst_excursion_rad'].values())
    for protocol, arms in results['runs'].items():
        base = arms['base']
        rows = {}
        for arm, r in arms.items():
            w = r['worst_excursion_rad']
            common = {k: v for k, v in w.items() if not k.startswith('mtp_angle')}
            ankle = max(w['ankle_angle_r'], w['ankle_angle_l'])
            b_ankle = max(base['worst_excursion_rad']['ankle_angle_r'],
                          base['worst_excursion_rad']['ankle_angle_l'])
            item = {'mtp_extrema_rad': {s: r['extrema_rad'].get('mtp_angle_' + s) for s in SIDES},
                    'mtp_excursion_past_declared_rad': {s: w.get('mtp_angle_' + s) for s in SIDES},
                    'ankle_excursion_rad': {s: w['ankle_angle_' + s] for s in SIDES},
                    'delta_ankle_vs_base_rad': ankle - b_ankle,
                    'ankle_changes_materially': abs(ankle - b_ankle) > Q2_MATERIAL_RAD,
                    'worst_excluding_mtp_rad': max(common.values()),
                    'worst_excluding_mtp_coordinate': max(common, key=common.get),
                    'wall_s_per_advance': r['wall_s_per_advance']}
            if protocol.startswith('crawl'):
                item.update({k: r[k] for k in ('pelvis_forward_travel_m', 'first_crawl_divergence',
                                               'stop_reason', 'failure', 'simulated_s')})
            rows[arm] = item
        out[protocol] = rows
    # Would F2 / G-S flip if their bar were measured on the other arms?
    flips = {}
    for arm in ('base_mtp_welded', 'base_toe_geometry_paths'):
        bar = max(v for k, v in results['runs']['supine_unstopped'][arm]['worst_excursion_rad'].items())
        flips[arm] = {'bar_rad': bar, 'bar_delta_vs_base_rad': bar - base_bar,
                      'F2_thoracic_extension_0.3940_still_over': 0.3940 > bar,
                      'G-S_ankles_1.1399_0.8709_still_over': 1.1399 > bar and 0.8709 > bar}
    out['gate_flips'] = flips
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--part', choices=('q1', 'q2', 'all'), default='all')
    a = ap.parse_args()
    if a.part in ('q1', 'all'):
        run_q1()
    if a.part in ('q2', 'all'):
        run_q2()


if __name__ == '__main__':
    main()
