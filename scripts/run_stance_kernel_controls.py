#!/usr/bin/env python3
"""Run docs/STANCE_KERNEL_CONTROLS.md: eight kernels x twenty-four pushes, then the statistic.

  --phase build    build the new control kernels and run the stance pipeline on each
  --phase eval     every policy x condition through evaluate_native_cortical_stance.py
  --phase collect  read every report back and write out/stance_kernel_controls.json

Each step is a subprocess writing its own directory, so a phase can be re-run and skips
what already exists. Parallelism is across subprocesses, one native body each.
"""
import argparse, hashlib, json, math, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PY = str(ROOT / '.venv/bin/python')
WORK = ROOT / 'data/runtime/motor-learning/stance-kernel-controls'
REG = 'data/models/engineering_stance_v1/registration.json'
TEACHER = 'data/research/locomotion_control/delay_design_e5s_lepu/candidate_state1.0_R1000.0.npz'
EXTRA = 'data/research/locomotion_control/patient_slow_lqr_delay4_com_push/trajectory.json'

EXISTING = {'fused': 'data/models/ibm_cortical_stance_v1/cortical_stance.pt',
            'c16': 'data/models/ibm_curriculum16_stance_v1/cortical_stance.pt',
            'perm-20260908': 'data/models/ibm_curriculum16_shuffled_stance_v1/cortical_stance.pt'}
NEW = {'perm-1': ('permute', 1), 'perm-2': ('permute', 2),
       'rand-11': ('random', 11), 'rand-12': ('random', 12), 'rand-13': ('random', 13)}
# Instrument: the existing permuted kernel pushed through THIS host's pipeline again.
REBUILD = {'perm-20260908-rebuilt': 'data/models/ibm_curriculum16_shuffled_kernel_v1/kernel.pt'}
LEARNED = ('fused', 'c16')
CONTROLS = ('perm-20260908', 'perm-1', 'perm-2', 'rand-11', 'rand-12', 'rand-13')


def conditions():
    out = []
    for x, z in ((5, 0), (-5, 0), (0, 5), (0, -5)): out.append(('A', x, z))
    for x, z in ((5, 5), (5, -5), (-5, 5), (-5, -5)): out.append(('B', x, z))
    for mag, tag in ((10, 'C'), (20, 'D')):
        for k in range(8):
            t = math.radians(45 * k)
            out.append((tag, round(mag * math.cos(t), 6) + 0.0, round(mag * math.sin(t), 6) + 0.0))
    return out


def cid(c):
    return f'{c[0]}_{c[1]:+g}_{c[2]:+g}'


def run(cmd, log):
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open('w') as f:
        r = subprocess.run(cmd, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT, env={'PYTHONPATH': str(ROOT), 'PATH': '/usr/bin:/bin', 'HOME': str(Path.home()), 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    if r.returncode: raise RuntimeError(f'{cmd[1]} failed; see {log}')


def pipeline(name, kernel):
    base = WORK / 'policies' / name
    final = base / 'stance'
    if (final / 'cortical_stance.pt').exists(): return str((final / 'cortical_stance.pt').relative_to(ROOT))
    s1, s2 = base / 'stage1', base / 'stage2'
    if not (s1 / 'cortical_stance.pt').exists():
        run([PY, 'scripts/train_cortical_stance_sequences.py', '--sites', '1024', '--encoder-kind', 'signed_identity',
             '--epochs', '100', '--lr', '0.001', '--reserve-weight', '--teacher', TEACHER, '--extra-trajectory', EXTRA,
             '--checkpoint', kernel, '--out', str(s1.relative_to(ROOT))], base / 'stage1.log')
    if not (s2 / 'cortical_stance.pt').exists():
        run([PY, 'scripts/fit_cortical_stance_local_readout.py', '--artifact', str((s1 / 'cortical_stance.pt').relative_to(ROOT)),
             '--output', str(s2.relative_to(ROOT))], base / 'stage2.log')
    run([PY, 'scripts/scale_cortical_stance_regime.py', '--artifact', str((s2 / 'cortical_stance.pt').relative_to(ROOT)),
         '--output', str(final.relative_to(ROOT)), '--scale', '100'], base / 'stage3.log')
    return str((final / 'cortical_stance.pt').relative_to(ROOT))


def policies():
    out = dict(EXISTING)
    for name in list(NEW) + list(REBUILD):
        p = WORK / 'policies' / name / 'stance' / 'cortical_stance.pt'
        if p.exists(): out[name] = str(p.relative_to(ROOT))
    return out


def build(jobs):
    tasks = {}
    for name, (mode, seed) in NEW.items():
        kdir = WORK / 'kernels' / name
        if not (kdir / 'kernel.pt').exists():
            run([PY, 'scripts/build_control_kernel.py', '--mode', mode, '--seed', str(seed), '--output', str(kdir.relative_to(ROOT))], WORK / 'kernels' / f'{name}.log')
        tasks[name] = str((kdir / 'kernel.pt').relative_to(ROOT))
    tasks.update(REBUILD)
    with ThreadPoolExecutor(jobs) as pool:
        futures = {n: pool.submit(pipeline, n, k) for n, k in tasks.items()}
        for n, f in futures.items(): print(n, f.result(), flush=True)


def evaluate(jobs, only=None):
    pols = policies(); work = []
    for name, artifact in pols.items():
        if only and name not in only: continue
        for c in conditions():
            out = WORK / 'evals' / name / cid(c)
            if (out / 'report.json').exists(): continue
            if out.exists():
                import shutil; shutil.rmtree(out)
            work.append((name, artifact, c, out))
    print(f'{len(work)} evaluations on {jobs} workers', flush=True)
    t0 = time.time()

    def one(item):
        name, artifact, c, out = item
        run([PY, 'scripts/evaluate_native_cortical_stance.py', '--artifact', artifact, '--registration', REG,
             '--seconds', '10', '--perturb-x', str(c[1]), '--perturb-z', str(c[2]), '--output', str(out.relative_to(ROOT))],
            out.parent / (out.name + '.log'))
        # Big per-frame trajectories stay out of the way; the report holds what is scored.
        for big in ('trajectory.json', 'partial.json'):
            p = out / big
            if p.exists(): p.unlink()
        return name, cid(c)

    with ThreadPoolExecutor(jobs) as pool:
        for i, (n, c) in enumerate(pool.map(one, work)):
            print(f'[{i + 1}/{len(work)} {time.time() - t0:.0f}s] {n} {c}', flush=True)


def collect():
    pols = policies(); rows = {}
    for name in pols:
        for c in conditions():
            rp = WORK / 'evals' / name / cid(c) / 'report.json'
            if not rp.exists(): continue
            r = json.loads(rp.read_text())
            full, sever = r['arms']['full'], r['arms']['sever']
            rows.setdefault(name, {})[cid(c)] = dict(
                set=c[0], force_n=r['perturbation_force_n'], artifact_sha256=r['artifact_sha256'],
                report_sha256=hashlib.sha256(rp.read_bytes()).hexdigest(),
                full=dict(completed=full['completed_horizon'], error=full['error'], peak_m=full['peak_com_displacement_m'],
                          final_m=full.get('final_com_displacement_m'), elapsed_s=full.get('elapsed_s')),
                sever=dict(completed=sever.get('completed_horizon'), peak_m=sever.get('peak_com_displacement_m'), elapsed_s=sever.get('elapsed_s'), error=sever.get('error')))
    conds = [cid(c) for c in conditions()]
    complete = [c for c in conds if all(n in rows and c in rows[n] and rows[n][c]['full']['completed'] for n in LEARNED + CONTROLS)]
    delta = []
    for c in complete:
        L = np.mean([math.log(rows[n][c]['full']['peak_m']) for n in LEARNED])
        K = [math.log(rows[n][c]['full']['peak_m']) for n in CONTROLS]
        delta.append(dict(condition=c, delta=L - np.mean(K), control_range=max(K) - min(K),
                          perm_mean=np.mean(K[:3]) - L, rand_mean=np.mean(K[3:]) - L))
    d = np.array([x['delta'] for x in delta])
    rng = np.random.default_rng(0)
    boot = np.array([rng.choice(d, len(d)).mean() for _ in range(10000)]) if len(d) else np.array([np.nan])
    lo, hi = np.percentile(boot, [2.5, 97.5]) if len(d) else (np.nan, np.nan)
    mean = float(d.mean()) if len(d) else float('nan')
    if len(d) == 0: verdict = 'no complete conditions'
    elif mean <= -0.05 and hi < 0: verdict = 'learned content HELPS (>=5%)'
    elif -0.05 < lo and hi < 0.05: verdict = 'EQUIVALENCE at the 5% margin'
    elif mean >= 0.05 and lo > 0: verdict = 'learned kernels WORSE'
    else: verdict = 'inconclusive at this n'
    recovered = {n: sum(1 for c in conds if c in rows.get(n, {}) and rows[n][c]['full']['completed']) for n in rows}
    learned_more = all(recovered.get(l, 0) > recovered.get(k, 0) for l in LEARNED for k in CONTROLS)
    out = dict(schema='ihm.stance-kernel-controls.v1', preregistration='docs/STANCE_KERNEL_CONTROLS.md',
               policies=pols, learned=LEARNED, controls=CONTROLS, conditions=conds,
               complete_conditions=len(complete), mean_delta_ln_peak=mean, ci95=[float(lo), float(hi)],
               verdict=verdict, recovered_of=len(conds), recovered=recovered, learned_recover_more_than_every_control=learned_more,
               median_control_range_ln=float(np.median([x['control_range'] for x in delta])) if delta else None,
               per_condition=delta, rows=rows)
    path = ROOT / 'out/stance_kernel_controls.json'
    path.write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps({k: out[k] for k in ('complete_conditions', 'mean_delta_ln_peak', 'ci95', 'verdict', 'recovered', 'median_control_range_ln')}, indent=1, default=float))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--phase', choices=('build', 'eval', 'collect'), required=True)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--only', nargs='*')
    a = ap.parse_args()
    {'build': lambda: build(a.jobs), 'eval': lambda: evaluate(a.jobs, a.only), 'collect': collect}[a.phase]()


if __name__ == '__main__':
    main()
