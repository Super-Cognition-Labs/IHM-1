#!/usr/bin/env python3
"""What one soft-tissue solve costs, change by change, and how far that is from the plant's 10 ms step.

MEASURED, not gated (the gates are S1-S6 and RT1 in scripts/verify_soft_tissue.py).  Each
configuration runs in its OWN child process, one after another, so the peak RSS printed is that
configuration's and nothing else's.  Fixture: calcn_l, fitted surface,
confined reading, in the battery's HEEL FIXTURE (scripts/verify_soft_tissue.py heel_pose: the least
rotation about the frame z axis that keeps the forefoot half clear at 6 mm).  The bundle frame pose is
NOT used: there the skin's lowest point is a seam wedge on the calcn/toes cut, which the fitted layers
do not reach at 1-2 mm, so every solve would be unloaded and cost nothing.

  components  one Newton iteration taken apart at the 2 mm solution: gradient (shipped np.add.at vs
              the cached scatter), Hessian (element_hessians + COO assembly vs a matrix-free
              product), linear solve (SuperLU COLAMD factor, MMD symmetric-mode factor, and a
              triangular solve with the cached rest factor)
  newton      method='newton', cold, at each penetration   (the layer as first built)
  fast        method='fast', cold, at each penetration     (cached factor + Newton-CG)
  ramp        method='fast', warm-started, a pose ramp at a fixed step per 10 ms plant step
  condensed   the SMALL-STRAIN condensed surface compliance (the rest stiffness condensed onto the
              candidate contact nodes, contact as a bound-constrained QP): its cost, and its force
              error against the nonlinear solve on the same ramp.  It is linear, so it cannot pass
              K2/K3 at finite strain and is NOT offered as a layer; it is here to say what that
              reduction would buy and cost.

    PYTHONPATH=. prlimit --as=4294967296 .venv/bin/python -u scripts/measure_soft_tissue_speed.py [--depth local|segment_median]
"""
import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly import soft_tissue_layer as stl        # noqa: E402

HEEL = 'calcn_l'
OUT = ROOT / 'data/derived/soft-tissue-layer-v1/speed.json'


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


ROTATION = np.eye(3)


def fixture(depth):
    global ROTATION
    sys.path.insert(0, str(ROOT / 'scripts'))
    from verify_soft_tissue import heel_pose
    record = [r for r in json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_text())['records'] if r['body'] == HEEL][0]
    ROTATION, low, _ = heel_pose(stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])[0])
    began = time.perf_counter()
    layer = stl.segment_layer(ROOT, HEEL, surface='fitted', depth=depth)
    return layer, low, time.perf_counter() - began


def timed(fn, repeat):
    best = np.inf
    for _ in range(repeat):
        t = time.perf_counter()
        out = fn()
        best = min(best, time.perf_counter() - t)
    return best, out


def components(layer, low, penetration):
    from scipy.sparse.linalg import splu
    from ihm.assembly.mechanics_backend import DeformableRegion
    r = layer.solve(rotation=ROTATION, plane_axis=1, plane_value_m=low + penetration, plane_sign=1, method='fast', return_positions=True)
    free_dof = ~np.repeat(layer.base[:, None], 3, axis=1)
    active = free_dof.any(axis=1)[layer.tets].any(axis=1)
    tets = layer.tets[active]
    nodes = np.unique(tets)
    index = np.full(len(layer.local), -1)
    index[nodes] = np.arange(len(nodes))
    region = DeformableRegion(layer.local[nodes], index[tets], mu_pa=layer.mu[active], lambda_pa=layer.lam[active], density_kg_m3=1.)
    t = time.perf_counter()
    kernel = layer._kernel(active, nodes, index[tets], free_dof[nodes])
    posed = kernel.posed(np.eye(3))
    y = r['positions_m'][nodes]
    out = {'kernel_build_s_first_call_or_cached': time.perf_counter() - t, 'factor_nnz': kernel.factor_nnz,
           'free_dof': int(kernel.free.sum())}
    out['gradient_shipped_s'], _ = timed(lambda: region.energy_gradient(y), 20)
    out['gradient_scatter_s'], _ = timed(lambda: posed.energy_gradient(y), 20)
    out['hessian_elements_s'], h = timed(lambda: stl.element_hessians(region, y, project=False), 5)
    out['hessian_elements_projected_s'], _ = timed(lambda: stl.element_hessians(region, y, project=True), 2)
    out['hessian_assemble_s'], k = timed(lambda: stl.assemble(region, h), 5)
    v = np.random.default_rng(0).normal(size=y.shape)
    out['hessian_prepare_s'], _ = timed(lambda: posed.prepare(y), 20)
    out['hessian_vector_s'], _ = timed(lambda: posed.hessian_vector(v), 20)
    idx = np.flatnonzero(kernel.free)
    kf = k[idx][:, idx].tocsc()
    out['factor_colamd_s'], lu = timed(lambda: splu(kf, permc_spec='COLAMD'), 3)
    out['factor_colamd_nnz'] = int(lu.L.nnz + lu.U.nnz)
    out['factor_mmd_symmetric_s'], lu2 = timed(lambda: splu(kf, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.,
                                                            options=dict(SymmetricMode=True)), 3)
    out['factor_mmd_symmetric_nnz'] = int(lu2.L.nnz + lu2.U.nnz)
    rhs = np.ones(len(idx))
    out['solve_cached_factor_s'], _ = timed(lambda: kernel.factor.solve(rhs), 20)
    out['precondition_rotated_s'], _ = timed(lambda: posed.precondition(v), 20)
    return out


def solves(layer, low, method, depths):
    rows = []
    for d in depths:
        try:
            r = layer.solve(rotation=ROTATION, plane_axis=1, plane_value_m=low + d, plane_sign=1, method=method)
        except (ValueError, RuntimeError) as error:
            rows.append({'penetration_mm': d * 1e3, 'not_computable': str(error)[:200]})
            continue
        rows.append({'penetration_mm': d * 1e3, 'force_n': r['segment_force_n'][1], 'iterations': r['iterations'],
                     'wall_s': r['wall_seconds'], 'converged': r['converged'], 'message': r['solver_message'],
                     'peak_mb': peak_mb()})
    return rows


def ramp(layer, low, top, step, tolerance=None):
    state, rows = None, []
    for k in range(1, int(round(top / step)) + 1):
        t = time.perf_counter()
        r = layer.solve(rotation=ROTATION, plane_axis=1, plane_value_m=low + k * step, plane_sign=1, method='fast',
                        force_tolerance_n=tolerance,
                        warm_start_local_m=None if state is None else state['positions_local_m'])
        rows.append({'penetration_mm': k * step * 1e3, 'wall_s': time.perf_counter() - t, 'force_n': r['segment_force_n'][1],
                     'iterations': r['iterations'], 'converged': r['converged']})
        state = r['state']
    walls = np.array([r['wall_s'] for r in rows])
    return {'step_mm': step * 1e3, 'rows': rows, 'median_ms': float(np.median(walls) * 1e3),
            'max_ms': float(walls.max() * 1e3), 'warm_median_ms': float(np.median(walls[1:]) * 1e3), 'peak_mb': peak_mb()}


def condensed(layer, low, top, step):
    """Small-strain surface compliance: C = K0^-1 on the candidate nodes' vertical DOF, contact as
    min 1/2 l'Cl - p'l, l >= 0 (the linear-elastic Signorini problem), solved by a primal active set
    (Lawson-Hanson's, on C directly; scipy's nnls took 0.23 s on 154 unknowns, which measured scipy and not
    the reduction)."""

    def active_set(c, p):
        lam = np.zeros(len(p))
        active = np.zeros(len(p), bool)
        for _ in range(4 * len(p)):
            w = p - c @ lam
            free = ~active
            if not free.any() or w[free].max() <= 1e-12 * max(1.0, np.abs(p).max()):
                return lam
            active[np.flatnonzero(free)[np.argmax(w[free])]] = True
            while True:
                z = np.zeros(len(p))
                idx = np.flatnonzero(active)
                z[idx] = np.linalg.solve(c[np.ix_(idx, idx)], p[idx])
                if (z[idx] > 0).all():
                    lam = z
                    break
                neg = idx[z[idx] <= 0]
                t = np.min(lam[neg] / (lam[neg] - z[neg]))
                lam = lam + t * (z - lam)
                active &= lam > 1e-15
        raise RuntimeError('active set did not terminate')
    free_dof = ~np.repeat(layer.base[:, None], 3, axis=1)
    active = free_dof.any(axis=1)[layer.tets].any(axis=1)
    tets = layer.tets[active]
    nodes = np.unique(tets)
    index = np.full(len(layer.local), -1)
    index[nodes] = np.arange(len(nodes))
    kernel = layer._kernel(active, nodes, index[tets], free_dof[nodes])
    y = (layer.local[nodes] @ ROTATION.T)[:, 1]
    candidates = np.flatnonzero((y < low + top) & ~layer.base[nodes])
    t = time.perf_counter()
    position = -np.ones(kernel.free.size, int)
    position[np.flatnonzero(kernel.free)] = np.arange(int(kernel.free.sum()))
    # the world-vertical DOF of a node is a direction n = R^T e_y in the segment frame, where K0 lives
    n = ROTATION.T[:, 1]
    cols = []
    rhs_all = np.zeros((len(candidates), int(kernel.free.sum())))
    for j, node in enumerate(candidates):
        rhs_all[j, position[node * 3 + np.arange(3)]] = n
    for j in range(len(candidates)):
        z = kernel.factor.solve(rhs_all[j])
        cols.append(np.array([z[position[m * 3 + np.arange(3)]] @ n for m in candidates]))
    c = np.array(cols).T
    c = 0.5 * (c + c.T)
    setup = time.perf_counter() - t
    rows = []
    for k in range(1, int(round(top / step)) + 1):
        t = time.perf_counter()
        pen = np.maximum(low + k * step - y[candidates], 0.)
        lam = active_set(c, pen)
        rows.append({'penetration_mm': k * step * 1e3, 'force_n': float(lam.sum()), 'wall_s': time.perf_counter() - t})
    return {'candidate_nodes': int(len(candidates)), 'setup_s': setup, 'rows': rows, 'peak_mb': peak_mb()}


def child(mode, depth):
    layer, low, build = fixture(depth)
    out = {'mode': mode, 'depth_rule': depth, 'build_s': build, 'dof': layer.dof, 'tetrahedra': int(len(layer.tets))}
    depths = (0.002, 0.006, 0.012) if depth == 'segment_median' else (0.001, 0.002, 0.004)
    top = 0.012 if depth == 'segment_median' else 0.004
    if mode == 'components':
        out['at'] = components(layer, low, 0.002)
    elif mode in ('newton', 'fast'):
        out['rows'] = solves(layer, low, mode, depths)
    elif mode == 'ramp':
        out['ramps'] = [ramp(layer, low, top, 0.0005), ramp(layer, low, top, 0.002)]
        # the same 0.5 mm ramp at looser force tolerances: what a caller trading accuracy for time would get,
        # with the force error that costs measured against the default tolerance's answer
        out['tolerance_sweep'] = [dict(ramp(layer, low, top, 0.0005, tol), force_tolerance_n=tol) for tol in (1e-6, 1e-4)]
    elif mode == 'condensed':
        out['condensed'] = condensed(layer, low, top, 0.0005)
        out['nonlinear'] = ramp(layer, low, top, 0.0005)
    out['peak_mb'] = peak_mb()
    print('RESULT ' + json.dumps(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--child', default=None)
    ap.add_argument('--depth', default='local', choices=stl.DEPTH_RULES)
    args = ap.parse_args()
    if args.child:
        child(args.child, args.depth)
        return
    results = {}
    for mode in ('components', 'newton', 'fast', 'ramp', 'condensed'):
        p = subprocess.run([sys.executable, '-u', __file__, '--child', mode, '--depth', args.depth],
                           capture_output=True, text=True, cwd=ROOT)
        line = [ln for ln in p.stdout.splitlines() if ln.startswith('RESULT ')]
        if p.returncode or not line:
            print(f'{mode}: FAILED\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}', flush=True)
            results[mode] = {'failed': p.stderr[-2000:]}
            continue
        results[mode] = json.loads(line[0][7:])
        print(f'{mode}: ' + json.dumps(results[mode], indent=None)[:3000], flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out = OUT.with_name(f'speed-{args.depth}.json')
    out.write_text(json.dumps(results, indent=1) + '\n')
    print('wrote', out.relative_to(ROOT))


if __name__ == '__main__':
    main()
