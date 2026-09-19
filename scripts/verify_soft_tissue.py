#!/usr/bin/env python3
"""Known answers for the deformable soft-tissue layer (ihm/assembly/soft_tissue_layer.py).

Every bar below was written into this header BEFORE this file was first run, except
the two marked MEASURED and the one marked RECORDED, which say why.

KNOWN ANSWERS -- each has a value fixed by something other than this solver
  K1  The Hessian is the derivative of the SHIPPED gradient (DeformableRegion.energy_gradient):
      central differences over a perturbed 2x2x2 box, max |H - H_fd| / max |H| <= 1e-7.
  K2  Confined uniaxial compression IS the shipped foundation law.  Lateral DOF held, a
      frictionless plane pushed in 10/30/50 %.  Force must equal
      `supine_contact.foundation` (called, not re-typed) on the same area: rel <= 1e-8.
  K3  Unconfined uniaxial compression IS compressible neo-Hookean with zero lateral stress:
      solve mu b + (lam ln(a b^2) - mu)/b = 0 for b, force = [mu a + (lam ln(a b^2) - mu)/a] A0,
      rel <= 1e-8 at the same three strains.
  K4  THE NULL CASE (the input on which the extra freedom must produce nothing): a support
      that does not reach the layer gives EXACTLY zero force, zero displacement, and no solve.
  K5  Load reaches the segment.  Heel pressed 2 mm into a floor: force on the segment equals
      the contact force, and moments agree, to |residual| / |contact| <= 1e-6.
  K6  CALL IT TWICE (CLAUDE.md).  The same solve twice, box and heel, and the same mesh
      built twice: bitwise identical.
  K7  Energy is not created.  A compressed box released with backward Euler: kinetic +
      elastic energy non-increasing at every one of 40 steps (tolerance 1e-9 of the start).
  K9  Two implementations agree.  Projected Newton against DeformableRegion's own L-BFGS-B
      on the unconfined box: rel <= 1e-6.
 K11  The unanchored list in plant_options is re-derived: exactly the segments whose layer
      has no core cell at 5 mm.
 K12  plant_options: `None` is the historical plant (kwargs {}); adding `soft_tissue`
      never changes the plant's kwargs; a tampered selection cannot build a layer; the
      selection is JSON (articulated.py writes it to disk).

CONTROLS THAT MUST FAIL -- a gate whose pass looks like its absence is not evidence
  C1  Lame parameters swapped: K2 must FAIL.
  C2  A solve starved to one iteration: the K5 balance gate must FAIL.
      (Run 3 recorded C2 FAILED: the starved solve read 0.00e+00 because the relative
      residual divided by a contact force that was zero and defaulted to 0.  The metric was
      the defect -- it now divides by the largest force in the balance -- and the bar is
      unchanged.  Run 1 recorded K5 FAILED at 0 N: a failed Newton solve reported its
      START, not its last iterate.  Both logs are kept in logs/.)
  C3  An integrator KNOWN to create energy (explicit symplectic Euler above its stability
      limit, same state and dt as K7): the K7 gate must FAIL.

MEASURED -- reported, not gated, because the answer is not known in advance
  M1  Cost: DOF, peak memory, wall seconds per solve, Newton iterations, heel at 2/6/12 mm.
  M2  The heel against the shipped 1-D models at the SAME physical penetration (measured from
      the skin mesh's lowest point, not the voxel surface): linear k = E/h (what the
      skin_layer_map contact bundle does), the confined neo-Hookean column at the declared h
      and at the voxel model's own column depth, and the 3-D layer under both modulus readings.
      MEASURED, not pre-registered: development runs had already shown the 3-D heel is roughly
      an order of magnitude softer, so no bar is set on it here.
  M3  Mesh convergence of the heel force at 6 mm over spacing 6/5/4 mm.

RECORDED
  R1  The Saint-Venant crop must put <= 1% of the load on the truncation at radius 3h, heel,
      6 mm.  This bar was written AFTER development runs showed 17-100% on this segment, so it
      records an observed failure rather than testing a prediction; it is printed as FAILED
      if it fails and excluded from the exit status for that reason, and only that one.

    PYTHONPATH=. prlimit --as=4294967296 .venv/bin/python -u scripts/verify_soft_tissue.py

Memory budget: 4 GiB address space (prlimit), measured peak is printed at the end.
One heavy job at a time.
"""
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.assembly import soft_tissue_layer as stl                      # noqa: E402
from ihm.assembly.mechanics_backend import DeformableRegion, tetra_box  # noqa: E402
from ihm.assembly.supine_contact import foundation                      # noqa: E402

REPORT = ROOT / 'data/derived/soft-tissue-layer-v1/report.json'
RECORDED_FAILURES = {'R1'}
HEEL = 'calcn_l'
results, failures = {}, []


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def gate(key, name, ok, detail=''):
    status = 'PASS' if ok else ('FAILED (recorded)' if key in RECORDED_FAILURES else 'FAIL')
    print(f'  {status:17s} {key:4s} {name}{"  " + detail if detail else ""}', flush=True)
    results[key] = {'name': name, 'pass': bool(ok), 'detail': detail}
    if not ok and key not in RECORDED_FAILURES:
        failures.append(key)
    return bool(ok)


# ---------------------------------------------------------------------------- fixtures
LX, LY, H = 0.02, 0.02, 0.01
E_HEEL, NU = 192550.0, 0.45            # the heel's own layer-map modulus and the bundle's Poisson ratio
MU, LAM, _ = stl.lame_from_layer_modulus(E_HEEL, NU, 'confined')


def box(mu=MU, lam=LAM, base='bottom'):
    x, t = tetra_box([LX, LY, H], (4, 4, 4))
    bottom = np.isclose(x[:, 2], 0)
    anchor = bottom if base == 'bottom' else np.all(np.isclose(x, 0), axis=1)
    layer = stl.SoftTissueLayer({'nodes_m': x, 'tetrahedra': t, 'base': anchor, 'spacing_m': H / 4},
                                mu_pa=mu, lambda_pa=lam, density_kg_m3=950.)
    return layer, x, bottom


def confined_held(x):
    held = np.zeros_like(x, bool)
    held[:, :2] = True
    return held


def roller_held(x, bottom):
    held = np.zeros_like(x, bool)
    held[bottom, 2] = True
    held[np.isclose(x[:, 1], 0) & bottom, 1] = True
    return held


def shipped_confined_force(stretch, mu, lam):
    """supine_contact.foundation itself, on one quadrature point of area A0."""
    thickness = 1.0
    penetration = (1 - stretch) * thickness
    out = foundation(np.array([[0., 0., 0.]]), np.zeros((1, 3)), np.array([LX * LY]), np.array([0]),
                     np.zeros((1, 3)), penetration,
                     dict(total_layer_thickness_m=thickness, shear_modulus_pa=mu, lame_lambda_pa=lam,
                          minimum_thickness_ratio=0.2, dissipation_s_m=0., dynamic_friction=0.,
                          transition_velocity_m_s=0.1, viscous_friction=0.))
    return float(out['body_forces_n'][0][0])


def unconfined_force(a, mu, lam):
    b = brentq(lambda b: mu * b + (lam * np.log(a * b * b) - mu) / b, 0.3, 5.0)
    return -(mu * a + (lam * np.log(a * b * b) - mu) / a) * LX * LY, b


# ----------------------------------------------------------------------------- battery
def k1_hessian():
    print('\nK1  the Hessian is the derivative of the shipped gradient', flush=True)
    from scipy import sparse
    x, t = tetra_box([LX, LY, H], (2, 2, 2))
    region = DeformableRegion(x, t, mu_pa=MU, lambda_pa=LAM, density_kg_m3=1.)
    y = x + np.random.default_rng(0).normal(scale=4e-4, size=x.shape)
    k = stl.assemble(region, stl.element_hessians(region, y, project=False)).toarray()
    fd = np.zeros_like(k)
    eps = 1e-8
    for i in range(x.size):
        d = np.zeros(x.size)
        d[i] = eps
        fd[:, i] = (region.energy_gradient(y + d.reshape(-1, 3))[1].ravel()
                    - region.energy_gradient(y - d.reshape(-1, 3))[1].ravel()) / (2 * eps)
    err = float(np.abs(k - fd).max() / np.abs(k).max())
    gate('K1', 'Hessian vs central differences of energy_gradient', err <= 1e-7, f'rel {err:.2e}')
    return err


def k2_confined(mu=MU, lam=LAM, key='K2', label='confined compression = supine_contact.foundation'):
    layer, x, _ = box(mu, lam)
    worst = 0.0
    for strain in (0.1, 0.3, 0.5):
        r = layer.solve(plane_axis=2, plane_value_m=H * (1 - strain), plane_sign=-1, held=confined_held(x))
        want = shipped_confined_force(1 - strain, MU, LAM)     # ALWAYS the true material
        worst = max(worst, abs(-r['segment_force_n'][2] - want) / want)
    return gate(key, label, worst <= 1e-8, f'worst rel {worst:.2e} over strain 0.1/0.3/0.5'), worst


def k3_unconfined():
    print('\nK3  unconfined compression = compressible neo-Hookean, zero lateral stress', flush=True)
    layer, x, bottom = box(base='corner')
    worst, rows = 0.0, []
    for strain in (0.1, 0.3, 0.5):
        r = layer.solve(plane_axis=2, plane_value_m=H * (1 - strain), plane_sign=-1,
                        held=roller_held(x, bottom), return_positions=True)
        want, b = unconfined_force(1 - strain, MU, LAM)
        got = -r['segment_force_n'][2]
        lateral = float(r['positions_m'][:, 0].max() / LX)
        worst = max(worst, abs(got - want) / want, abs(lateral - b) / b)
        rows.append((strain, got, want, lateral, b, r['iterations']))
    for s, got, want, lateral, b, it in rows:
        print(f'      strain {s:.1f}: F {got:.10f} N vs {want:.10f}; lateral stretch {lateral:.8f} vs {b:.8f}; {it} it')
    gate('K3', 'unconfined force and lateral stretch vs analytic', worst <= 1e-8, f'worst rel {worst:.2e}')
    confined = shipped_confined_force(0.7, MU, LAM)
    print(f'      at 30% strain the tissue free to bulge carries {rows[1][2] / confined:.3f} of the confined force')
    return rows


def k4_null(heel):
    print('\nK4  the null case: a support that does not reach the layer', flush=True)
    low = heel.local[:, 1].min()
    r = heel.solve(plane_axis=1, plane_value_m=low - 0.001, plane_sign=1, return_positions=True)
    zero = (not any(r['segment_force_n']) and not any(r['segment_moment_nm'])
            and r['maximum_displacement_m'] == 0.0 and r['iterations'] == 0)
    gate('K4', 'exactly zero force, moment and displacement, no solve', zero,
         f'{r["solver_message"]}; F {r["segment_force_n"]}')
    touching = heel.solve(plane_axis=1, plane_value_m=low, plane_sign=1)
    gate('K4b', 'a support exactly touching carries exactly nothing',
         not any(touching['segment_force_n']), f'F {touching["segment_force_n"]}')


def k5_transmission(heel, mesh_low):
    print('\nK5  load reaches the segment (heel, 2 mm into the floor)', flush=True)
    r = heel.solve(plane_axis=1, plane_value_m=mesh_low + 0.002, plane_sign=1)
    contact = np.linalg.norm(r['contact_force_on_tissue_n'])
    moment = np.linalg.norm(r['balance_moment_residual_nm']) / max(
        np.linalg.norm(r['contact_moment_about_segment_nm']), 1e-300)
    gate('K5', 'segment force = contact force', r['balance_force_relative'] <= 1e-6 and contact > 0,
         f'contact {contact:.6f} N, relative residual {r["balance_force_relative"]:.2e}')
    gate('K5b', 'segment moment = contact moment', moment <= 1e-6, f'relative residual {moment:.2e}')
    return r


def k6_idempotent(heel, mesh_low, first_heel):
    print('\nK6  call it twice at the same input', flush=True)
    layer, x, bottom = box(base='corner')
    a = layer.solve(plane_axis=2, plane_value_m=0.007, plane_sign=-1, held=roller_held(x, bottom), return_positions=True)
    b = layer.solve(plane_axis=2, plane_value_m=0.007, plane_sign=-1, held=roller_held(x, bottom), return_positions=True)
    same_box = np.array_equal(a['positions_m'], b['positions_m']) and a['segment_force_n'] == b['segment_force_n']
    gate('K6', 'box solve twice: bitwise identical', same_box)
    again = heel.solve(plane_axis=1, plane_value_m=mesh_low + 0.002, plane_sign=1)
    gate('K6b', 'heel solve twice: bitwise identical',
         again['segment_force_n'] == first_heel['segment_force_n']
         and again['segment_moment_nm'] == first_heel['segment_moment_nm'])
    record = [r for r in json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_text())['records'] if r['body'] == HEEL][0]
    v, f = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
    m1 = stl.layer_mesh(v, f, record['layer']['thickness_m'])
    m2 = stl.layer_mesh(v, f, record['layer']['thickness_m'])
    gate('K6c', 'heel mesh built twice: bitwise identical',
         all(np.array_equal(m1[k], m2[k]) for k in ('nodes_m', 'tetrahedra', 'base')))


def released_box():
    """A box pressed 30% by a plane, then released: the starting state of K7 and C3."""
    layer, x, _ = box()
    pressed = layer.solve(plane_axis=2, plane_value_m=H * 0.7, plane_sign=-1, return_positions=True)
    return layer, {'positions_m': pressed['positions_m'], 'velocities_m_s': np.zeros_like(x)}


def total_energy(layer, state):
    e, _ = layer.energy_gradient(state['positions_m'])
    return e + 0.5 * float(np.sum(layer.nodal_mass[:, None] * state['velocities_m_s'] ** 2))


def k7_energy(dt=1e-3, steps=40):
    print('\nK7  energy is not created (backward Euler, released box)', flush=True)
    layer, state = released_box()
    energies = [total_energy(layer, state)]
    for _ in range(steps):
        r = layer.solve(dt_s=dt, state=state)
        state = r['state']
        energies.append(total_energy(layer, state))
    rises = np.diff(energies)
    worst = float(rises.max() / energies[0])
    gate('K7', 'kinetic + elastic non-increasing at every step', worst <= 1e-9,
         f'E0 {energies[0]:.6e} J -> {energies[-1]:.6e} J; worst step rise {worst:+.2e} of E0')
    return energies


def c3_energy_creating(dt=1e-3, steps=40):
    layer, state = released_box()
    x, v = state['positions_m'].copy(), state['velocities_m_s'].copy()
    free = ~layer.base
    energies = [total_energy(layer, state)]
    for _ in range(steps):
        try:
            _, g = layer.energy_gradient(x)
        except ValueError:
            energies.append(np.inf)          # inverted: it has already exploded
            break
        v[free] -= dt * g[free] / layer.nodal_mass[free, None]
        x[free] += dt * v[free]
        if np.linalg.det(layer.deformation(x)).min() <= 0:
            energies.append(np.inf)          # an element inverted: the energy has diverged
            break
        energies.append(total_energy(layer, {'positions_m': x, 'velocities_m_s': v}))
    worst = float(np.max(np.diff(energies)) / energies[0])
    caught = not (worst <= 1e-9)
    gate('C3', 'the K7 gate FAILS an integrator known to create energy', caught,
         f'explicit symplectic Euler, same state and dt: worst step rise {worst:+.2e} of E0')


def k9_cross():
    print('\nK9  projected Newton vs DeformableRegion\'s own L-BFGS-B', flush=True)
    layer, x, bottom = box(base='corner')
    kw = dict(plane_axis=2, plane_value_m=H * 0.7, plane_sign=-1, held=roller_held(x, bottom))
    n = layer.solve(method='newton', **kw)
    b = layer.solve(method='lbfgsb', **kw)
    rel = abs(n['segment_force_n'][2] - b['segment_force_n'][2]) / abs(n['segment_force_n'][2])
    gate('K9', 'same force from two minimisers', rel <= 1e-6,
         f'rel {rel:.2e}; Newton {n["iterations"]} it {n["wall_seconds"]:.2f}s, '
         f'L-BFGS-B {b["iterations"]} it {b["wall_seconds"]:.2f}s')
    return {'newton': [n['iterations'], n['wall_seconds']], 'lbfgsb': [b['iterations'], b['wall_seconds']]}


def controls():
    print('\nC1/C2  controls that must fail', flush=True)
    ok_swapped, worst = k2_confined(mu=LAM, lam=MU, key='C1raw', label='(swapped Lame run, expected to fail K2)')
    results.pop('C1raw')
    if 'C1raw' in failures:
        failures.remove('C1raw')      # the swapped run is SUPPOSED to fail; C1 is the gate on that
    gate('C1', 'swapping the Lame parameters FAILS K2', not ok_swapped, f'worst rel {worst:.2e}')
    layer, x, bottom = box(base='corner')
    starved = layer.solve(plane_axis=2, plane_value_m=H * 0.7, plane_sign=-1, held=roller_held(x, bottom), maxiter=1)
    gate('C2', 'a solve starved to one iteration FAILS the K5 balance gate',
         not (starved['balance_force_relative'] <= 1e-6),
         f'relative residual {starved["balance_force_relative"]:.2e}, converged={starved["converged"]}')


def k11_unanchored():
    print('\nK11  which segments have no rigid core at 5 mm (re-derived)', flush=True)
    from ihm.assembly.plant_options import SOFT_TISSUE_UNANCHORED
    manifest = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_text())
    census, empty = {}, []
    for record in manifest['records']:
        began = time.perf_counter()
        v, f = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
        mesh = stl.layer_mesh(v, f, record['layer']['thickness_m'], name=record['body'])
        n_base = int(mesh['base'].sum())
        census[record['body']] = {
            'thickness_mm': round(record['layer']['thickness_m'] * 1e3, 2),
            'cells_through_thickness': round(mesh['cells_through_thickness'], 2),
            'nodes': int(len(mesh['nodes_m'])), 'tetrahedra': int(len(mesh['tetrahedra'])),
            'base_nodes': n_base, 'dof': int((~mesh['base']).sum()) * 3 if n_base else 0,
            'layer_volume_l': round(mesh['layer_volume_m3'] * 1e3, 3),
            'unanchored_islands': mesh['unanchored_islands'],
            'unanchored_island_tetrahedra': mesh['unanchored_island_tetrahedra'],
            'build_seconds': round(time.perf_counter() - began, 2)}
        if n_base == 0:
            empty.append(record['body'])
        c = census[record['body']]
        print(f'      {record["body"]:10s} h {c["thickness_mm"]:5.1f} mm ({c["cells_through_thickness"]:.1f} cells) '
              f'DOF {c["dof"]:7d} tets {c["tetrahedra"]:7d} base {n_base:6d} islands {c["unanchored_islands"]} '
              f'({c["unanchored_island_tetrahedra"]} tets) build {c["build_seconds"]:.1f}s', flush=True)
    gate('K11', 'unanchored segments match plant_options.SOFT_TISSUE_UNANCHORED',
         sorted(empty) == sorted(SOFT_TISSUE_UNANCHORED), f'derived {sorted(empty)}')
    total = {'dof': sum(c['dof'] for c in census.values()),
             'tetrahedra': sum(c['tetrahedra'] for c in census.values()),
             'layer_volume_l': round(sum(c['layer_volume_l'] for c in census.values()), 3)}
    print(f'      whole body at 5 mm: {total["dof"]} DOF, {total["tetrahedra"]} tetrahedra, '
          f'{total["layer_volume_l"]} L of layer', flush=True)
    return census, total


def k12_plant_options():
    print('\nK12  plant_options wiring', flush=True)
    from ihm.assembly.plant_options import resolve_fidelity, SOFT_TISSUE_LAYERS
    kwargs, selection = resolve_fidelity(ROOT, None)
    gate('K12', 'None is the historical plant', kwargs == {} and selection['soft_tissue'] is None)
    same = True
    for base in ({}, {'joint_stops': True}, {'joint_stops': 'firm', 'tissue_ligaments': 'admissible'}):
        for env in ('upright', 'supine'):
            a, _ = resolve_fidelity(ROOT, dict(base) or None, environment=env) if base else ({}, None)
            b, sel = resolve_fidelity(ROOT, {**base, 'soft_tissue': 'layer_map_confined'}, environment=env)
            same &= json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)
            json.dumps(sel)
    gate('K12b', 'soft_tissue never changes the plant kwargs, and the selection is JSON', same)
    refused = 0
    for bad in ({'soft_tissue': 'data/derived/segment-contact-meshes/skin'}, {'soft_tissue': True}):
        try:
            resolve_fidelity(ROOT, bad, environment='upright')
        except ValueError:
            refused += 1
    try:
        resolve_fidelity(ROOT, {'soft_tissue': 'layer_map_confined'}, environment='free')
    except ValueError:
        refused += 1
    _, sel = resolve_fidelity(ROOT, {'soft_tissue': 'layer_map_confined'}, environment='upright')
    tampered = json.loads(json.dumps(sel))
    tampered['soft_tissue']['mapping'] = 'unconfined'
    try:
        stl.build_selected_layers(ROOT, tampered, [HEEL])
    except ValueError:
        refused += 1
    try:
        stl.build_selected_layers(ROOT, sel, ['patella_l'])
    except ValueError:
        refused += 1
    gate('K12c', 'path-as-id, boolean, free environment, tampered mapping, unanchored segment: all refused',
         refused == 5, f'{refused}/5 refused')
    built = stl.build_selected_layers(ROOT, sel, [HEEL])[HEEL]
    record = [r for r in json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_text())['records'] if r['body'] == HEEL][0]
    gate('K12d', 'the built layer carries the bundle\'s own measured depth and modulus',
         built.meta['thickness_m'] == record['layer']['thickness_m']
         and built.meta['apparent_modulus_pa'] == record['layer']['apparent_modulus_pa']
         and built.meta['mapping'] == SOFT_TISSUE_LAYERS['layer_map_confined']['mapping'])


def m_heel(heel, heel_unconfined, mesh, record):
    print('\nM1/M2  the heel: cost, and the shipped 1-D models at the same physical penetration', flush=True)
    v, f = mesh
    tri = v[f]
    centroid = tri.mean(axis=1)
    area = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    low = float(v[:, 1].min())
    voxel_low = float(heel.local[:, 1].min())
    column = float(heel.local[heel.base, 1].min() - voxel_low)
    h, k = record['layer']['thickness_m'], record['layer']['stiffness_pa_per_m']
    print(f'      skin mesh lowest y {low*1e3:.2f} mm; voxel surface lowest {voxel_low*1e3:.2f} mm '
          f'(staircase offset {(voxel_low-low)*1e3:+.2f} mm); voxel tissue column under the heel '
          f'{column*1e3:.1f} mm against the declared {h*1e3:.1f} mm')

    def column_force(depth, thickness):
        pen = np.maximum(0., depth - (centroid[:, 1] - low))
        s = 1 - pen / thickness
        ok = s > 0.2
        return float((stl.confined_column_pressure(np.where(ok, s, 1.), MU, LAM) * area * ok).sum()), int((~ok & (pen > 0)).sum())

    rows = []
    for depth in (0.002, 0.006, 0.012):
        began_mb = peak_mb()
        r = heel.solve(plane_axis=1, plane_value_m=low + depth, plane_sign=1)
        u = heel_unconfined.solve(plane_axis=1, plane_value_m=low + depth, plane_sign=1)
        linear = float((k * np.maximum(0., depth - (centroid[:, 1] - low)) * area).sum())
        declared, over_d = column_force(depth, h)
        own, over_o = column_force(depth, column)
        row = {'penetration_mm': depth * 1e3, 'layer_3d_confined_n': r['segment_force_n'][1],
               'layer_3d_unconfined_n': u['segment_force_n'][1], 'linear_k_e_over_h_n': linear,
               'nh_column_declared_h_n': declared, 'nh_column_voxel_depth_n': own,
               'columns_past_domain': [over_d, over_o],
               'contact_nodes': r['contact_nodes'], 'solved_dof': r['solved_dof'],
               'iterations': r['iterations'], 'wall_seconds': round(r['wall_seconds'], 2),
               'wall_seconds_unconfined': round(u['wall_seconds'], 2),
               'minimum_jacobian': r['minimum_jacobian'], 'balance_relative': r['balance_force_relative'],
               'message': r['solver_message'], 'peak_mb': round(max(began_mb, peak_mb()), 1)}
        rows.append(row)
        print(f'      {depth*1e3:4.1f} mm | 3-D confined {row["layer_3d_confined_n"]:8.3f} N '
              f'({row["iterations"]} it, {row["wall_seconds"]:.1f} s) | 3-D unconfined {row["layer_3d_unconfined_n"]:8.3f} N '
              f'({row["wall_seconds_unconfined"]:.1f} s) | 1-D linear k=E/h {linear:8.2f} N | '
              f'NH column h={h*1e3:.1f} {declared:8.2f} N | NH column h={column*1e3:.1f} {own:8.2f} N', flush=True)
    return {'mesh_lowest_y_m': low, 'voxel_lowest_y_m': voxel_low, 'voxel_column_m': column,
            'declared_thickness_m': h, 'rows': rows}


def m3_convergence(record, mesh):
    print('\nM3  mesh convergence, heel at 6 mm', flush=True)
    v, f = mesh
    low = float(v[:, 1].min())
    out = []
    for spacing in (0.006, 0.005, 0.004):
        layer = stl.segment_layer(ROOT, HEEL, spacing_m=spacing)
        r = layer.solve(plane_axis=1, plane_value_m=low + 0.006, plane_sign=1)
        out.append({'spacing_mm': spacing * 1e3, 'dof': layer.dof, 'force_n': r['segment_force_n'][1],
                    'voxel_lowest_y_mm': float(layer.local[:, 1].min() * 1e3),
                    'wall_seconds': round(r['wall_seconds'], 2), 'iterations': r['iterations']})
        print(f'      spacing {spacing*1e3:.0f} mm: DOF {layer.dof:6d}  F {r["segment_force_n"][1]:8.3f} N  '
              f'voxel surface {out[-1]["voxel_lowest_y_mm"]:.2f} mm  {r["iterations"]} it  {r["wall_seconds"]:.1f} s', flush=True)
    return out


def r1_crop(heel, low):
    print('\nR1  the Saint-Venant crop (bar written after development runs; see header)', flush=True)
    h = heel.meta['thickness_m']
    full = heel.solve(plane_axis=1, plane_value_m=low + 0.006, plane_sign=1)
    crop = heel.solve(plane_axis=1, plane_value_m=low + 0.006, plane_sign=1, crop_radius_m=3 * h)
    err = crop['segment_force_n'][1] / full['segment_force_n'][1] - 1
    gate('R1', 'crop at 3h puts <= 1% of the load on the truncation', crop['crop_load_fraction'] <= 0.01,
         f'truncation carries {crop["crop_load_fraction"]:.1%}; segment force error {err:+.1%}; '
         f'solved DOF {crop["solved_dof"]} of {full["solved_dof"]}')
    return {'load_on_truncation': crop['crop_load_fraction'], 'segment_force_error': err,
            'solved_dof': [crop['solved_dof'], full['solved_dof']]}


def main() -> int:
    began = time.perf_counter()
    report = {'bars_fixed_before_run': True, 'recorded_failures': sorted(RECORDED_FAILURES)}
    report['K1_hessian_rel'] = k1_hessian()
    print('\nK2  confined compression is the shipped foundation law', flush=True)
    report['K2_worst_rel'] = k2_confined()[1]
    report['K3'] = k3_unconfined()
    record = [r for r in json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_text())['records'] if r['body'] == HEEL][0]
    mesh = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
    low = float(mesh[0][:, 1].min())
    build = time.perf_counter()
    heel = stl.segment_layer(ROOT, HEEL)
    report['heel_build_seconds'] = time.perf_counter() - build
    heel_unconfined = stl.segment_layer(ROOT, HEEL, mapping='unconfined')
    report['heel'] = {'dof': heel.dof, 'nodes': int(len(heel.local)), 'tetrahedra': int(len(heel.tets)),
                      'base_nodes': int(heel.base.sum()), 'meta': heel.meta}
    k4_null(heel)
    first = k5_transmission(heel, low)
    report['K5'] = {k: first[k] for k in ('segment_force_n', 'contact_force_on_tissue_n', 'balance_force_relative',
                                          'balance_moment_residual_nm', 'iterations', 'wall_seconds')}
    k6_idempotent(heel, low, first)
    report['K7_energies'] = k7_energy()
    c3_energy_creating()
    report['K9'] = k9_cross()
    controls()
    k12_plant_options()
    report['M1_M2'] = m_heel(heel, heel_unconfined, mesh, record)
    report['R1'] = r1_crop(heel, low)
    report['M3'] = m3_convergence(record, mesh)
    report['census'], report['census_total'] = k11_unanchored()
    report['results'] = results
    report['failures'] = failures
    report['peak_memory_mb'] = peak_mb()
    report['wall_seconds'] = time.perf_counter() - began
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(report, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)) + '\n')
    print(f'\npeak memory {report["peak_memory_mb"]:.0f} MB, wall {report["wall_seconds"]:.0f} s; report {REPORT.relative_to(ROOT)}')
    print(f'{"ALL GATES PASS" if not failures else "FAILED: " + ", ".join(failures)}'
          f'{"" if not any(not results[k]["pass"] for k in RECORDED_FAILURES if k in results) else "  (recorded failures: " + ", ".join(k for k in RECORDED_FAILURES if k in results and not results[k]["pass"]) + ")"}')
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
