"""A deformable soft-tissue layer over a rigid segment: 3-D neo-Hookean, contact-local.

WHY THIS EXISTS.  docs/ACTUATION_STAGES.md asks for a `fully present participant`
body in which soft tissue deforms and mediates contact.  Before this module every
compliant element in the stack was a 1-D law on a rigid carrier:
`supine_contact.foundation` is an independent confined neo-Hookean COLUMN under each
quadrature point, and the segment contact meshes are an elastic foundation --
a spring per triangle.  Neither has a continuum: a column cannot bulge, cannot
share load with its neighbour, and cannot tell a heel-sized footprint from a
fingertip-sized one.  This module is the smallest thing that can.

WHAT IT IS.  Per scaffold segment, the soft tissue between the skin and the depth
at which this body's own depth map says bone or muscle begins, discretised as
constant-strain tetrahedra and carried by the segment:

* geometry -- the segment's own closed skin surface
  (`data/derived/segment-contact-meshes/skin-layer-map-v1`), voxelised by
  `material_domains.voxel_partition` (the repo's own mesher: six tets per cube, so
  no element can be born inverted).  A cell whose centre lies deeper than the
  segment's MEASURED soft-tissue depth `h` (the layer map's own per-segment value)
  is rigid core; the rest is the layer.  Every node of a core cell is a BASE node,
  rigidly carried by the segment.  Everything else is free.
* material -- `DeformableRegion`'s compressible neo-Hookean
  `W = mu/2 (I1-3) - mu ln J + lam/2 (ln J)^2`, the SAME law as the shipped
  foundation column, so the two can be compared exactly (see `LAYER_MODULUS_MAPPINGS`).
* contact -- a rigid frictionless half-space as an exact bound on every free node
  (L-BFGS-B box bound, as DeformableRegion's own indenter), not a penalty: there is
  no contact stiffness to choose.
* coupling back -- the reaction at the base nodes, summed, IS the force and moment
  the soft tissue transmits to the segment.  Internal forces of an isotropic
  hyperelastic energy sum to zero and carry no net moment, so at equilibrium this
  equals the contact resultant up to the solve residual; `balance_*` reports both.
* time -- quasi-static by default (the layer is massless and the segment carries
  the mass, as today), or backward-Euler incremental potential with lumped mass
  (`dt_s`), which can dissipate energy and cannot create it.

THE REDUCTION, AND WHY IT IS HONEST.  A whole body of this layer at 5 mm is
855,468 DOF over 1.51 M tetrahedra (scripts/verify_soft_tissue.py K11; docs/SOFT_BODY.md).  Two reductions,
each checkable:

1. A quasi-static layer under no load is exactly its rigidly carried rest state, so
   a segment out of contact costs NOTHING and returns an exact zero.  That is the
   null case of the battery, not an assumption.
2. `crop_radius_m`: only free nodes within that distance of the contact region
   move; the rest are held at the rigid rest state.  This is a Saint-Venant
   truncation, and its error is MEASURED against the uncropped solve, never assumed.
   MEASURED TO FAIL on the heel: at 3h the truncation carries 67.9% of the load,
   because the depth-defined core sits 30 mm above the plantar skin and the crop
   never reaches the thing that carries the load.  `crop_load_fraction` says so on
   every cropped call; read it before trusting a cropped force.

WHAT IT IS NOT.
* Not muscle.  Everything deeper than `h` is rigid core, which is also the shipped
  layer map's assumption (`k = E/h` over a rigid substrate).
* Not bone-anchored.  The bone meshes would be the right core, but the skin bundle
  records `bone_vertices_inside_skin` 0.444 and `segments_enclosing_their_bone` 0:
  skin and bone are not co-registered, so a bone core would be built on a known
  mismatch.  The core is the depth-defined one.
* Not continuous across segments.  Each segment's layer is independent, and every
  segment boundary is a seam, exactly as in the skin bundle it is built from.
* Not in the native integrator.  The plant does not call this; `plant_options`
  resolves it and hands the caller a layer.  See docs/SOFT_BODY.md.
* The surface is a VOXEL surface (staircase at the cell size), not the triangle
  skin.  At 5 mm that is a +-2.5 mm geometric error on where contact begins.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy import ndimage, sparse
from scipy.optimize import minimize
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu

from .mechanics_backend import DeformableRegion
from .material_domains import voxel_partition

BUNDLE = 'data/derived/segment-contact-meshes/skin-layer-map-v1'
BUNDLE_SCHEMA = 'ihm.segment-contact-meshes.v1'
DEFAULT_SPACING_M = 0.005

# The in-vivo source card (data/sources/in-vivo-soft-tissue-compression.json) states
# the heel modulus is pressure over THICKNESS strain -- an apparent LAYER modulus, so
# `k = E/h` with no confinement factor on top.  A 3-D material needs a Young's
# modulus, and the layer measurement does not fix one: it brackets it.
#
# * `confined` (default): choose E_young so the 3-D layer, loaded over a footprint
#   much wider than it is thick, has exactly the measured layer stiffness E_app/h.
#   Then (lam + 2 mu) = E_app, the 3-D layer reproduces the shipped foundation column
#   in the limit where that column is exact, and is SOFTER wherever tissue can bulge.
#   It is the softest reading of the measurement.
# * `unconfined`: E_young = E_app.  Correct only if the in-vivo footprint let the
#   tissue expand freely; 1/0.264 = 3.79x stiffer than `confined` under a wide load at
#   nu = 0.45.  The stiffest reading.
#
# The truth is between them and depends on the in-vivo footprint-to-thickness ratio,
# which the source card does not report.  Both are offered; neither is a fit.
LAYER_MODULUS_MAPPINGS = {
    'confined': 'E_young = E_app (1+nu)(1-2nu)/(1-nu): lam+2mu = E_app, so a wide load '
                'reproduces the measured layer stiffness E_app/h exactly (softest reading).',
    'unconfined': 'E_young = E_app: the apparent modulus read as a Young\'s modulus '
                  '(stiffest reading; 3.79x confined under a wide load at nu 0.45).',
}


def lame_from_layer_modulus(apparent_pa, poisson, mapping='confined'):
    """(mu, lam, young) for a layer modulus under a named, disclosed mapping."""
    e, nu = float(apparent_pa), float(poisson)
    if not np.isfinite(e) or e <= 0 or not 0 <= nu < 0.5:
        raise ValueError('Positive modulus and Poisson ratio in [0, 0.5) required')
    if mapping == 'confined':
        young = e * (1 + nu) * (1 - 2 * nu) / (1 - nu)
    elif mapping == 'unconfined':
        young = e
    else:
        raise ValueError('Unknown layer modulus mapping')
    return young / (2 * (1 + nu)), young * nu / ((1 + nu) * (1 - 2 * nu)), young


def confined_column_pressure(stretch, mu, lam):
    """The shipped 1-D law, `supine_contact.foundation`'s, as a compressive pressure."""
    s = np.asarray(stretch, float)
    return -(mu * (s - 1 / s) + lam * np.log(s) / s)


def element_hessians(region, positions, project=True):
    """12x12 Hessian of `DeformableRegion`'s own energy, per tetrahedron.

    With F = Ds Dm^-1, node a's shape gradient B_a, g_a = F^-T B_a and c = lam ln J - mu,
        H[a m, b n] = V [ mu (B_a . B_b) delta_mn + lam g_a[m] g_b[n] - c g_b[m] g_a[n] ].
    Verified against central differences of `energy_gradient` in
    scripts/verify_soft_tissue.py.  `project` clamps each element's negative
    eigenvalues to zero (Teran et al. 2005) so Newton always has a descent direction;
    it changes the step, never the energy being minimised.
    """
    f = region.deformation(positions)
    inverse_t = np.swapaxes(np.linalg.inv(f), 1, 2)
    c = region.lam * np.log(np.linalg.det(f)) - region.mu
    shape = np.concatenate([-region.inverse.sum(axis=1, keepdims=True), region.inverse], axis=1)
    g = np.einsum('tmj,taj->tam', inverse_t, shape)
    gram = np.einsum('taj,tbj->tab', shape, shape)
    h = (region.mu[:, None, None, None, None] * gram[:, :, None, :, None] * np.eye(3)[None, None, :, None, :]
         + region.lam[:, None, None, None, None] * g[:, :, :, None, None] * g[:, None, None, :, :]
         - c[:, None, None, None, None] * np.einsum('tbm,tan->tambn', g, g))
    h = h.reshape(-1, 12, 12) * region.volumes[:, None, None]
    h = 0.5 * (h + np.swapaxes(h, 1, 2))
    if project:
        w, v = np.linalg.eigh(h)
        h = np.einsum('tij,tj,tkj->tik', v, np.maximum(w, 0.), v)
    return h


def assemble(region, element):
    dof = (region.tets[:, :, None] * 3 + np.arange(3)).reshape(-1, 12)
    rows = np.repeat(dof[:, :, None], 12, axis=2)
    n = region.reference.size
    return sparse.coo_matrix((element.ravel(), (rows.ravel(), np.swapaxes(rows, 1, 2).ravel())),
                             shape=(n, n)).tocsc()


def load_obj(path):
    vertices, faces = [], []
    with open(path) as handle:
        for line in handle:
            if line.startswith('v '):
                vertices.append([float(t) for t in line.split()[1:4]])
            elif line.startswith('f '):
                faces.append([int(t.split('/')[0]) - 1 for t in line.split()[1:4]])
    return np.asarray(vertices, float), np.asarray(faces, int)


def layer_mesh(vertices, faces, thickness_m, spacing_m=DEFAULT_SPACING_M, name='segment'):
    """Tetrahedral soft-tissue shell of a closed surface: depth <= thickness is layer.

    Depth of a cell is its centre's Euclidean distance to the nearest OUTSIDE cell
    centre less half a cell -- the distance to the voxel boundary.  Returns node
    coordinates in the surface's own frame, tets, and the base (rigid) node mask.
    """
    h, sp = float(thickness_m), float(spacing_m)
    if not np.isfinite(h) or h <= 0 or not np.isfinite(sp) or sp <= 0:
        raise ValueError('Positive thickness and spacing required')
    vp = voxel_partition([(name, np.asarray(vertices, float), np.asarray(faces, int))], spacing_m=sp)
    cells = vp['cell_indices']
    shape = cells.max(axis=0) + 3
    occupied = np.zeros(tuple(shape), bool)
    occupied[cells[:, 0] + 1, cells[:, 1] + 1, cells[:, 2] + 1] = True
    depth = ndimage.distance_transform_edt(occupied, sampling=sp)[
        cells[:, 0] + 1, cells[:, 1] + 1, cells[:, 2] + 1] - sp / 2
    shell = depth <= h
    tets_by_cell = vp['tetrahedra'].reshape(-1, 6, 4)       # voxel_partition: six tets per cell, in cell order
    if len(tets_by_cell) != len(cells):
        raise ValueError('voxel_partition no longer emits six tets per cell in cell order')
    layer = tets_by_cell[shell].reshape(-1, 4)
    if not len(layer):
        raise ValueError('No layer cells')
    used = np.unique(layer)
    remap = np.full(len(vp['vertices_m']), -1)
    remap[used] = np.arange(len(used))
    core_nodes = np.intersect1d(np.unique(tets_by_cell[~shell].reshape(-1)), used)
    base = np.zeros(len(used), bool)
    base[remap[core_nodes]] = True
    tets = remap[layer]
    # a connected piece of layer that touches no core node is carried by nothing: it
    # cannot transmit load and would make the stiffness singular, so it is removed
    # and COUNTED, never silently pinned
    edges = np.concatenate([tets[:, [i, j]] for i in range(4) for j in range(i + 1, 4)])
    graph = sparse.coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(len(used),) * 2)
    count, label = connected_components(graph, directed=False)
    anchored = np.zeros(count, bool)
    anchored[np.unique(label[base])] = True
    keep_tet = anchored[label[tets[:, 0]]]
    islands = int(count - anchored.sum())
    island_tets = int((~keep_tet).sum())
    tets = tets[keep_tet]
    used2 = np.unique(tets)
    remap2 = np.full(len(used), -1)
    remap2[used2] = np.arange(len(used2))
    nodes = vp['vertices_m'][used][used2]
    base = base[used2]
    tets = remap2[tets]
    return {'unanchored_islands': islands, 'unanchored_island_tetrahedra': island_tets,
            'nodes_m': nodes, 'tetrahedra': tets, 'base': base,
            'spacing_m': sp, 'declared_thickness_m': h,
            'cells_through_thickness': h / sp,
            'occupied_cells': int(len(cells)), 'layer_cells': int(shell.sum()),
            'core_cells': int((~shell).sum()),
            'layer_volume_m3': float(len(tets)) * sp ** 3 / 6,
            'enclosed_volume_m3': float(len(cells)) * sp ** 3,
            'boundary_discretization_m': sp / 2}


class SoftTissueLayer(DeformableRegion):
    """One segment's soft tissue.  `solve` is a pure function of its arguments.

    Reference nodes are in the SEGMENT frame.  A pose `(rotation, translation)` maps
    them to world as `x_world = rotation @ x_local + translation`.  The energy is
    frame-indifferent, so the solve is done in world on the posed reference, and the
    half-space can then be an exact axis-aligned bound.
    """

    def __init__(self, mesh, *, mu_pa, lambda_pa, density_kg_m3, meta=None):
        base = np.asarray(mesh['base'], bool)
        if not base.any():
            raise ValueError('Unanchored layer: no node is deeper than the declared thickness, '
                             'so nothing carries it and it cannot transmit load to the segment')
        super().__init__(mesh['nodes_m'], np.asarray(mesh['tetrahedra']),
                         mu_pa=mu_pa, lambda_pa=lambda_pa, density_kg_m3=density_kg_m3)
        self.local = self.reference.copy()
        self.base = base
        self.mesh = {k: v for k, v in mesh.items() if k not in ('nodes_m', 'tetrahedra', 'base')}
        self.meta = dict(meta or {})
        self.dof = int((~base).sum()) * 3

    # -- pure solve -------------------------------------------------------------
    def rest_world(self, rotation, translation):
        r = np.asarray(rotation, float)
        t = np.asarray(translation, float)
        if r.shape != (3, 3) or t.shape != (3,) or not np.isfinite(r).all() or not np.isfinite(t).all():
            raise ValueError('Finite 3x3 rotation and 3-vector translation required')
        if not np.allclose(r @ r.T, np.eye(3), atol=1e-9) or np.linalg.det(r) < 0.999999:
            raise ValueError('Proper rotation required')
        return self.local @ r.T + t

    def solve(self, *, rotation=np.eye(3), translation=np.zeros(3), plane_axis=0, plane_value_m=None,
              plane_sign=1.0, dt_s=None, state=None, gravity_m_s2=(0., 0., 0.), held=None,
              crop_radius_m=None, method='newton', force_tolerance_n=None, maxiter=None,
              gtol=1e-12, return_positions=False):
        """Equilibrium (or one backward-Euler step) of the layer at a rigid segment pose.

        The support is the half-space `plane_sign * x[plane_axis] >= plane_sign * plane_value_m`.
        `state` is `None` (rest, at rest) or a previous result's `state`; it is only read.
        `held` is an optional (N, 3) bool of extra DOF carried rigidly with the segment
        (rollers, symmetry planes); their reactions count as the segment's.
        `method` is 'newton' (projected Newton on the analytic Hessian, the default) or
        'lbfgsb' (DeformableRegion's own minimiser, kept as the independent cross-check).
        Returns forces ON THE SEGMENT, the contact resultant ON THE TISSUE, energies,
        the balance residuals and the measured cost of this call.
        """
        began = time.perf_counter()
        rest = self.rest_world(rotation, translation)
        origin = np.asarray(translation, float)
        gravity = np.asarray(gravity_m_s2, float)
        if gravity.shape != (3,) or not np.isfinite(gravity).all():
            raise ValueError('Finite gravity required')
        if plane_axis not in (0, 1, 2) or plane_sign not in (1, -1, 1.0, -1.0):
            raise ValueError('Axis-aligned half-space required')
        dynamic = dt_s is not None
        if dynamic and (not np.isfinite(dt_s) or dt_s <= 0):
            raise ValueError('Positive finite dt required')
        if state is None:
            previous = rest.copy()
            velocity = np.zeros_like(rest)
        else:
            previous = np.asarray(state['positions_m'], float)
            velocity = np.asarray(state['velocities_m_s'], float)
            if previous.shape != rest.shape or velocity.shape != rest.shape:
                raise ValueError('State does not belong to this layer')
        carried = np.repeat(self.base[:, None], 3, axis=1)
        if held is not None:
            held = np.asarray(held, bool)
            if held.shape != rest.shape:
                raise ValueError('held must be (N, 3) bool')
            carried = carried | held
        sign = float(plane_sign)

        # --- which DOF move -------------------------------------------------------
        penetrating = np.zeros(len(rest), bool)
        if plane_value_m is not None:
            gap = sign * (rest[:, plane_axis] - plane_value_m)       # >0 clear, <0 in the support
            if (gap[carried[:, plane_axis]] < 0).any():
                raise ValueError('Bottomed out: the rigid core penetrates the support by %.4f m; '
                                 'the soft tissue cannot carry this pose'
                                 % float(-gap[carried[:, plane_axis]].min()))
            penetrating = gap < 0
        if not dynamic and not penetrating.any() and not gravity.any():
            # quasi-static, unloaded: the rigidly carried rest state is the exact minimum
            zero = np.zeros(3)
            result = self._receipt(rest, rest, zero, zero, zero, zero, zero, zero, 0.0, 0.0, 0, 0.0,
                                   time.perf_counter() - began, 'unloaded: exact rest state, no solve',
                                   0, 0, return_positions)
            result.update({'minimum_jacobian': 1.0, 'converged': True, 'active_tetrahedra': 0,
                           'crop_held_nodes': 0})
            result['state'] = {'positions_m': rest, 'velocities_m_s': np.zeros_like(rest)}
            return result
        cropped = np.zeros(len(rest), bool)
        if crop_radius_m is not None:
            if not penetrating.any():
                raise ValueError('crop_radius_m needs a contact region to crop around')
            from scipy.spatial import cKDTree
            distance, _ = cKDTree(rest[penetrating]).query(rest, k=1)
            cropped = ~self.base & (distance > float(crop_radius_m))
        crop_dof = np.repeat(cropped[:, None], 3, axis=1) & ~carried
        free_dof = ~carried & ~crop_dof
        active_tets = free_dof.any(axis=1)[self.tets].any(axis=1)
        tets = self.tets[active_tets]
        nodes = np.unique(tets)
        index = np.full(len(rest), -1)
        index[nodes] = np.arange(len(nodes))
        sub = DeformableRegion(rest[nodes], index[tets], mu_pa=self.mu[active_tets],
                               lambda_pa=self.lam[active_tets], density_kg_m3=1.0)
        free = free_dof[nodes]
        carry = carried[nodes]
        crop = crop_dof[nodes]
        fixed = ~free
        mass = self.nodal_mass[nodes]           # the FULL layer's lumped mass, restricted

        # --- bounds -----------------------------------------------------------------
        lo = np.full((len(nodes), 3), -np.inf)
        hi = np.full((len(nodes), 3), np.inf)
        lo[fixed] = 0.
        hi[fixed] = 0.
        bounded = free[:, plane_axis]
        if plane_value_m is not None:
            limit = plane_value_m - rest[nodes, plane_axis]
            if sign > 0:
                lo[bounded, plane_axis] = np.maximum(lo[bounded, plane_axis], limit[bounded])
            else:
                hi[bounded, plane_axis] = np.minimum(hi[bounded, plane_axis], limit[bounded])

        # --- objective: E(y) + sum m/(2dt^2)|y - yhat|^2 - sum m g.y  -------------------
        ref = rest[nodes]
        external = mass[:, None] * gravity
        predicted = previous[nodes] + dt_s * velocity[nodes] - ref if dynamic else None
        inertia = mass[:, None] / dt_s ** 2 if dynamic else None
        scale = float(self.mesh.get('spacing_m') or np.min(np.ptp(ref, axis=0)))
        energy_scale = float(np.mean(sub.mu)) * scale ** 3

        def objective(z):
            u = (z * scale).reshape(-1, 3)
            energy, gradient = sub.energy_gradient(ref + u, trial=True)
            energy -= float(np.sum(external * u))
            gradient = gradient - external
            if dynamic:
                delta = u - predicted
                energy += 0.5 * float(np.sum(inertia * delta ** 2))
                gradient = gradient + inertia * delta
            gradient = np.where(fixed, 0., gradient)
            return energy / energy_scale, gradient.ravel() * scale / energy_scale

        # deterministic cold start: the rigid rest state (or the inertial prediction),
        # projected onto the bounds -- never the last call's answer
        unclipped = np.where(fixed, 0., predicted) if dynamic else np.zeros_like(ref)
        if dynamic and np.linalg.det(sub.deformation(ref + unclipped)).min() <= 0:
            # the inertial prediction itself inverts an element: start from the last accepted
            # positions instead, which are valid by construction (still deterministic)
            unclipped = np.where(fixed, 0., previous[nodes] - ref)
        start = np.clip(unclipped, lo, hi)
        if force_tolerance_n is None:
            force_tolerance_n = 1e-8 * float(np.mean(sub.mu)) * scale ** 2
        if method == 'lbfgsb':
            result = minimize(objective, start.ravel() / scale, jac=True, method='L-BFGS-B',
                              bounds=list(zip(lo.ravel() / scale, hi.ravel() / scale)),
                              options={'ftol': 1e-16, 'gtol': gtol, 'maxiter': int(maxiter or 20000),
                                       'maxls': 50, 'maxcor': 30})
            u = result.x.reshape(-1, 3) * scale
            iterations, success, message = int(result.nit), bool(result.success), str(result.message)
        elif method == 'newton':
            # A start that violates the bound by more than about a cell inverts the elements
            # it is clipped into.  So the bound is brought in by load steps, each warm-started
            # from the last (prescribed_deformation.py's precedent): the whole way at once
            # first, halving the increment only when the clipped start is actually inverted.
            # The path is fixed by the arguments, so the result is still a function of them.
            relax_lo = np.minimum(unclipped - lo, 0.)
            relax_hi = np.maximum(unclipped - hi, 0.)
            u, iterations, success, message = unclipped.copy(), 0, True, ''
            remaining, increment, steps, halvings = 1.0, 1.0, 0, 0
            while remaining > 0:
                target = max(remaining - increment, 0.)
                lo_k, hi_k = lo + target * relax_lo, hi + target * relax_hi
                last = target == 0.
                trial, used, success, message = self._projected_newton(
                    sub, ref, np.clip(u, lo_k, hi_k), lo_k, hi_k, fixed, external, predicted, inertia,
                    float(force_tolerance_n) * (1. if last else 1e3), int(maxiter or 200))
                iterations += used
                if success is None:                 # the clipped start inverted: smaller increment
                    increment /= 2
                    halvings += 1
                    if increment < 1.0 / 1024:
                        success, message = False, 'load step increment below 1/1024 with an inverted start'
                        break
                    continue
                if not success:
                    # report the LAST ITERATE, flagged non-converged -- never the start, which
                    # is a state no solver produced (the first K5 run read 0 N because of this)
                    u = trial
                    break
                u, remaining, steps = trial, target, steps + 1
            message = '%d load step(s), %d halving(s); %s' % (steps, halvings, message)
        else:
            raise ValueError('Unknown solve method')
        y = ref + u
        jac = np.linalg.det(sub.deformation(y))
        if jac.min() <= 0.2:
            raise RuntimeError('Left the constitutive domain: minimum J %.3f <= 0.2' % jac.min())
        elastic, g = sub.energy_gradient(y)
        g = g - external
        if dynamic:
            g = g + inertia * (u - predicted)

        # --- who pushes whom ---------------------------------------------------------
        # g = dPhi/dy.  At a prescribed DOF the constraint applies +g to the tissue, so the
        # tissue applies -g to whatever holds it.  At an active bound the KKT sign is
        # sign * g_axis >= 0: the plane pushes the tissue out of the support.
        on_plane = np.zeros(len(nodes), bool)
        if plane_value_m is not None:
            bound = lo[:, plane_axis] if sign > 0 else hi[:, plane_axis]
            on_plane = bounded & np.isclose(u[:, plane_axis], bound, rtol=0, atol=1e-12)
        pushed = on_plane & (sign * g[:, plane_axis] > 0)
        contact_nodal = np.zeros_like(g)
        contact_nodal[pushed, plane_axis] = g[pushed, plane_axis]
        contact_force = contact_nodal.sum(axis=0)
        contact_moment = np.cross(y - origin, contact_nodal).sum(axis=0)
        residual = np.where(fixed, 0., g - contact_nodal)
        on_segment = np.where(carry, -g, 0.)
        segment_force = on_segment.sum(axis=0)
        segment_moment = np.cross(y - origin, on_segment).sum(axis=0)
        # cropped DOF are held by the truncation, NOT by the segment; their reaction is
        # the crop's error and is reported as such, never added to what the segment feels
        on_crop = np.where(crop, -g, 0.)
        crop_force = on_crop.sum(axis=0)
        crop_moment = np.cross(y - origin, on_crop).sum(axis=0)
        # the layer's own gravity and inertia (zero when quasi-static and weightless)
        body_nodal = -external + (inertia * (u - predicted) if dynamic else 0.)
        # Internal forces of a frame-indifferent energy sum to zero and carry no moment, so
        #   segment + crop - contact + body = sum of the free residual   (force and moment).
        balance_force = segment_force + crop_force - contact_force + body_nodal.sum(axis=0)
        balance_moment = (segment_moment + crop_moment - contact_moment
                          + np.cross(y - origin, body_nodal).sum(axis=0))

        full = rest.copy()
        full[nodes] = y
        if dynamic:
            new_velocity = (rest - previous) / dt_s          # outside the active set: rigid carry
            new_velocity[nodes] = (y - previous[nodes]) / dt_s
            kinetic = 0.5 * float(np.sum(self.nodal_mass[:, None] * new_velocity ** 2))
        else:
            new_velocity = np.zeros_like(rest)
            kinetic = 0.0
        receipt = self._receipt(full, rest, segment_force, segment_moment, contact_force,
                                contact_moment, balance_force, balance_moment, elastic, kinetic,
                                iterations, float(np.max(np.abs(residual))),
                                time.perf_counter() - began, message,
                                int(pushed.sum()), int(free.any(axis=1).sum()), return_positions)
        receipt.update({'minimum_jacobian': float(jac.min()),
                        'crop_reaction_n': crop_force.tolist(),
                        'crop_reaction_moment_nm': crop_moment.tolist(),
                        # the share of the load the TRUNCATION carried instead of the segment:
                        # a crop that does not enclose the load path to the core puts it all here
                        'crop_load_fraction': float(np.linalg.norm(crop_force) / np.linalg.norm(contact_force))
                        if np.linalg.norm(contact_force) > 0 else 0.0,
                        'layer_body_force_n': body_nodal.sum(axis=0).tolist(),
                        'crop_held_nodes': int(crop.any(axis=1).sum()),
                        'active_tetrahedra': int(active_tets.sum()),
                        'solved_dof': int(free.sum()),
                        'converged': bool(success), 'method': method,
                        'force_tolerance_n': float(force_tolerance_n)})
        receipt['state'] = {'positions_m': full, 'velocities_m_s': new_velocity}
        return receipt

    @staticmethod
    def _projected_newton(sub, ref, start, lo, hi, fixed, external, predicted, inertia, tolerance, maxiter):
        """Minimise Phi = E + inertia - external over the box [lo, hi], Newton on the free set.

        The active set is every prescribed DOF plus every DOF sitting on its bound with the
        gradient pushing it further out.  Each step solves the PSD-projected Newton system on
        the rest and backtracks along the PROJECTED path until Phi decreases (Armijo), so
        Phi is monotone and no accepted iterate has an inverted element.
        """
        dynamic = inertia is not None

        def phi(u, need_gradient=True):
            e, g = sub.energy_gradient(ref + u)                 # raises on inversion
            e -= float(np.sum(external * u))
            g = g - external
            if dynamic:
                d = u - predicted
                e += 0.5 * float(np.sum(inertia * d ** 2))
                g = g + inertia * d
            return e, g

        u = start.copy()
        try:
            value, g = phi(u)
        except ValueError:
            return u, 0, None, 'start inverted'
        span = np.maximum(np.abs(ref).max(), 1.0)
        for iteration in range(1, maxiter + 1):
            at_lo = u <= lo + 1e-14 * span
            at_hi = u >= hi - 1e-14 * span
            active = fixed | (at_lo & (g > 0)) | (at_hi & (g < 0))
            projected = np.where(active, 0., g)
            if np.max(np.abs(projected)) <= tolerance:
                return u, iteration - 1, True, 'projected gradient below %.3g N' % tolerance
            free = ~active.ravel()
            rhs = g.ravel()[free]
            step = None
            for project in (False, True):
                k = assemble(sub, element_hessians(sub, ref + u, project=project))
                if dynamic:
                    k = k + sparse.diags(np.broadcast_to(inertia, ref.shape).ravel())
                kf = k[free][:, free]
                # a free DOF with no stiffness gets a tiny diagonal so the factor exists;
                # it cannot change the minimiser, only whether the system is solvable
                diagonal = kf.diagonal()
                kf = kf + sparse.diags(np.where(diagonal > 0, 0., max(float(diagonal.max()), 1.) * 1e-12))
                try:
                    candidate = -splu(kf.tocsc(), permc_spec='COLAMD').solve(rhs)
                except RuntimeError:
                    continue
                if np.isfinite(candidate).all() and float(candidate @ rhs) < 0:
                    step = np.zeros(u.size)
                    step[free] = candidate
                    step = step.reshape(-1, 3)
                    break
            if step is None:
                return u, iteration, False, 'no descent direction at projected gradient %.3g N' % np.max(np.abs(projected))
            alpha = 1.0
            while True:
                trial = np.clip(u + alpha * step, lo, hi)
                try:
                    trial_value, trial_g = phi(trial)
                except ValueError:
                    trial_value = np.inf
                if trial_value <= value + 1e-4 * float(np.sum(g * (trial - u))) and trial_value < value:
                    break
                alpha *= 0.5
                if alpha < 1e-12:
                    # Phi can no longer be resolved in floating point.  Converged only if the
                    # projected gradient is already within 100x the tolerance; the achieved value
                    # is reported either way (receipt `free_residual_n`), never hidden.
                    level = float(np.max(np.abs(projected)))
                    return (u, iteration, level <= 100 * tolerance,
                            'stagnated at the floating-point floor, projected gradient %.3g N' % level)
            u, value, g = trial, trial_value, trial_g
        return u, maxiter, False, 'Newton iteration limit reached'

    def _receipt(self, positions, rest, segment_force, segment_moment, contact_force,
                 contact_moment, balance_force, balance_moment, elastic, kinetic, iterations,
                 residual, wall, message, contact_nodes, free_nodes, return_positions):
        magnitude = np.linalg.norm(positions - rest, axis=1)
        contact = float(np.linalg.norm(contact_force))
        out = {'segment_force_n': np.asarray(segment_force, float).tolist(),
               'segment_moment_nm': np.asarray(segment_moment, float).tolist(),
               'contact_force_on_tissue_n': np.asarray(contact_force, float).tolist(),
               'contact_moment_about_segment_nm': np.asarray(contact_moment, float).tolist(),
               'balance_force_residual_n': np.asarray(balance_force, float).tolist(),
               'balance_moment_residual_nm': np.asarray(balance_moment, float).tolist(),
               # relative to the LARGEST force in the balance, not to the contact alone: a solve
               # that never reached contact has zero contact force, and dividing by it (or
               # defaulting to 0) certified a starved, unconverged solve as balanced (C2, run 3)
               'balance_force_relative': _relative(balance_force, contact_force, segment_force),
               'elastic_energy_j': float(elastic), 'kinetic_energy_j': float(kinetic),
               'maximum_displacement_m': float(magnitude.max(initial=0.)),
               'contact_nodes': int(contact_nodes), 'free_nodes_solved': int(free_nodes),
               'layer_dof': self.dof, 'solved_dof': int(free_nodes) * 3,
               'iterations': int(iterations), 'free_residual_n': float(residual),
               'wall_seconds': float(wall), 'solver_message': message}
        if return_positions:
            out['positions_m'] = positions
        return out


def _relative(residual, *scales):
    size = max(float(np.linalg.norm(v)) for v in scales)
    norm = float(np.linalg.norm(residual))
    if size > 0:
        return norm / size
    return 0.0 if norm == 0 else float('inf')


def _verified_record(root, bundle, body):
    root = Path(root)
    path = (root / bundle).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError('Owned bundle required')
    manifest = json.loads((path / 'manifest.json').read_bytes())
    if manifest.get('schema') != BUNDLE_SCHEMA or manifest.get('layer') != 'skin':
        raise ValueError('Not a skin segment-contact bundle')
    if 'layer_map' not in manifest:
        raise ValueError('Bundle carries no measured per-segment layer depth')
    records = [r for r in manifest['records'] if r['body'] == body]
    if len(records) != 1:
        raise ValueError('Segment has no admitted skin record: ' + str(body))
    record = records[0]
    mesh_path = path / 'meshes' / record['mesh_file']
    if hashlib.sha256(mesh_path.read_bytes()).hexdigest() != record['written_sha256']:
        raise ValueError('Skin mesh does not match the bundle record: ' + record['mesh_file'])
    return manifest, record, mesh_path


def segment_layer(root, body, *, bundle=BUNDLE, spacing_m=DEFAULT_SPACING_M, mapping='confined'):
    """Build one segment's layer from the measured layer map, and say what it rests on."""
    from .tissue_materials import DENSITY
    manifest, record, mesh_path = _verified_record(root, bundle, body)
    layer = record['layer']
    poisson = float(manifest['skin_material']['poissons_ratio'])
    mu, lam, young = lame_from_layer_modulus(layer['apparent_modulus_pa'], poisson, mapping)
    rho = DENSITY['adipose'][0]
    vertices, faces = load_obj(mesh_path)
    mesh = layer_mesh(vertices, faces, layer['thickness_m'], spacing_m, name=body)
    meta = {'body': body, 'bundle': bundle, 'mesh_sha256': record['written_sha256'],
            'thickness_m': layer['thickness_m'], 'depth_points': layer['depth_points'],
            'apparent_modulus_pa': layer['apparent_modulus_pa'],
            'modulus_basis': layer['modulus_basis'], 'poisson_ratio': poisson,
            'poisson_basis': 'the bundle\'s declared skin_material poissons_ratio',
            'mapping': mapping, 'mapping_basis': LAYER_MODULUS_MAPPINGS[mapping],
            'young_modulus_pa': young, 'mu_pa': mu, 'lambda_pa': lam,
            'density_kg_m3': rho, 'density_basis': 'tissue_materials.DENSITY adipose (ICRU-44)',
            'layer_map_stiffness_pa_per_m': layer['stiffness_pa_per_m']}
    return SoftTissueLayer(mesh, mu_pa=mu, lambda_pa=lam, density_kg_m3=rho, meta=meta)


def build_selected_layers(root, selection, bodies=None):
    """Build the layers a resolved `plant_options` selection names, and nothing else.

    `selection` is the second return of `resolve_fidelity` (or a body's
    `mechanical_fidelity`).  The caller cannot pass a bundle, a modulus or a mesh: every
    one of them comes from the server-owned identity the selection carries.
    """
    chosen = (selection or {}).get('soft_tissue')
    if not chosen:
        raise ValueError('The selection does not ask for a soft tissue layer')
    from .plant_options import SOFT_TISSUE_LAYERS
    spec = SOFT_TISSUE_LAYERS.get(chosen.get('id'))
    if spec is None or chosen.get('bundle') != spec['bundle'] or chosen.get('mapping') != spec['mapping']:
        raise ValueError('Soft tissue selection does not match a server-owned identity')
    wanted = chosen['segments'] if bodies is None else list(bodies)
    unknown = set(wanted) - set(chosen['segments'])
    if unknown:
        raise ValueError('Not an anchored segment of this selection: ' + ', '.join(sorted(unknown)))
    return {body: segment_layer(root, body, bundle=spec['bundle'], spacing_m=spec['spacing_m'],
                                mapping=spec['mapping']) for body in wanted}
