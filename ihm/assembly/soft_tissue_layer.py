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
   MEASURED TO FAIL on calcn_l: at 3h the truncation carries 67.9% of the load.  The
   fixture presses the seam wedge on the calcn/toes cut, whose nearest core is at the
   midfoot, more than 30 mm away, so the crop never reaches the thing that carries the load.
   (This was first written as "the core sits 30 mm above the plantar skin", which measured
   between two places 55 mm apart; docs/SOFT_BODY.md.)  `crop_load_fraction` says so on
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
* By default (`surface='voxel'`) the surface is a VOXEL surface (staircase at the cell
  size), not the triangle skin.  At 5 mm that is a +-2.5 mm geometric error on where
  contact begins.  `surface='fitted'` cuts the lattice by the skin and by the declared
  depth instead (`fitted_layer_mesh`), and `depth='local'` reads the depth map under each
  point of the skin instead of the segment's median.  The fitted heel force is monotone in
  the spacing but not converged: first order, ~8-15% from its extrapolated limit at the
  finest spacing affordable (docs/SOFT_BODY.md).
* `method='fast'` is the same minimisation to the same tolerance with a factorisation
  reused across iterations, steps and poses: 0.17-0.44 s per cold heel solve and ~0.14 s
  per warm-started step at 7,296 DOF.  That is not real-time; the plant step is 10 ms.
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
# Per-vertex depth from this body's depth map, in each segment's frame
# (scripts/build_soft_tissue_local_depth.py, which gates the map on two known answers).
LOCAL_DEPTH = 'data/derived/soft-tissue-local-depth-v1'
LOCAL_DEPTH_SCHEMA = 'ihm.soft-tissue-local-depth.v1'
SURFACES = ('voxel', 'fitted')
# When a list, `_projected_newton_cg` appends one record per iteration (diagnostics only).
TRACE = None
DEPTH_RULES = ('segment_median', 'local')

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


# ------------------------------------------------------------------ body-fitted boundary
# The voxel surface is a staircase: on calcn_l its lowest point sits 0.85 / 2.85 / 2.85 mm
# above the skin's at 5 / 4 / 6 mm cells, and the contact force followed where the staircase
# landed rather than the spacing (M3: 0.457 / 0.657 / 1.113 N, non-monotone).  The fitted mesh
# keeps the voxel CONNECTIVITY (six tets per cube, so nothing is born inverted) and moves the
# NODES: the lattice is cut by the skin and by the declared depth, isosurface-stuffing style
# (`fitted_layer_mesh`).  Two earlier versions MOVED the voxel surface nodes instead, and both
# failed on the same object, recorded in docs/SOFT_BODY.md: a tet whose four nodes all land on one
# smooth surface is flat, no interior motion can fix it, and its volume ratio forms a continuum
# (calcn_l, 5 mm: -0.17 to 0.2 at the median) with no gap to cut at.  (1) A penalty on distance,
# raised in decades over a mesh-quality energy, stalled at 0.14 mm / 1.07 mm after 353 s.
# (2) Projecting exactly and relaxing the interior reached the surfaces exactly but kept
# near-flat slivers (volume ratio < 5e-7), which would stiffen the contact, and could not carry
# the local rule's thin layers (43-49% of the load).

def _closest_on_triangles(p, a, b, c):
    """Closest point on triangle (a, b, c) to p, row-wise (Ericson, Real-Time Collision Detection 5.1.5)."""
    ab, ac = b - a, c - a
    ap, bp, cp = p - a, p - b, p - c
    d1, d2 = np.einsum('ij,ij->i', ab, ap), np.einsum('ij,ij->i', ac, ap)
    d3, d4 = np.einsum('ij,ij->i', ab, bp), np.einsum('ij,ij->i', ac, bp)
    d5, d6 = np.einsum('ij,ij->i', ab, cp), np.einsum('ij,ij->i', ac, cp)
    va, vb, vc = d3 * d6 - d5 * d4, d5 * d2 - d1 * d6, d1 * d4 - d3 * d2
    with np.errstate(divide='ignore', invalid='ignore'):
        denom = va + vb + vc
        out = a + (vb / denom)[:, None] * ab + (vc / denom)[:, None] * ac          # face interior
        w = (d4 - d3) / ((d4 - d3) + (d5 - d6))
        out = np.where(((va <= 0) & (d4 - d3 >= 0) & (d5 - d6 >= 0))[:, None], b + w[:, None] * (c - b), out)
        w = d2 / (d2 - d6)
        out = np.where(((vb <= 0) & (d2 >= 0) & (d6 <= 0))[:, None], a + w[:, None] * ac, out)
        out = np.where(((d6 >= 0) & (d5 <= d6))[:, None], c, out)
        v = d1 / (d1 - d3)
        out = np.where(((vc <= 0) & (d1 >= 0) & (d3 <= 0))[:, None], a + v[:, None] * ab, out)
    out = np.where(((d3 >= 0) & (d4 <= d3))[:, None], b, out)
    out = np.where(((d1 <= 0) & (d2 <= 0))[:, None], a, out)
    return out


def closest_points(points, vertices, faces, chunk=4096):
    """EXACT closest point on a triangle set for every query point: (distance, point, face).

    Candidates are pruned by a bound that cannot drop the answer: the nearest vertex gives an
    upper bound D on the distance, and a face whose bounding sphere lies farther than D cannot
    hold the nearest point.  No face-size assumption, so a long facet tail costs time, never
    correctness.  Ties are broken by the lowest face index, so the result is a function of
    its arguments.
    """
    from scipy.spatial import cKDTree
    p = np.asarray(points, float)
    tri = np.asarray(vertices, float)[np.asarray(faces)]
    centre = tri.mean(axis=1)
    radius = np.linalg.norm(tri - centre[:, None], axis=2).max(axis=1)
    used = np.unique(np.asarray(faces))
    upper, _ = cKDTree(np.asarray(vertices, float)[used]).query(p)
    # faces far larger than the median (joint caps are centroid fans tens of mm across) would
    # widen every query; they are few, so they are tested against every point instead
    large = radius > 4 * np.median(radius)
    small = np.flatnonzero(~large)
    tree = cKDTree(centre[small])
    reach = float(radius[small].max(initial=0.))
    big = np.flatnonzero(large)
    distance = np.empty(len(p))
    nearest = np.empty_like(p)
    face = np.empty(len(p), int)
    for s in range(0, len(p), chunk):
        q, ub = p[s:s + chunk], upper[s:s + chunk]
        lists = tree.query_ball_point(q, ub + reach)
        counts = np.fromiter((len(c) for c in lists), int, len(lists)) + len(big)
        pi = np.repeat(np.arange(len(q)), counts)
        ti = np.fromiter((t for c in lists for t in (*small[c], *big)), int, int(counts.sum()))
        keep = np.linalg.norm(q[pi] - centre[ti], axis=1) - radius[ti] <= ub[pi] * (1 + 1e-12) + 1e-15
        pi, ti = pi[keep], ti[keep]
        cp = _closest_on_triangles(q[pi], tri[ti, 0], tri[ti, 1], tri[ti, 2])
        d2 = np.einsum('ij,ij->i', cp - q[pi], cp - q[pi])
        order = np.lexsort((ti, d2, pi))
        first = order[np.r_[True, pi[order][1:] != pi[order][:-1]]]
        if len(first) != len(q):
            raise RuntimeError('closest_points lost a query point; the pruning bound is wrong')
        distance[s:s + chunk] = np.sqrt(d2[first])
        nearest[s:s + chunk] = cp[first]
        face[s:s + chunk] = ti[first]
    return distance, nearest, face


def _barycentric(points, a, b, c):
    v0, v1, v2 = b - a, c - a, points - a
    d00, d01, d11 = (v0 * v0).sum(1), (v0 * v1).sum(1), (v1 * v1).sum(1)
    d20, d21 = (v2 * v0).sum(1), (v2 * v1).sum(1)
    den = d00 * d11 - d01 * d01
    v = (d11 * d20 - d01 * d21) / den
    w = (d00 * d21 - d01 * d20) / den
    return np.stack([1 - v - w, v, w], axis=1)


class DepthSurface:
    """Where the layer ends: `depth(x)` is the distance from x to the skin, `h(x)` the declared
    layer thickness at x's nearest skin point.  `h` is one number (`segment_median`, the bundle's
    per-segment median, measured to the CAPPED mesh exactly as the voxel rule does) or the depth
    map's own value interpolated on the nearest skin triangle (`local`, measured to the SKIN
    faces only: a joint cap is a cut through the body, not skin, and the depth map has no point
    on it)."""

    def __init__(self, vertices, faces, thickness_m=None, vertex_depth_m=None, skin_faces=None):
        self.vertices = np.asarray(vertices, float)
        faces = np.asarray(faces)
        if (thickness_m is None) == (vertex_depth_m is None):
            raise ValueError('Exactly one of a uniform thickness and a per-vertex depth')
        if vertex_depth_m is None:
            self.faces, self.uniform, self.vertex_depth = faces, float(thickness_m), None
            self.rule = 'segment_median'
        else:
            mask = np.asarray(skin_faces, bool)
            self.faces = faces[mask]
            self.vertex_depth = np.asarray(vertex_depth_m, float)
            if not np.isfinite(self.vertex_depth[np.unique(self.faces)]).all():
                raise ValueError('A skin face has a vertex with no measured depth')
            self.uniform, self.rule = None, 'local'

    def query(self, points):
        distance, nearest, face = closest_points(points, self.vertices, self.faces)
        if self.uniform is not None:
            return distance, nearest, np.full(len(distance), self.uniform)
        f = self.faces[face]
        weights = _barycentric(nearest, *(self.vertices[f[:, k]] for k in range(3)))
        return distance, nearest, np.einsum('ij,ij->i', np.clip(weights, 0., 1.), self.vertex_depth[f]) \
            / np.clip(weights, 0., 1.).sum(axis=1)


def _boundary_nodes(tetrahedra, n_vertices):
    faces = np.sort(np.concatenate([tetrahedra[:, [0, 1, 2]], tetrahedra[:, [0, 1, 3]],
                                    tetrahedra[:, [0, 2, 3]], tetrahedra[:, [1, 2, 3]]]), axis=1)
    unique, counts = np.unique(faces, axis=0, return_counts=True)
    mask = np.zeros(n_vertices, bool)
    mask[unique[counts == 1].ravel()] = True
    return mask


def _level_set(points, depth, sweeps=50):
    """Move points onto {depth(x) = h}: repeat x <- cp + h (x - cp)/|x - cp| until it stops moving.
    One sweep is exact for a single smooth sheet; where another part of the skin is nearer
    (near a cap, across a thin section) the nearest point changes and it takes more."""
    x = np.asarray(points, float).copy()
    for _ in range(sweeps):
        d, cp, h = depth.query(x)
        y = cp + (x - cp) / d[:, None] * h[:, None]
        moved = float(np.abs(y - x).max(initial=0.))
        x = y
        if moved <= 1e-12:
            break
    return x


# Labelle & Shewchuk, "Isosurface stuffing: fast tetrahedral meshes with good dihedral angles",
# ACM TOG 26(3) 2007: a lattice node whose cut point on an edge lies within ALPHA of the edge's
# length is moved onto the surface before cutting, so no cut makes a sliver.  0.24969 is their
# alpha_long for the BCC lattice.  This is a Kuhn lattice, for which they prove nothing: the value is
# BORROWED, not derived, and the mesh quality it gives here is measured and reported on every build
# (`fit['minimum_quality']`), which is what would show it wrong.
STUFFING_ALPHA = 0.24969


def _kuhn_lattice(cells, origin, spacing):
    """Nodes and six-tet Kuhn cubes over integer cells (the same pattern as voxel_partition)."""
    import itertools
    corners = np.array(list(itertools.product((0, 1), repeat=3)))
    raw = (cells[:, None, :] + corners).reshape(-1, 3)
    unique, index = np.unique(raw, axis=0, return_inverse=True)
    vertex_map = index.reshape(-1, 8)
    pattern = []
    for order in itertools.permutations(range(3)):
        v = np.zeros(3, int)
        tet = [0]
        for axis in order:
            v = v.copy()
            v[axis] += 1
            tet.append(int(np.flatnonzero(np.all(corners == v, axis=1))[0]))
        pattern.append(tet)
    tets = vertex_map[:, np.array(pattern)].reshape(-1, 4)
    nodes = origin + unique * spacing
    det = np.linalg.det(np.swapaxes(nodes[tets[:, 1:]] - nodes[tets[:, 0, None]], 1, 2))
    tets[det < 0] = tets[det < 0][:, [0, 2, 1, 3]]
    return nodes, tets


def _edges(tets):
    e = np.sort(np.concatenate([tets[:, [i, j]] for i in range(4) for j in range(i + 1, 4)]), axis=1)
    return np.unique(e, axis=0)


def _warp(x, phi, tets, movable, alpha):
    """Stuffing's warp: a node within alpha of a cut point on one of its edges moves to the nearest such cut
    point and becomes a surface node (phi = 0).  Returns the moved positions and which nodes moved."""
    e = _edges(tets)
    a, b = e[:, 0], e[:, 1]
    crossing = (phi[a] * phi[b] < 0)
    a, b = a[crossing], b[crossing]
    t = phi[a] / (phi[a] - phi[b])
    cut = x[a] + t[:, None] * (x[b] - x[a])
    node = np.concatenate([a, b])
    frac = np.concatenate([t, 1 - t])
    point = np.concatenate([cut, cut])
    ok = (frac < alpha) & movable[node]
    node, frac, point = node[ok], frac[ok], point[ok]
    order = np.lexsort((frac, node))
    first = order[np.r_[True, node[order][1:] != node[order][:-1]]] if len(order) else order
    moved = np.zeros(len(x), bool)
    y = x.copy()
    y[node[first]] = point[first]
    moved[node[first]] = True
    return y, moved


def _cut(tets, n_vertices, phi):
    """Keep {phi >= 0} of each tet (marching tetrahedra with consistent prism/pyramid splits).

    Cut points are keyed by their edge, so neighbours share them; every quad is split along the
    diagonal through its lowest-indexed vertex, which makes the split of a shared face the same
    from both sides (Dompierre et al. 1999).  A tet with all four vertices at phi = 0 is discarded
    (stuffing's rule).  Returns the new tets, and for each NEW vertex its edge (a, b) and t.
    """
    sign = np.sign(phi)
    s = sign[tets]
    pos, neg = (s > 0).sum(1), (s < 0).sum(1)
    whole = tets[(neg == 0) & (pos > 0)]
    cutting = np.flatnonzero((neg > 0) & (pos > 0))
    key = {}
    new_edges, new_t = [], []

    def cp(i, j):
        k = (i, j) if i < j else (j, i)
        if k not in key:
            key[k] = n_vertices + len(new_edges)
            new_edges.append(k)
            new_t.append(phi[k[0]] / (phi[k[0]] - phi[k[1]]))
        return key[k]

    def pyramid(q, apex):
        m = int(np.argmin(q))
        if m in (0, 2):
            return [(q[0], q[1], q[2], apex), (q[0], q[2], q[3], apex)]
        return [(q[1], q[2], q[3], apex), (q[1], q[3], q[0], apex)]

    def prism(bottom, top):
        allv = list(bottom) + list(top)
        m = int(np.argmin(allv))
        if m >= 3:
            bottom, top = top, bottom
            m -= 3
        r = [m, (m + 1) % 3, (m + 2) % 3]
        A, B, C = (bottom[i] for i in r)
        D, E, F = (top[i] for i in r)
        if min(B, F) < min(C, E):
            return [(A, B, C, F), (A, B, F, E), (A, E, F, D)]
        return [(A, B, C, E), (A, E, C, F), (A, E, F, D)]

    out = []
    for k in cutting:
        tet = tets[k]
        sg = sign[tet]
        P = [int(v) for v, g in zip(tet, sg) if g > 0]
        Z = [int(v) for v, g in zip(tet, sg) if g == 0]
        N = [int(v) for v, g in zip(tet, sg) if g < 0]
        if len(P) == 1:
            out.append((P[0], *Z, *(cp(P[0], n) for n in N)))
        elif len(P) == 2 and len(N) == 2:
            out += prism((P[0], cp(P[0], N[0]), cp(P[0], N[1])), (P[1], cp(P[1], N[0]), cp(P[1], N[1])))
        elif len(P) == 2:
            out += pyramid((P[0], P[1], cp(P[1], N[0]), cp(P[0], N[0])), Z[0])
        else:
            out += prism(tuple(P), tuple(cp(p, N[0]) for p in P))
    tets_out = np.concatenate([whole[~(sign[whole] == 0).all(1)] if len(whole) else whole.reshape(0, 4),
                               np.asarray(out, int).reshape(-1, 4)])
    return tets_out, np.asarray(new_edges, int).reshape(-1, 2), np.asarray(new_t, float)


def _orient(x, tets):
    v = np.linalg.det(np.swapaxes(x[tets[:, 1:]] - x[tets[:, 0, None]], 1, 2))
    tets = tets.copy()
    tets[v < 0] = tets[v < 0][:, [0, 2, 1, 3]]
    return tets, np.abs(v) / 6


def fitted_layer_mesh(vertices, faces, spacing_m=DEFAULT_SPACING_M, name='segment', *,
                      thickness_m=None, vertex_depth_m=None, skin_faces=None, verbose=False):
    """The layer with a BODY-FITTED boundary, by isosurface stuffing on the voxel lattice.

    Two level sets: phi_s, the signed distance to the (capped) skin, inside positive; and
    phi_c = depth(x) - h(x), zero on the layer/core interface (`DepthSurface`: the segment-median
    rule measures depth to the capped mesh exactly as the voxel rule did; the local rule to the
    skin faces, with the depth map's own value).  The lattice is the voxel rule's occupied cells
    grown by one, so it covers the skin.  For each level set in turn: warp nodes near a cut onto
    the surface (STUFFING_ALPHA), then cut every tet and keep the tissue side.  Every vertex made
    on a surface is then projected onto it EXACTLY (closest point on the skin; `_level_set` for
    the interface), so the boundary is the skin and the declared depth, not a linear
    interpolation of either.  Interface vertices are the BASE (carried by the segment).
    """
    from .material_domains import WindingHierarchy
    sp = float(spacing_m)
    if not np.isfinite(sp) or sp <= 0:
        raise ValueError('Positive spacing required')
    vertices, faces = np.asarray(vertices, float), np.asarray(faces, int)
    depth = DepthSurface(vertices, faces, thickness_m=thickness_m, vertex_depth_m=vertex_depth_m,
                         skin_faces=skin_faces)
    vp = voxel_partition([(name, vertices, faces)], spacing_m=sp)
    cells = vp['cell_indices']
    shape = cells.max(axis=0) + 5
    occupied = np.zeros(tuple(shape), bool)
    occupied[cells[:, 0] + 2, cells[:, 1] + 2, cells[:, 2] + 2] = True
    grown = np.argwhere(ndimage.binary_dilation(occupied, structure=np.ones((3, 3, 3), bool))) - 2
    x, tets = _kuhn_lattice(grown, vp['origin_m'], sp)
    winding = WindingHierarchy(vertices, faces)

    def phi_skin(p):
        d, _, _ = closest_points(p, vertices, faces)
        inside = np.abs(winding(p)) > 0.5
        return np.where(inside, d, -d)

    def phi_core(p):
        d, _, h = depth.query(p)
        return d - h

    # --- pass 1: the skin -------------------------------------------------------------------
    ps = phi_skin(x)
    x, warped = _warp(x, ps, tets, np.ones(len(x), bool), STUFFING_ALPHA)
    if warped.any():
        x[warped] = closest_points(x[warped], vertices, faces)[1]
    ps[warped] = 0.
    tets, edges, t = _cut(tets, len(x), ps)
    if len(edges):
        cutx = x[edges[:, 0]] + t[:, None] * (x[edges[:, 1]] - x[edges[:, 0]])
        x = np.concatenate([x, closest_points(cutx, vertices, faces)[1]])
        ps = np.concatenate([ps, np.zeros(len(edges))])
    on_skin = ps == 0
    # --- pass 2: the interface (keep phi_c <= 0; skin vertices never move) ---------------------
    used = np.unique(tets)
    pc = np.full(len(x), np.nan)
    pc[used] = phi_core(x[used])
    keep = -np.nan_to_num(pc, nan=1.0)
    y, warped = _warp(x, keep, tets, ~on_skin, STUFFING_ALPHA)
    if warped.any():
        y[warped] = _level_set(y[warped], depth)
    x = y
    keep[warped] = 0.
    tets, edges, t = _cut(tets, len(x), keep)
    if len(edges):
        cutx = x[edges[:, 0]] + t[:, None] * (x[edges[:, 1]] - x[edges[:, 0]])
        # a cut on an edge lying IN the skin (where the interface meets a cap or the skin) is a point
        # of the skin: it is projected onto the skin, not onto the level set, and it is both
        born_on_skin = on_skin[edges].all(axis=1)
        proj = cutx.copy()
        if (~born_on_skin).any():
            proj[~born_on_skin] = _level_set(cutx[~born_on_skin], depth)
        if born_on_skin.any():
            proj[born_on_skin] = closest_points(cutx[born_on_skin], vertices, faces)[1]
        x = np.concatenate([x, proj])
        keep = np.concatenate([keep, np.zeros(len(edges))])
        on_skin = np.concatenate([on_skin, born_on_skin])
    on_interface = keep == 0
    used = np.unique(tets)
    remap = np.full(len(x), -1)
    remap[used] = np.arange(len(used))
    nodes, tets = x[used], remap[tets]
    base, on_skin = on_interface[used], on_skin[used]
    tets, volume = _orient(nodes, tets)
    flat = volume <= 0
    tets, volume = tets[~flat], volume[~flat]
    # islands carried by nothing are removed and counted, as in the voxel rule
    edges_all = np.concatenate([tets[:, [i, j]] for i in range(4) for j in range(i + 1, 4)])
    graph = sparse.coo_matrix((np.ones(len(edges_all)), (edges_all[:, 0], edges_all[:, 1])), shape=(len(nodes),) * 2)
    count, label = connected_components(graph, directed=False)
    anchored = np.zeros(count, bool)
    anchored[np.unique(label[base])] = True
    keep_tet = anchored[label[tets[:, 0]]]
    islands, island_tets = int(count - anchored.sum()), int((~keep_tet).sum())
    tets, volume = tets[keep_tet], volume[keep_tet]
    used = np.unique(tets)
    remap = np.full(len(nodes), -1)
    remap[used] = np.arange(len(used))
    nodes, base, on_skin, tets = nodes[used], base[used], on_skin[used], remap[tets]
    # quality: volume over the cube of the longest edge, relative to a lattice Kuhn tet's (1/6 / 3^1.5)
    edge = np.stack([np.linalg.norm(nodes[tets[:, i]] - nodes[tets[:, j]], axis=1)
                     for i in range(4) for j in range(i + 1, 4)], axis=1)
    quality = (volume / edge.max(axis=1) ** 3) / ((1 / 6) / 3 ** 1.5)
    skin_res = closest_points(nodes[on_skin], vertices, faces)[0] if on_skin.any() else np.zeros(0)
    d_if, _, h_if = depth.query(nodes[base & ~on_skin]) if (base & ~on_skin).any() else (np.zeros(0), None, np.zeros(0))
    fit = {'skin_nodes': int(on_skin.sum()), 'interface_nodes': int(base.sum()),
           'skin_residual_max_m': float(skin_res.max(initial=0.)),
           'interface_residual_max_m': float(np.abs(d_if - h_if).max(initial=0.)),
           'flat_tetrahedra_removed': int(flat.sum()), 'islands': islands, 'island_tetrahedra': island_tets,
           'minimum_volume_m3': float(volume.min()), 'minimum_quality': float(quality.min()),
           'quality_percentiles_1_5_50': [float(q) for q in np.percentile(quality, [1, 5, 50])],
           'alpha': STUFFING_ALPHA, 'alpha_basis': 'Labelle & Shewchuk 2007 alpha_long for BCC, borrowed'}
    if verbose:
        print('        stuffing:', fit, flush=True)
    return {'unanchored_islands': islands, 'unanchored_island_tetrahedra': island_tets,
            'nodes_m': nodes, 'tetrahedra': tets, 'base': base, 'on_skin': on_skin,
            'rigid_skin_nodes': int((on_skin & base).sum()), 'spacing_m': sp, 'surface': 'fitted',
            'depth_rule': depth.rule, 'declared_thickness_m': None if thickness_m is None else float(thickness_m),
            'occupied_cells': int(len(cells)), 'lattice_cells': int(len(grown)),
            'layer_volume_m3': float(volume.sum()), 'enclosed_volume_m3': float(len(cells)) * sp ** 3,
            'fit': fit, 'boundary_discretization_m': max(fit['skin_residual_max_m'], fit['interface_residual_max_m'])}


# ------------------------------------------------------------------------- the fast path
# method='fast' minimises the SAME Phi over the SAME box as method='newton', to the same
# projected-gradient tolerance.  What changes is only how each Newton step is computed:
#   * the gradient is one sparse product with a scatter matrix built once, not four np.add.at;
#   * the Hessian is never assembled: Hessian-vector products are formed per element from F;
#   * the linear system is solved by truncated preconditioned CG (stopping on negative
#     curvature), preconditioned by the REST stiffness factorised ONCE per layer and free set
#     and reused across every iteration, load step, call and pose (a rotation only rotates it).
# The factorisation is a deterministic function of the layer, so caching it cannot change a
# result: the battery clears the cache and demands the same bits.

class _Kernel:
    """Element data of one sub-mesh in the SEGMENT frame, built once and posed per call.

    Arrays are held component-major -- (3, 3, T) rather than (T, 3, 3) -- because numpy's batched
    3x3 products are ~10x slower than the same sums written over contiguous (T,) components
    (measured on 12,000 tets: 1.33 ms matmul, 0.14 ms component-major; docs/SOFT_BODY.md).
    """

    def __init__(self, reference_local, tets, mu, lam, free_dof):
        region = DeformableRegion(reference_local, tets, mu_pa=mu, lambda_pa=lam, density_kg_m3=1.0)
        self.tets, self.mu, self.lam = region.tets, region.mu, region.lam
        self.volumes, self.inverse_local = region.volumes, region.inverse
        n, t = len(reference_local), len(tets)
        # local force of node a, coordinate i, tet e sits at column (a*3 + i)*T + e
        rows = (self.tets.T[:, None, :] * 3 + np.arange(3)[None, :, None]).ravel()
        self.scatter = sparse.csr_matrix((np.ones(12 * t), (rows, np.arange(12 * t))), shape=(3 * n, 12 * t))
        self.corner = [np.ascontiguousarray(self.tets[:, a]) for a in range(4)]
        # Ds[i, k, e] = y[tet_e[k+1], i] - y[tet_e[0], i] as ONE sparse product with the flattened positions
        e = np.arange(t)
        rows, cols, vals = [], [], []
        for i in range(3):
            for k in range(3):
                r = (i * 3 + k) * t + e
                rows += [r, r]
                cols += [self.tets[:, k + 1] * 3 + i, self.tets[:, 0] * 3 + i]
                vals += [np.ones(t), -np.ones(t)]
        self.gather = sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                                        shape=(9 * t, 3 * n))
        self.free = np.asarray(free_dof, bool).ravel()
        rest = assemble(region, element_hessians(region, region.reference, project=False))
        index = np.flatnonzero(self.free)
        kf = rest[index][:, index].tocsc()
        diagonal = kf.diagonal()
        kf = kf + sparse.diags(np.where(diagonal > 0, 0., max(float(diagonal.max()), 1.) * 1e-12))
        self.factor = splu(kf, permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0., options=dict(SymmetricMode=True))
        self.factor_nnz = int(self.factor.L.nnz + self.factor.U.nnz)
        self.n = n

    def posed(self, rotation):
        return _Posed(self, np.asarray(rotation, float))


def _mm(a, b):
    """(3,3,T) @ (3,3,T), component-major."""
    out = np.empty(a.shape)
    for i in range(3):
        for k in range(3):
            out[i, k] = a[i, 0] * b[0, k] + a[i, 1] * b[1, k] + a[i, 2] * b[2, k]
    return out


def _mmt(a, b):
    """a @ b^T, component-major."""
    out = np.empty(a.shape)
    for i in range(3):
        for k in range(3):
            out[i, k] = a[i, 0] * b[k, 0] + a[i, 1] * b[k, 1] + a[i, 2] * b[k, 2]
    return out


class _Posed:
    """`DeformableRegion`'s energy and gradient (same formulas, verified against it: battery S1) and the
    Hessian-vector product of `element_hessians(project=False)` (battery S2), on a posed kernel."""

    def __init__(self, kernel, rotation):
        self.k, self.rotation = kernel, rotation
        inverse = kernel.inverse_local @ rotation.T               # Dm_world = R Dm_local
        self.inverse = inverse                                    # (T,3,3), for callers that want it
        self.minv = np.ascontiguousarray(inverse.transpose(1, 2, 0))
        self.mu, self.lam, self.volumes, self.tets = kernel.mu, kernel.lam, kernel.volumes, kernel.tets

    def _grad(self, y):
        """F = Ds Dm^-1, component-major, from (N,3) positions."""
        ds = (self.k.gather @ np.asarray(y, float).ravel()).reshape(3, 3, -1)          # ds[i, k]
        return _mm(ds, self.minv)

    def deformation(self, positions):
        return np.ascontiguousarray(self._grad(positions).transpose(2, 0, 1))

    def _scatter(self, p):
        h = _mmt(p, self.minv) * self.volumes                    # h[i, k] = V sum_j P[i, j] Minv[k, j]
        local = np.concatenate([-h.sum(axis=1)[None], h.transpose(1, 0, 2)], axis=0)   # (4, 3, T)
        return (self.k.scatter @ local.ravel()).reshape(-1, 3)

    @staticmethod
    def _det_cof(f):
        cof = np.empty(f.shape)
        for i in range(3):
            for j in range(3):
                cof[i, j] = (f[(i + 1) % 3, (j + 1) % 3] * f[(i + 2) % 3, (j + 2) % 3]
                             - f[(i + 1) % 3, (j + 2) % 3] * f[(i + 2) % 3, (j + 1) % 3])
        det = f[0, 0] * cof[0, 0] + f[0, 1] * cof[0, 1] + f[0, 2] * cof[0, 2]
        return det, cof

    def energy_gradient(self, positions, *, trial=False):
        f = self._grad(positions)
        j, cofactor = self._det_cof(f)
        if not trial and np.any(j <= 0):
            raise ValueError('Inverted mechanical element')
        floor = .2
        safe = np.maximum(j, floor) if trial else j
        logj = np.log(safe)
        volumetric = -self.mu * logj + self.lam / 2 * logj ** 2
        derivative = (-self.mu + self.lam * logj) / safe
        if trial:
            delta = np.minimum(j - floor, 0.)
            second = (self.mu + self.lam * (1 - np.log(floor))) / floor ** 2
            volumetric += derivative * delta + .5 * second * delta ** 2
            derivative += second * delta
        density = self.mu / 2 * (np.sum(f * f, axis=(0, 1)) - 3) + volumetric
        piola = self.mu * f + derivative * cofactor
        # np.sum, not `@`: a float64 vector dot goes to OpenBLAS, whose thread start-up on this shared
        # 20-core machine costs ~5 ms for a 24,000-element dot against 7 us for np.sum (measured)
        return float(np.sum(self.volumes * density)), self._scatter(piola)

    def energy_difference(self, positions, delta):
        """E(positions + delta) - E(positions), computed from dF itself so that no two large totals are
        subtracted: |F+dF|^2 - |F|^2 = dF:(2F + dF), det(F+dF) - det F = cof F:dF + cof dF:F + det dF
        (exact for 3x3), and ln J1 - ln J0 = log1p(dJ / J0).  Near convergence the energy change of a
        Newton step falls below the resolution of the total (the reason `_projected_newton` stagnates
        at a floating-point floor); this difference does not.  Returns inf if an element inverts."""
        f0 = self._grad(positions)
        df = self._grad(delta)
        j0, cof0 = self._det_cof(f0)
        dj_, cofd = self._det_cof(df)
        dj = np.sum(cof0 * df, axis=(0, 1)) + np.sum(cofd * f0, axis=(0, 1)) + dj_
        ratio = dj / j0
        if np.any(ratio <= -1):
            return np.inf
        dlog = np.log1p(ratio)
        logj0 = np.log(j0)
        dvol = dlog * (-self.mu + self.lam / 2 * (2 * logj0 + dlog))
        ddev = self.mu / 2 * np.sum(df * (2 * f0 + df), axis=(0, 1))
        return float(np.sum(self.volumes * (ddev + dvol)))

    def prepare(self, positions):
        """Per-iterate data the Hessian-vector product needs (F^-T and c = lam ln J - mu)."""
        f = self._grad(positions)
        j, cofactor = self._det_cof(f)
        self._inv_t = cofactor / j                                # F^-T = cof(F) / det F
        self._c = self.lam * np.log(j) - self.mu

    def hessian_vector(self, v):
        """H v with H[a m, b n] = V [mu (B_a.B_b) d_mn + lam g_a[m] g_b[n] - c g_b[m] g_a[n]] (element_hessians'
        own formula, unprojected), without forming H: for Gv = sum_b v_b B_b^T and A = F^-T Gv^T,
        (H v)_a = V [mu Gv + lam tr(A) F^-T - c A F^-T] B_a."""
        gv = self._grad(v)
        a = _mmt(self._inv_t, gv)
        p = self.mu * gv + (self.lam * (a[0, 0] + a[1, 1] + a[2, 2])) * self._inv_t - self._c * _mm(a, self._inv_t)
        return self._scatter(p)

    def precondition(self, r):
        """Apply the rest stiffness inverse, rotated into the world: z = R K0^-1 R^T r on the free DOF."""
        local = (np.asarray(r).reshape(-1, 3) @ self.rotation).ravel()
        z = np.zeros_like(local)
        z[self.k.free] = self.k.factor.solve(local[self.k.free])
        return (z.reshape(-1, 3) @ self.rotation.T)

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
              gtol=1e-12, return_positions=False, warm_start_local_m=None):
        """Equilibrium (or one backward-Euler step) of the layer at a rigid segment pose.

        The support is the half-space `plane_sign * x[plane_axis] >= plane_sign * plane_value_m`.
        `state` is `None` (rest, at rest) or a previous result's `state`; it is only read.
        `held` is an optional (N, 3) bool of extra DOF carried rigidly with the segment
        (rollers, symmetry planes); their reactions count as the segment's.
        `method` is 'newton' (projected Newton on the analytic Hessian, the default) or
        'lbfgsb' (DeformableRegion's own minimiser, kept as the independent cross-check).
        'fast' minimises the same Phi to the same tolerance with Newton-CG preconditioned by the
        cached rest stiffness (see `_Kernel`); it takes no `held` and no crop.
        `warm_start_local_m` (N, 3), optional, is a start in the SEGMENT frame -- a previous
        result's `state['positions_local_m']` -- so a caller stepping a pose can begin from the last
        step's shape.  It is an argument like any other: the result is still a function of the
        arguments, and it converges to the same tolerance as a cold start.
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
                           'crop_held_nodes': 0, 'force_tolerance_n': 0.0, 'method': method})
            result['state'] = {'positions_m': rest, 'velocities_m_s': np.zeros_like(rest),
                               'positions_local_m': self.local.copy()}
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
        if method == 'fast':
            if held is not None or crop_radius_m is not None:
                raise ValueError("method='fast' takes no held DOF and no crop")
            sub = self._kernel(active_tets, nodes, index[tets], free_dof[nodes]).posed(rotation)
        else:
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
        if warm_start_local_m is not None:
            warm = np.asarray(warm_start_local_m, float)
            if warm.shape != rest.shape or not np.isfinite(warm).all():
                raise ValueError('warm_start_local_m must be (N, 3) and finite')
            if dynamic:
                raise ValueError('A dynamic step starts from its inertial prediction; no warm start')
            unclipped = np.where(fixed, 0., (warm @ np.asarray(rotation, float).T + origin)[nodes] - ref)
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
        elif method in ('newton', 'fast'):
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
                inner = self._projected_newton if method == 'newton' else self._projected_newton_cg
                trial, used, success, message = inner(
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
        receipt['state'] = {'positions_m': full, 'velocities_m_s': new_velocity,
                            'positions_local_m': (full - origin) @ np.asarray(rotation, float)}
        return receipt

    def _kernel(self, active_tets, nodes, tets, free_dof):
        """The cached `_Kernel` for this sub-mesh and free set (a deterministic function of both)."""
        cache = self.__dict__.setdefault('_kernels', {})
        key = hashlib.sha256(np.packbits(active_tets).tobytes() + np.packbits(free_dof).tobytes()).hexdigest()
        if key not in cache:
            cache[key] = _Kernel(self.local[nodes], tets, self.mu[active_tets], self.lam[active_tets], free_dof)
        return cache[key]

    def clear_cache(self):
        self.__dict__.pop('_kernels', None)

    @staticmethod
    def _projected_newton_cg(sub, ref, start, lo, hi, fixed, external, predicted, inertia, tolerance, maxiter):
        """`_projected_newton`'s iteration, active set and line search exactly; only the step differs:
        truncated CG on the free DOF, preconditioned by the cached rest factorisation, stopping on
        negative curvature (Steihaug) or at a relative residual of min(0.5, sqrt(|r|/|r0|)) (the
        Eisenstat-Walker forcing term, choice 2's square-root form), which keeps Newton's local
        superlinear rate without solving each step exactly."""
        dynamic = inertia is not None

        def phi(u):
            e, g = sub.energy_gradient(ref + u)
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
        first = None
        cg_total = 0
        for iteration in range(1, maxiter + 1):
            at_lo = u <= lo + 1e-14 * span
            at_hi = u >= hi - 1e-14 * span
            active = fixed | (at_lo & (g > 0)) | (at_hi & (g < 0))
            projected = np.where(active, 0., g)
            level = float(np.max(np.abs(projected)))
            if level <= tolerance:
                return u, iteration - 1, True, 'projected gradient below %.3g N (%d CG)' % (tolerance, cg_total)
            first = level if first is None else first
            sub.prepare(ref + u)
            keep = ~active

            def hv(v):
                out = sub.hessian_vector(v)
                if dynamic:
                    out = out + inertia * v
                return np.where(keep, out, 0.)

            r = -projected
            z = np.where(keep, sub.precondition(r), 0.)
            p = z.copy()
            x = np.zeros_like(u)
            rz = float(np.sum(r * z))
            r0 = float(np.sqrt(np.sum(r * r)))
            eta = min(0.5, np.sqrt(level / first))
            for _ in range(200):
                hp = hv(p)
                curvature = float(np.sum(p * hp))
                if curvature <= 0:
                    if not x.any():
                        x = p                     # the preconditioned steepest descent is a descent direction
                    break
                alpha = rz / curvature
                x = x + alpha * p
                r = r - alpha * hp
                cg_total += 1
                if float(np.sqrt(np.sum(r * r))) <= eta * r0:
                    break
                z = np.where(keep, sub.precondition(r), 0.)
                rz_new = float(np.sum(r * z))
                p = z + (rz_new / rz) * p
                rz = rz_new
            step = x
            if TRACE is not None:
                TRACE.append({'iteration': iteration, 'projected_n': level, 'active': int(active.sum() - fixed.sum()),
                              'cg': _ + 1, 'curvature_stop': curvature <= 0})
            if not float(np.sum(step * projected)) < 0:
                return u, iteration, False, 'no descent direction at projected gradient %.3g N' % level
            alpha = 1.0
            while True:
                trial = np.clip(u + alpha * step, lo, hi)
                change = trial - u
                # the same Armijo test on the same Phi, with the change computed accurately
                delta = sub.energy_difference(ref + u, change)
                delta -= float(np.sum(external * change))
                if dynamic:
                    delta += 0.5 * float(np.sum(inertia * change * (trial + u - 2 * predicted)))
                if np.isfinite(delta) and delta <= 1e-4 * float(np.sum(g * change)) and delta < 0:
                    try:
                        trial_value, trial_g = phi(trial)
                        break
                    except ValueError:
                        pass
                alpha *= 0.5
                if alpha < 1e-12:
                    return (u, iteration, level <= 100 * tolerance,
                            'stagnated at the floating-point floor, projected gradient %.3g N' % level)
            if TRACE is not None:
                TRACE[-1]['alpha'] = alpha
            u, value, g = trial, trial_value, trial_g
        return u, maxiter, False, 'Newton iteration limit reached'

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
                    # the Hessian is symmetric: a symmetric fill-reducing order and no row pivoting
                    # (SuperLU's SymmetricMode).  Was COLAMD: 0.21 s and 4.74 M factor entries per
                    # heel iteration against 0.14 s and 3.36 M (docs/SOFT_BODY.md, speed).  A pivot
                    # breakdown on the unprojected Hessian falls through to the projected one below.
                    candidate = -splu(kf.tocsc(), permc_spec='MMD_AT_PLUS_A', diag_pivot_thresh=0.,
                                      options=dict(SymmetricMode=True)).solve(rhs)
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


def local_depth(root, body, record, bundle=BUNDLE):
    """(vertex_depth_m, skin_face_mask) for one segment, refused unless it was built from THIS mesh."""
    path = Path(root) / LOCAL_DEPTH / 'manifest.json'
    if not path.exists():
        raise ValueError('No local depth artefact; run scripts/build_soft_tissue_local_depth.py')
    manifest = json.loads(path.read_bytes())
    if manifest.get('schema') != LOCAL_DEPTH_SCHEMA or manifest.get('bundle') != bundle:
        raise ValueError('Local depth artefact is not for this bundle')
    entry = manifest['bodies'].get(body)
    if entry is None or entry['mesh_sha256'] != record['written_sha256']:
        raise ValueError('Local depth was not built from this skin mesh: ' + str(body))
    if not (entry['A1_max_gap_m'] <= 1e-6 and entry['A2_median_equals_declared']):
        raise ValueError('Local depth failed its own known answers: ' + str(body))
    depth = np.array([np.nan if x is None else x for x in entry['vertex_depth_m']], float)
    return depth, ~np.asarray(entry['cap_face'], bool)


def segment_layer(root, body, *, bundle=BUNDLE, spacing_m=DEFAULT_SPACING_M, mapping='confined',
                  surface='voxel', depth='segment_median'):
    """Build one segment's layer from the measured layer map, and say what it rests on.

    `surface='voxel', depth='segment_median'` is the layer as first built (commit a1b2812), bit
    for bit.  `surface='fitted'` puts the boundary nodes on the skin and on the declared depth
    (`fitted_layer_mesh`); `depth='local'` takes the depth map's own value under each point of
    the skin instead of the segment's median (fitted surface only).
    """
    from .tissue_materials import DENSITY
    if surface not in SURFACES or depth not in DEPTH_RULES:
        raise ValueError('Unknown surface or depth rule')
    if surface == 'voxel' and depth != 'segment_median':
        raise ValueError('The local depth rule is built on the fitted surface only')
    manifest, record, mesh_path = _verified_record(root, bundle, body)
    layer = record['layer']
    poisson = float(manifest['skin_material']['poissons_ratio'])
    mu, lam, young = lame_from_layer_modulus(layer['apparent_modulus_pa'], poisson, mapping)
    rho = DENSITY['adipose'][0]
    vertices, faces = load_obj(mesh_path)
    if surface == 'voxel':
        mesh = layer_mesh(vertices, faces, layer['thickness_m'], spacing_m, name=body)
    elif depth == 'segment_median':
        mesh = fitted_layer_mesh(vertices, faces, spacing_m, name=body, thickness_m=layer['thickness_m'])
    else:
        vertex_depth, skin_faces = local_depth(root, body, record, bundle)
        mesh = fitted_layer_mesh(vertices, faces, spacing_m, name=body, vertex_depth_m=vertex_depth,
                                 skin_faces=skin_faces)
    meta = {'body': body, 'bundle': bundle, 'mesh_sha256': record['written_sha256'],
            'thickness_m': layer['thickness_m'], 'depth_points': layer['depth_points'],
            'apparent_modulus_pa': layer['apparent_modulus_pa'],
            'modulus_basis': layer['modulus_basis'], 'poisson_ratio': poisson,
            'poisson_basis': 'the bundle\'s declared skin_material poissons_ratio',
            'mapping': mapping, 'mapping_basis': LAYER_MODULUS_MAPPINGS[mapping],
            'young_modulus_pa': young, 'mu_pa': mu, 'lambda_pa': lam,
            'density_kg_m3': rho, 'density_basis': 'tissue_materials.DENSITY adipose (ICRU-44)',
            'layer_map_stiffness_pa_per_m': layer['stiffness_pa_per_m'],
            'surface': surface, 'depth_rule': depth}
    if surface == 'fitted':
        meta['fit'] = {k: v for k, v in (mesh.get('fit') or {}).items() if k != 'stages'}
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
    refused = sorted(set(wanted) & set(spec.get('refused_segments', ())))
    if refused:
        # never dropped silently: a caller asking for every segment is told which ones this identity cannot build
        raise ValueError('This identity refuses ' + ', '.join(refused) + '; pass `bodies` without them')
    return {body: segment_layer(root, body, bundle=spec['bundle'], spacing_m=spec['spacing_m'],
                                mapping=spec['mapping'], surface=spec.get('surface', 'voxel'),
                                depth=spec.get('depth', 'segment_median')) for body in wanted}


# ---------------------------------------------------------------------------------------
# COUPLING THE LAYER INTO THE PLANT'S INTEGRATION LOOP
# ---------------------------------------------------------------------------------------
# docs/ACTUATION_STAGES.md's "fully present participant" mode is this coupling: the layer
# deforms under contact and what it transmits reaches the thing that integrates.  Until now
# nothing it computed reached the plant (docs/WORKBENCH_AUTHENTICITY.md 2.1).
#
# THE COST THAT SHAPES IT.  A loaded, warm-started heel step costs 136-159 ms against the
# plant's 10 ms step (docs/SOFT_BODY.md, RT1h recorded FAILED).  One solve per plant step is
# therefore not available, and this class is the honest coupling that is:
#
#   1. ONLY WHERE IT MATTERS.  Every plant step, each coupled segment's layer is tested for
#      contact by one matrix-vector product over its own nodes (microseconds).  A segment
#      whose skin is clear of the support returns EXACTLY zero and emits NO force port, so a
#      plant with nothing in contact is the historical plant to the bit.  This is not an
#      approximation: it is the same predicate `solve` uses to decide there is nothing to do.
#   2. SUB-CYCLING.  A loaded segment is re-solved when `resolve_interval_s` has elapsed
#      since its last solve, and the reaction is HELD in between.  What is held is the world
#      force vector (the support's normal direction is fixed by the support) and the centre
#      of pressure as a station in the SEGMENT frame, so the station follows the segment and
#      the moment it delivers changes as the segment moves.  What staleness costs is measured
#      in scripts/verify_soft_tissue_coupled.py, and the cadence is chosen from that.
#   3. WARM STARTING.  Each solve starts from the last solve's shape in the segment frame,
#      which `SoftTissueLayer.solve` takes as an argument like any other.
#
# THE TWO GATES THAT ARE NOT SUB-CYCLED.  Contact onset and BOTTOMING OUT are checked every
# plant step, at full rate, because both are free.  A held reaction that sailed past the step
# where the rigid core entered the support would be a coupling that hid its own failure.
#
# WHAT IT REFUSES TO SWALLOW.  `solve` raises when the rigid core would enter the support, and
# when an element leaves the constitutive domain.  Both are re-raised here with the segment
# and the depth named; neither is caught by the caller in `articulated.py`, so the plant step
# does not happen and the plant is left exactly where it was.
#
# THE WRENCH, AND THE ONE COMPONENT A POINT FORCE CANNOT CARRY.  The engine's force port is a
# point force on a body (`NativeMechanicalStream.advance`), so the layer's (force, moment)
# pair is delivered as a force at its centre of pressure, p = (F x M)/|F|^2, which satisfies
# p x F = M EXACTLY when M is perpendicular to F.  Frictionless contact against a half-space
# puts every nodal contact force along the support normal, so M.F is zero up to the solve
# residual -- and that is checked on every emit rather than assumed.  A wrench with a real
# axial moment (friction, when it exists) CANNOT be delivered this way, and this raises
# instead of silently dropping it.
#
# NOT REAL-TIME BY CONSTRUCTION, and not converged: no force from this layer is good to
# better than about 10% (docs/SOFT_BODY.md CV7-CV10).  A cadence cannot fix either.

COUPLING_SCHEMA = 'ihm.soft-tissue-coupling.v1'


class SoftTissueBottomedOut(ValueError):
    """The rigid core of a coupled segment would enter the support.

    The tissue cannot carry the pose.  Raised BEFORE the plant step is taken, so the
    plant is left where it was; the caller decides what to do about it.
    """

    def __init__(self, body, depth_m):
        self.body = str(body)
        self.depth_m = float(depth_m)
        super().__init__('Bottomed out: the rigid core of %s penetrates the support by '
                         '%.4f m; the soft tissue cannot carry this pose' % (self.body, self.depth_m))


class SoftTissueLeftDomain(RuntimeError):
    """A coupled layer's solve left the constitutive domain (minimum J <= 0.2)."""

    def __init__(self, body, message):
        self.body = str(body)
        super().__init__('%s: %s' % (self.body, message))


class SoftTissueWrenchNotDeliverable(ValueError):
    """The layer's wrench cannot be delivered as a point force on the segment.

    Either it carries a moment about its own force axis (which a point force cannot
    produce), or its centre of pressure falls outside the segment's own layer.  Raised
    rather than dropping the part that does not fit.
    """


class SoftTissueCoupling:
    """Pose the coupled layers at the plant's own segment transforms and return force ports.

    `advance_state` is a PURE function of `(state, transforms, dt_s)`: it returns a new
    state and never touches the one it was given, and called twice at the same input it
    returns bitwise the same thing.  All state a caller has to carry is the returned dict,
    which is what makes checkpoint and rollback exact.

    Frames: `transforms` are the plant's own `bodies[*]['transform_ground']` (4x4, native
    SOURCE frame), the support is the half-space `sign * x[axis] >= sign * plane`, and the
    emitted ports are `{'body', 'point_m', 'force_n'}` in that same source frame.
    """

    def __init__(self, layers, *, plane_axis, plane_sign, plane_value_m, resolve_interval_s,
                 method='fast', force_tolerance_n=None, minimum_force_n=1e-9,
                 axial_moment_relative=1e-6):
        if not layers:
            raise ValueError('A coupling needs at least one layer')
        if plane_axis not in (0, 1, 2) or plane_sign not in (1, -1, 1.0, -1.0):
            raise ValueError('Axis-aligned half-space required')
        if not np.isfinite(plane_value_m) or not np.isfinite(resolve_interval_s) or resolve_interval_s <= 0:
            raise ValueError('Finite support plane and a positive resolve interval required')
        if method not in ('fast', 'newton'):
            raise ValueError("Coupled solves use method='fast' or 'newton'")
        self.layers = dict(layers)
        self.bodies = sorted(self.layers)
        self.axis = int(plane_axis)
        self.sign = float(plane_sign)
        self.plane = float(plane_value_m)
        self.interval_s = float(resolve_interval_s)
        self.method = str(method)
        self.force_tolerance_n = force_tolerance_n
        self.minimum_force_n = float(minimum_force_n)
        self.axial_moment_relative = float(axial_moment_relative)
        # the farthest any node of the layer sits from the segment origin: a centre of
        # pressure outside it is not a point of this segment's tissue
        self.radius_m = {b: float(np.linalg.norm(L.local, axis=1).max()) for b, L in self.layers.items()}

    # -- identity ----------------------------------------------------------------------
    def identity(self):
        return {'schema': COUPLING_SCHEMA, 'bodies': list(self.bodies),
                'support': {'axis': self.axis, 'sign': self.sign, 'plane_value_m': self.plane},
                'resolve_interval_s': self.interval_s, 'method': self.method,
                'force_tolerance_n': self.force_tolerance_n,
                'minimum_force_n': self.minimum_force_n,
                'axial_moment_relative': self.axial_moment_relative,
                'layer_dof': {b: int(L.dof) for b, L in self.layers.items()},
                'layer_nodes': {b: int(len(L.local)) for b, L in self.layers.items()},
                'layer_meta': {b: {k: v for k, v in L.meta.items() if k != 'fit'}
                               for b, L in self.layers.items()},
                'basis': 'Contact gate and bottoming-out gate every step; the SOLVE is '
                         'sub-cycled at resolve_interval_s and the reaction is held between '
                         'solves as a world force at a segment-frame station. The wrench is '
                         'delivered as a point force at the centre of pressure. No force from '
                         'this layer is converged to better than about 10% (docs/SOFT_BODY.md).'}

    def new_state(self):
        return {b: None for b in self.bodies}

    # -- the cheap per-step predicate ----------------------------------------------------
    def gaps(self, body, rotation, translation):
        """(skin gap, core gap) of this layer against the support, in metres, positive clear.

        One matrix-vector product over the layer's own nodes -- the same quantity
        `solve` computes to decide whether there is anything to do.
        """
        L = self.layers[body]
        r = np.asarray(rotation, float)
        t = np.asarray(translation, float)
        height = L.local @ r[self.axis, :3] + t[self.axis]
        gap = self.sign * (height - self.plane)
        return float(gap.min()), float(gap[L.base].min())

    # -- the step ------------------------------------------------------------------------
    def advance_state(self, state, transforms, dt_s):
        """(new state, force ports, receipts) for one plant step at `transforms`.

        `state` is not modified.  `dt_s` only ages the held reactions, so the cadence is a
        statement about TIME and not about how a caller happens to chop it up.
        """
        dt = float(dt_s)
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError('Positive finite step required')
        new_state, ports, receipts = {}, [], []
        for body in self.bodies:
            T = np.asarray(transforms[body], float)
            if T.shape != (4, 4) or not np.isfinite(T).all():
                raise ValueError('Finite 4x4 segment transform required: ' + body)
            rotation, origin = T[:3, :3], T[:3, 3]
            skin_gap, core_gap = self.gaps(body, rotation, origin)
            # NOT sub-cycled: the tissue saying it cannot carry the pose is checked every step
            if core_gap < 0:
                raise SoftTissueBottomedOut(body, -core_gap)
            if skin_gap >= 0:
                new_state[body] = None
                receipts.append({'body': body, 'in_contact': False, 'skin_gap_m': skin_gap,
                                 'core_gap_m': core_gap, 'solved': False, 'held_age_s': 0.0,
                                 'force_n': [0.0, 0.0, 0.0], 'moment_about_origin_nm': [0.0, 0.0, 0.0],
                                 'wall_seconds': 0.0})
                continue
            previous = state.get(body)
            due = previous is None or previous['age_s'] >= self.interval_s - 1e-15
            wall = 0.0
            if due:
                began = time.perf_counter()
                held = self._solve(body, rotation, origin, previous)
                wall = time.perf_counter() - began
            else:
                held = {k: v for k, v in previous.items()}
            emitted = self._emit(body, held, rotation, origin)
            if emitted is not None:
                ports.append(emitted)
            receipts.append({'body': body, 'in_contact': True, 'skin_gap_m': skin_gap,
                             'core_gap_m': core_gap, 'solved': bool(due),
                             'held_age_s': float(held['age_s']),
                             'depth_change_since_solve_m': float(-skin_gap - held['solved_depth_m']),
                             'force_n': held['force_world_n'].tolist(),
                             'moment_about_origin_nm': [] if emitted is None else
                                 np.cross(np.asarray(emitted['point_m']) - origin,
                                          np.asarray(emitted['force_n'])).tolist(),
                             'solved_moment_about_origin_nm': held['moment_nm'].tolist(),
                             'iterations': int(held['iterations']),
                             'contact_nodes': int(held['contact_nodes']),
                             'converged': bool(held['converged']),
                             'free_residual_n': float(held['free_residual_n']),
                             'balance_force_relative': float(held['balance_force_relative']),
                             'wall_seconds': float(wall)})
            held = {k: v for k, v in held.items()}
            held['age_s'] = float(held['age_s'] + dt)
            new_state[body] = held
        return new_state, ports, receipts

    def _solve(self, body, rotation, origin, previous):
        layer = self.layers[body]
        warm = None if previous is None else previous['positions_local_m']
        try:
            result = layer.solve(rotation=rotation, translation=origin, plane_axis=self.axis,
                                 plane_value_m=self.plane, plane_sign=self.sign,
                                 method=self.method, force_tolerance_n=self.force_tolerance_n,
                                 warm_start_local_m=warm)
        except ValueError as error:
            if 'Bottomed out' in str(error):
                raise SoftTissueBottomedOut(body, float('nan')) from error
            raise
        except RuntimeError as error:
            raise SoftTissueLeftDomain(body, str(error)) from error
        force = np.asarray(result['segment_force_n'], float)
        moment = np.asarray(result['segment_moment_nm'], float)
        magnitude = float(np.linalg.norm(force))
        station = np.zeros(3)
        if magnitude > self.minimum_force_n:
            axial = float(force @ moment) / magnitude
            scale = max(float(np.linalg.norm(moment)), magnitude * self.radius_m[body])
            if abs(axial) > self.axial_moment_relative * scale:
                raise SoftTissueWrenchNotDeliverable(
                    '%s: the layer transmits %.3e N.m about its own force axis, which a point '
                    'force cannot deliver (%.3e of the wrench scale, bar %.1e)'
                    % (body, axial, abs(axial) / scale if scale else float('inf'),
                       self.axial_moment_relative))
            lever = np.cross(force, moment) / magnitude ** 2
            if float(np.linalg.norm(lever)) > self.radius_m[body]:
                raise SoftTissueWrenchNotDeliverable(
                    '%s: the centre of pressure is %.4f m from the segment origin, outside the '
                    'layer\'s own %.4f m' % (body, float(np.linalg.norm(lever)), self.radius_m[body]))
            station = np.asarray(rotation, float).T @ lever
        return {'force_world_n': force, 'moment_nm': moment, 'station_local_m': station,
                'zero': magnitude <= self.minimum_force_n,
                'positions_local_m': result['state']['positions_local_m'],
                'age_s': 0.0, 'solved_depth_m': -self.gaps(body, rotation, origin)[0],
                'iterations': int(result['iterations']),
                'contact_nodes': int(result['contact_nodes']),
                'converged': bool(result['converged']),
                'free_residual_n': float(result['free_residual_n']),
                'balance_force_relative': float(result['balance_force_relative'])}

    def _emit(self, body, held, rotation, origin):
        if held['zero']:
            return None
        point = np.asarray(origin, float) + np.asarray(rotation, float) @ held['station_local_m']
        return {'body': body, 'point_m': point.tolist(), 'force_n': held['force_world_n'].tolist()}
