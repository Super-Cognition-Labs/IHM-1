"""The contact law the deformable layer implies: an ELASTIC FOUNDATION fitted to it.

WHY THIS EXISTS.  `ihm/assembly/soft_tissue_layer.py` builds a 3-D neo-Hookean layer per
segment and `SoftTissueCoupling` puts it in the plant's loop at one solve per step
(docs/SOFT_BODY.md).  That coupling is EXPLICIT: the reaction is computed at the pose the
step starts from and held through it.  The same document measures why that cannot be made
accurate by running it more often -- the reaction grows at **1,472.86 N/s**, so a 10% bar
allows a hold of **1.0 ms**, a tenth of the plant step.  One solve per step is already 10x
outside its own bar, and it costs **53.9 ms per 10 ms step** on one segment.

So the layer does not belong inside the loop; a LAW derived from it does.  The engine
already carries an implicit-friendly contact element -- `ElasticFoundationForce` over a
`ContactMesh`, solved inside the error-controlled integrator, parameterised by a Young
modulus, a Poisson ratio and a layer thickness.  This module is the bridge: the layer is
used OFFLINE as the ground truth, and the foundation's single stiffness is fitted to it.

WHAT IS AND IS NOT FITTED.  `ElasticFoundationForce` has exactly one elastic parameter per
contact mesh -- a stiffness in Pa/m -- and the repo's own reading of it
(`ihm/native/mechanical_stream.py`) is `k = (1-v)E/((1+v)(1-2v)h)`.  A fit therefore moves
ONE number per segment.  There is no second elastic degree of freedom to spend: a
"per-segment thickness/modulus map" is the same k written differently, which is why
`material_from_stiffness` below exists and why this module never claims to have fitted E and
h separately.  The dissipation, the friction coefficients and the transition velocity are
NOT touched: they are the shipped engineering constants, and a quasi-static layer carries no
information about any of them.

THE LAW, READ OUT OF THE ENGINE'S OWN SOURCE, NOT REMEMBERED.
`data/raw/mechanics/simbody/Simbody/src/ElasticFoundationForce.cpp`:

  * `setBodyParameters` (lines 76-88): `springPosition[i]` is the face's CENTROID and
    `springArea[i]` is `mesh.getFaceArea(i)`.
  * `processContact` (lines 143-192): for each face the broadphase returned, the nearest
    point on the other object is found; if the centroid is not INSIDE the other object the
    face is skipped; otherwise `f = stiffness * area * distance * (1 + dissipation*vnormal)`
    along the displacement direction, applied at the NEAREST POINT (not the centroid).
  * `areaScale` (line 107) is 1 unless BOTH surfaces carry parameters.  OpenSim's
    `ElasticFoundationForce::extendAddToSystem`
    (`data/raw/mechanics/opensim-core/.../ElasticFoundationForce.cpp:86-92`) registers
    parameters only for a `ContactMesh`, and the floor is a `ContactHalfSpace`, so the
    plant's own contacts always run at `areaScale = 1`.

Against an axis-aligned half-space the nearest point is the perpendicular projection, so the
distance is the centroid's penetration depth and the direction is the support normal.  The
force is therefore LINEAR IN THE STIFFNESS at a fixed pose:

    F(pose, k) = k * sum_{faces with centroid inside} area_f * depth_f * normal

which is why `foundation_geometry` below returns the k = 1 wrench and every fit here is a
one-parameter linear fit against a quantity that depends only on geometry.

WHAT THIS IS NOT.  It is a FIT to a layer that is itself converged to no better than about
10% (docs/SOFT_BODY.md CV7-CV10: monotone, first order, Richardson puts the finest heel
forces 8% and 15% above their limits).  Every number produced here inherits that bar and
says so.  It is not the layer, it does not deform, it has no lateral bulge, it shares no
load between neighbouring springs, and it is frictionless only because the layer it was
fitted to was.
"""
import numpy as np

# The elastic foundation's own reading of its stiffness, and the repo's: a uniform elastic
# LAYER of thickness h over a rigid substrate, loaded much wider than it is thick, has the
# confined (oedometric) modulus E(1-v)/((1+v)(1-2v)) and therefore this stiffness per metre.
# Identical expression to ihm/native/mechanical_stream.py and to the bundle manifest's
# `skin_material.stiffness_basis`; neither is edited from here.
STIFFNESS_BASIS = ("k=(1-v)E/((1+v)(1-2v)h) -- the elastic foundation's own law for a uniform "
                   "elastic layer of thickness h over a rigid substrate")

FOUNDATION_BASIS = (
    'SimTK::ElasticFoundationForce: one independent spring at every triangle CENTROID whose '
    'centroid lies inside the other object, f = k * area * depth, applied at the nearest '
    'point on the support. Read from data/raw/mechanics/simbody/Simbody/src/'
    'ElasticFoundationForce.cpp (setBodyParameters, processContact), not from memory.')

# The layer is quasi-static: `SoftTissueLayer.solve` takes no velocity in its default mode,
# so it carries NO information about a rate term.  The foundation's rate term is its
# dissipation coefficient, which is a shipped engineering constant.  A fit made from a
# quasi-static layer therefore calibrates the ELASTIC term and nothing else, and saying so
# is part of the disclosure rather than a footnote.
RATE_BASIS = ('The layer is quasi-static, so it carries no rate information and the '
              "foundation's dissipation coefficient is NOT fitted: it stays the shipped "
              'engineering constant. Only the elastic stiffness is calibrated.')

OBJECTIVES = ('log', 'linear')
OBJECTIVE_BASIS = {
    'log': 'least squares on log(force): every training penetration carries equal RELATIVE '
           'weight, which is the unit the substitution cost is reported in. The minimiser is '
           'the geometric mean of the per-point effective stiffness.',
    'linear': 'least squares on force: the deepest training penetration dominates. Reported '
              'beside the log fit as a sensitivity, never as the headline.',
}


def stiffness_from_material(youngs_modulus_pa, poissons_ratio, layer_thickness_m):
    """The elastic foundation's stiffness per metre, from E, v and h."""
    e = float(youngs_modulus_pa)
    v = float(poissons_ratio)
    h = float(layer_thickness_m)
    if not np.isfinite([e, v, h]).all() or e <= 0 or h <= 0:
        raise ValueError('Positive finite modulus and thickness required')
    if not 0 <= v < 0.5:
        raise ValueError('Poisson ratio must lie in [0, 0.5)')
    return (1 - v) * e / ((1 + v) * (1 - 2 * v) * h)


def material_from_stiffness(stiffness_pa_per_m, poissons_ratio, layer_thickness_m):
    """The Young modulus a stiffness implies at a DECLARED v and h -- the exact inverse.

    There is one elastic parameter and three names for it.  Holding v and h at the values the
    bundle declares is what makes the answer a modulus rather than a free choice, and it is
    stated at every call site rather than assumed.
    """
    k = float(stiffness_pa_per_m)
    v = float(poissons_ratio)
    h = float(layer_thickness_m)
    if not np.isfinite([k, v, h]).all() or k <= 0 or h <= 0:
        raise ValueError('Positive finite stiffness and thickness required')
    if not 0 <= v < 0.5:
        raise ValueError('Poisson ratio must lie in [0, 0.5)')
    return k * (1 + v) * (1 - 2 * v) * h / (1 - v)


def _posed(vertices, faces, rotation, translation):
    v = np.asarray(vertices, float)
    f = np.asarray(faces)
    r = np.asarray(rotation, float)
    t = np.asarray(translation, float)
    if v.ndim != 2 or v.shape[1] != 3 or f.ndim != 2 or f.shape[1] != 3:
        raise ValueError('An (N,3) vertex array and an (M,3) face array are required')
    if r.shape != (3, 3) or t.shape != (3,) or not np.isfinite(r).all() or not np.isfinite(t).all():
        raise ValueError('Finite 3x3 rotation and 3-vector translation required')
    if not np.allclose(r @ r.T, np.eye(3), atol=1e-9) or np.linalg.det(r) < 0.999999:
        raise ValueError('Proper rotation required')
    return v @ r.T + t, f


def foundation_geometry(vertices, faces, *, rotation=np.eye(3), translation=np.zeros(3),
                        plane_axis, plane_value_m, plane_sign=1.0, about=None):
    """The k = 1 elastic-foundation wrench of one contact mesh against a half-space.

    The support is `plane_sign * x[plane_axis] >= plane_sign * plane_value_m`, the same
    convention `SoftTissueLayer.solve` and `SoftTissueCoupling` use, so a pose can be handed
    to both without a frame change.

    Returns `integral_m3 = sum area_f * depth_f` over the faces whose CENTROID is inside the
    support, the contacting face count and area, and the unit wrench `(force, moment about
    `about`, default the segment origin `translation`)`.  The actual wrench at stiffness k is
    exactly k times this, because the law is linear in k at a fixed pose.

    A pure function of its arguments: no cache, no generator, no state.
    """
    if plane_axis not in (0, 1, 2):
        raise ValueError('Axis-aligned half-space required')
    if plane_sign not in (1, -1, 1.0, -1.0):
        raise ValueError('Axis-aligned half-space required')
    if not np.isfinite(float(plane_value_m)):
        raise ValueError('Finite support plane required')
    posed, f = _posed(vertices, faces, rotation, translation)
    axis = int(plane_axis)
    sign = float(plane_sign)
    plane = float(plane_value_m)
    reference = np.asarray(translation, float) if about is None else np.asarray(about, float)
    if reference.shape != (3,) or not np.isfinite(reference).all():
        raise ValueError('Finite 3-vector reference point required')

    tri = posed[f]
    centroid = tri.mean(axis=1)
    area = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
    # `distance` in the engine is |nearestPoint - springPosition|, and the engine SKIPS a
    # face whose centroid is not inside the other object and a face at distance exactly
    # zero.  Both are reproduced here rather than approximated by `maximum(0, ...)`.
    depth = sign * (plane - centroid[:, axis])
    inside = depth > 0
    normal = np.zeros(3)
    normal[axis] = sign
    weight = area[inside] * depth[inside]
    integral = float(weight.sum())
    # the force acts at the nearest point on the support: the centroid projected onto the
    # plane.  The projection moves the point ALONG the force direction, so the moment is
    # unchanged by it -- it is done anyway, because reproducing the engine is the point.
    station = centroid[inside].copy()
    station[:, axis] = plane
    force = integral * normal
    moment = np.cross(station - reference, weight[:, None] * normal[None, :]).sum(axis=0) \
        if inside.any() else np.zeros(3)
    # The contact FOOTPRINT, in the support plane.  It is reported because the confined
    # reading of the layer modulus (`soft_tissue_layer.LAYER_MODULUS_MAPPINGS`) is a
    # statement about a load much WIDER than the layer is thick, and the in-vivo source card
    # does not report the footprint-to-thickness ratio it was measured at.  The minor extent
    # is the one that decides it: a long narrow strip is not a wide load.
    footprint = (0.0, 0.0)
    if int(inside.sum()) >= 2:
        flat = np.delete(station, axis, axis=1)
        flat = flat - flat.mean(axis=0)
        basis = np.linalg.svd(flat, full_matrices=False)[2]
        spread = np.ptp(flat @ basis.T, axis=0)
        footprint = (float(spread.max()), float(spread.min()))
    return {'integral_m3': integral,
            'contact_faces': int(inside.sum()),
            'contact_area_m2': float(area[inside].sum()),
            'footprint_major_m': footprint[0],
            'footprint_minor_m': footprint[1],
            'total_area_m2': float(area.sum()),
            'maximum_depth_m': float(depth[inside].max()) if inside.any() else 0.0,
            'unit_force_n': force,
            'unit_moment_nm': moment,
            'reference_point_m': reference,
            'plane': {'axis': axis, 'sign': sign, 'value_m': plane},
            'basis': FOUNDATION_BASIS}


def foundation_resultant(vertices, faces, *, stiffness_pa_per_m, **pose):
    """The elastic foundation's force and moment on the segment at a pose.  `k` times the unit wrench."""
    k = float(stiffness_pa_per_m)
    if not np.isfinite(k) or k <= 0:
        raise ValueError('Positive finite foundation stiffness required')
    geometry = foundation_geometry(vertices, faces, **pose)
    return {**geometry, 'stiffness_pa_per_m': k,
            'force_n': k * geometry['unit_force_n'],
            'moment_nm': k * geometry['unit_moment_nm']}


def fit_stiffness(integral_m3, force_n, *, objective='log'):
    """The single foundation stiffness that best reproduces measured forces.

    `integral_m3` is `foundation_geometry(...)['integral_m3']` at each pose and `force_n` the
    LAYER's force magnitude at the same pose.  Both must be strictly positive: a pose the
    foundation does not touch, or one the layer does not load, carries no information about a
    stiffness and is refused rather than silently weighted as zero.

    `objective='log'` minimises sum (log k*G_i - log F_i)^2, whose minimiser is the geometric
    mean of F_i/G_i; `'linear'` minimises sum (k*G_i - F_i)^2, whose minimiser is
    sum(G_i F_i)/sum(G_i^2).  Both are closed form: there is no iteration, no seed and no
    generator anywhere in this module.
    """
    if objective not in OBJECTIVES:
        raise ValueError('Unknown fit objective: ' + str(objective))
    g = np.asarray(integral_m3, float).ravel()
    f = np.asarray(force_n, float).ravel()
    if g.shape != f.shape or g.size == 0:
        raise ValueError('One foundation integral per measured force is required')
    if not (np.isfinite(g).all() and np.isfinite(f).all()) or (g <= 0).any() or (f <= 0).any():
        raise ValueError('Every training point must have a positive foundation integral and a '
                         'positive measured force')
    if objective == 'log':
        return float(np.exp(np.mean(np.log(f / g))))
    return float((g @ f) / (g @ g))


def relative_errors(predicted_n, measured_n):
    """Signed relative error of a prediction against the layer, point by point."""
    p = np.asarray(predicted_n, float).ravel()
    m = np.asarray(measured_n, float).ravel()
    if p.shape != m.shape or p.size == 0:
        raise ValueError('One prediction per measurement is required')
    if (m == 0).any():
        raise ValueError('A relative error against a zero measurement is not defined')
    return (p - m) / m


def skill(predicted_n, measured_n, baseline_n):
    """1 - SSE(prediction)/SSE(baseline) in newtons, against an EXPLICIT baseline.

    A raw error is not a result (IBM-1 CLAUDE.md).  `baseline_n` is a scalar (predict zero,
    predict the mean) or a vector (another law's prediction on the same points).  1.0 is
    perfect, 0.0 is the baseline, negative is worse than the baseline.
    """
    p = np.asarray(predicted_n, float).ravel()
    m = np.asarray(measured_n, float).ravel()
    b = np.broadcast_to(np.asarray(baseline_n, float), m.shape)
    if p.shape != m.shape:
        raise ValueError('One prediction per measurement is required')
    denominator = float(np.sum((b - m) ** 2))
    if denominator == 0:
        return float('inf') if float(np.sum((p - m) ** 2)) > 0 else 0.0
    return 1.0 - float(np.sum((p - m) ** 2)) / denominator


def summarise(predicted_n, measured_n):
    """The substitution cost of one law on one set of poses, in newtons and as a fraction."""
    error = relative_errors(predicted_n, measured_n)
    absolute = np.asarray(predicted_n, float).ravel() - np.asarray(measured_n, float).ravel()
    m = np.asarray(measured_n, float).ravel()
    return {'n': int(error.size),
            'median_absolute_relative': float(np.median(np.abs(error))),
            'max_absolute_relative': float(np.abs(error).max()),
            'worst_at_force_n': float(m[int(np.argmax(np.abs(error)))]),
            'signed_relative': [float(v) for v in error],
            'max_absolute_n': float(np.abs(absolute).max()),
            'rms_n': float(np.sqrt(np.mean(absolute ** 2))),
            'skill_vs_zero': skill(predicted_n, measured_n, 0.0),
            'skill_vs_mean': skill(predicted_n, measured_n, float(np.mean(m)))}


# The disclosure a plant is handed with this law.  It is not an optional note: a caller that
# reads a force out of a plant carrying this bundle is reading a FIT to a layer that is
# itself good to about 10%, and the selection carries that sentence.
DISCLOSURE = (
    'The contact law the deformable soft-tissue layer implies. The elastic foundation over '
    "this segment's own skin mesh, with its stiffness FITTED OFFLINE to the 3-D neo-Hookean "
    'layer (ihm/assembly/soft_tissue_layer.py) at poses that actually load it, scored on '
    'penetrations the fit never saw. It is IMPLICIT -- solved inside the integrator, at no '
    'per-step cost beyond the segment contact the plant already carries -- and it is NOT the '
    'layer. THREE THINGS IT IS NOT. (1) It is a fit to a layer that is converged to no better '
    'than about 10 per cent (docs/SOFT_BODY.md CV7-CV10), so every fitted stiffness inherits '
    'that bar and no force from this law is good to better than it. (2) It does not deform: '
    'each spring is independent, there is no lateral bulge and no load sharing, which is the '
    'whole physical difference between a foundation and a continuum. (3) Only the ELASTIC '
    'term is fitted. The layer is quasi-static and carries no rate information, so the '
    "foundation's dissipation and friction stay the shipped engineering constants.")
