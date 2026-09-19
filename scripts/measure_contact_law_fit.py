#!/usr/bin/env python3
"""The contact law the soft-tissue layer implies: fit an ELASTIC FOUNDATION to the layer, and
measure what the substitution costs.

docs/SOFT_BODY.md put the 3-D neo-Hookean layer in the plant's loop and drew the conclusion
that shapes this file: the reaction grows at **1,472.86 N/s**, so a 10% bar allows a hold of
**1.0 ms** -- a tenth of the plant step.  One solve per plant step is therefore ALREADY 10x
outside its own bar, and it costs **53.9 ms per 10 ms step** on one segment.  An explicit
coupling cannot be fixed by running it more often; what it needs is to be IMPLICIT.

The engine already has an implicit contact element: `ElasticFoundationForce` over a
`ContactMesh`, solved inside the error-controlled integrator.  It has exactly ONE elastic
parameter per mesh, a stiffness in Pa/m, which `ihm/native/mechanical_stream.py` reads as
`k = (1-v)E/((1+v)(1-2v)h)`.  So: use the layer OFFLINE as the ground truth, fit that one
number per segment to it, and report honestly where the substitution tracks the layer and
where it departs.

    PYTHONPATH=. prlimit --as=4294967296 .venv/bin/python -u scripts/measure_contact_law_fit.py

================================================================================
PRE-REGISTRATION.  Everything below was written before any fit existed.  No bar here moves
afterwards, and a gate that fails is recorded FAILED.
================================================================================

WHAT WAS ALREADY MEASURED WHEN THIS WAS WRITTEN.  Three development probes preceded it; all
are disclosed, and NONE produced a fitted stiffness, an error or a bar:

  P1  An uncoupled 140-step collapse from `pelvis_ty = 1.03` puts `ulna_l`'s skin below the
      floor from step 81, reaching 23.07 mm by step 140.  20.0 s of wall clock, 391 MB.
      This reproduces docs/SOFT_BODY.md's own P3/P5 (contact at step 81) and is used ONLY to
      establish that fixture U has poses in it.
  P2  The `ulna_l` fitted-local layer builds in 2.7 s at 14,874 DOF and solves cold at those
      poses in 0.2-3.9 s, raising `Left the constitutive domain: minimum J 0.135 <= 0.2` at
      step 104.  Used ONLY to size this script.
  P3  The foundation INTEGRAL `sum area*depth` over `skin_calcn_l.obj` at the heel fixture,
      which is geometry and contains no force: 2.909e-08 m^3 at 0.5 mm rising to 7.310e-06
      at 6 mm, over 22 to 304 faces.  It reproduces the heel fixture's rotation,
      17.25522360129255 deg, to every printed digit of docs/SOFT_BODY.md's 17.2552, which is
      the only thing it was run for.

  A FOURTH thing was found and is a BUG, not a probe: `{'segment_contact': 'skin_layer_map'}`
  RAISED `Bundle carries a per-segment layer map; a caller E, p or h would override it` at
  plant construction.  `plant_options.resolve_fidelity` decided "per-record layer map?" from
  a top-level manifest key the bundle does not carry, while `NativeMechanicalStream` decides
  it from the records; so the resolver passed the uniform triple with every layer-map bundle
  and every such selection raised.  **The per-segment measured stiffness had never been
  reachable through the resolver.**  Fixed in `plant_options.py` (both sides now read the
  same records), gated in `scripts/verify_contact_law_fit.py`.

THE LAW BEING FITTED, read out of the engine's own source and not remembered
(`data/raw/mechanics/simbody/Simbody/src/ElasticFoundationForce.cpp`): one independent
spring at every triangle CENTROID that lies inside the support, `f = k * area * depth`,
applied at the nearest point on the support.  Against an axis-aligned half-space that makes
the wrench exactly linear in k at a fixed pose, so every fit below is a ONE-PARAMETER linear
fit against a purely geometric quantity, in closed form, with no seed and no generator.

THE FIXTURES.  Every one of them is a pose this repo already chose by measurement; none is
invented here.

  H  `calcn_l`, THE HEEL FIXTURE: the least rotation about the segment z axis that keeps
     every vertex of the forefoot half clear of the support at 6 mm, found by bisection on
     the skin mesh (`scripts/verify_soft_tissue.py:heel_pose`, 17.2552 deg).  docs/SOFT_BODY.md
     withdrew every earlier "heel" number for having pressed a 1.7 mm seam wedge at the
     calcn/toes cut instead; this is the fixture that replaced it.  Penetrations 0.5 to
     4.0 mm in 0.5 mm steps, measured from the ROTATED skin mesh's lowest point.  The 4.0 mm
     ceiling is the LAYER's, not a choice: the depth map puts 6.4-14.1 mm of tissue under the
     plantar heel and `verify_soft_tissue.py`'s S3h stops at 4 mm for that reason.
  S  `calcn_l`, THE FLAT STANCE FIXTURE (identity rotation), 0.5 to 2.5 mm in 0.5 mm steps.
     docs/SOFT_BODY.md's P2 measured this layer leaving its constitutive domain at 3.0 mm
     here, so 2.5 mm is the ceiling and it too is the layer's.
  U  `ulna_l`, THE PLANT'S OWN POSES.  A 140-step uncoupled collapse from `pelvis_ty = 1.03`,
     10 ms steps, no actuation: the fixture the coupling measurement chose BY MEASUREMENT
     (docs/SOFT_BODY.md P3 -- 69 steps in which its skin is below the floor while its own
     engine contact element is not, the longest such window in the body).  Every step whose
     skin is below the floor is a pose, and the layer is solved COLD at the plant's own
     transform, never warm-started from the previous step, so each measurement is a function
     of its own pose alone.

  **This is the 22-segment scaffold.  It is not the body, and nothing here is a statement
  about a human forearm or a human heel.**

THE LAYER IS THE GROUND TRUTH, at the identity that ships: `layer_fitted_local_confined`
(fitted surface, local depth, confined modulus reading) at its own 5 mm spacing, solved with
`method='fast'`.  It is converged to no better than about 10% (CV7-CV10: monotone, first
order, Richardson puts the finest heel forces 8% and 15% above their limits).  **Every number
this script produces inherits that bar and says so.**

THE SPLIT, fixed here, because a fit reported on its own training points is not a result.
  H, S:  TRAIN on the odd-numbered penetrations (0.5, 1.5, 2.5, 3.5 mm), SCORE on the even
         ones (1.0, 2.0, 3.0, 4.0 mm).  Interleaved, so the held-out points are
         INTERPOLATION; extrapolation is a different question and is asked separately.
  U:     TRAIN on the even-indexed loaded poses, SCORE on the odd-indexed ones.
  CROSS-FIXTURE, the hard one: the stiffness fitted on H alone is scored on EVERY point of S
         -- a different pose family of the same segment, which the fit never saw.
  A pose at which the LAYER carries exactly zero force, or at which the foundation touches no
  face, is dropped from both sides and REPORTED: it carries no information about a stiffness,
  and weighting it as zero would be a silent choice.

THE FIT.  One stiffness per segment, least squares on LOG force
(`contact_law.fit_stiffness(objective='log')`), so every training penetration carries equal
relative weight -- which is the unit the substitution cost is reported in.  The linear-force
fit is computed and printed beside it as a sensitivity and is never the headline.  Nothing
else is free: v and h stay the values the bundle declares, and the fitted k is reported as
the Young modulus it implies AT THOSE.

THE BASELINES, all four, on every scored set, because a raw error is not a result:
  B0  predict zero.
  B1  predict the mean of the layer's force over the scored points.
  B2  THE TRIVIAL ONE: the bundle's own UNFITTED per-segment stiffness, `k = E_app/h` from
      the layer map (`calcn_l` 1.0347e7 Pa/m, `ulna_l` 3.2826e6 Pa/m).
  B3  the uniform `skin_material` triple the bundle was built with -- E = 3000 Pa, v = 0.45,
      h = 6.6 mm, k = 1.7241e6 Pa/m -- which is the "shipped E/v/h triple, unfitted".

THE BARS, fixed here.
  F1  GATED.  The fit's MEDIAN absolute relative error on the HELD-OUT penetrations of its
      own fixture is <= 0.10.  10% is the layer's own convergence bar, and it is the only
      place this number comes from: a fit inside it is as good as the thing it was fitted to
      is known, and a tighter bar would be a claim about a number the layer does not have.
  F2  GATED.  The fit BEATS the trivial baseline B2 on the held-out points by more than the
      bar on BOTH sides -- sqrt(2) * 0.10 = 0.1414 of median relative error -- because B2 was
      measured against the same uncertain layer and a bare "beats the baseline" would be
      comparing one uncertain number with another.  If it does not, the honest result is that
      the shipped stiffness was already inside the bar and no fit was needed.  That is
      RECORDED as the outcome, never repaired.
  F3  RECORDED, never gated.  The cross-fixture error; and the per-pose EFFECTIVE stiffness
      `k_eff = F_layer / G`, with its spread.  A single k can follow the layer only in so far
      as k_eff is constant, so **the spread of k_eff IS the part of the layer the elastic
      foundation cannot represent**, and it is reported whatever it turns out to be.
  F4  RECORDED.  The layer is QUASI-STATIC, so it carries no rate information at all and the
      foundation's dissipation coefficient is NOT fitted.  Stated rather than measured,
      because there is nothing to measure: `SoftTissueLayer.solve` takes no velocity in the
      mode the coupling uses.

THREE OUTCOMES, all of them results, fixed before the data:
  (a) F1 and F2 pass -> the fitted law is offered as the layer's implied contact law.
  (b) F1 passes and F2 does not -> the foundation ALREADY tracks the layer at the shipped
      stiffness.  The fit is offered with its measured (small) improvement, and the FINDING
      is that the bundle's own k was right.
  (c) F1 fails -> the elastic foundation cannot represent this layer over this range.  The
      k_eff spread says by how much and where, and THAT is the deliverable.

================================================================================
AMENDMENT, written after run 1 and before run 2.  **NO BAR MOVES AND NO GATE IS RESCORED**:
run 1's verdicts stand (`logs/measure_contact_law_fit.run1.log`, F1H PASS, F2H FAILED, F1U
FAILED, F2U PASS) and run 2 must reproduce them gate for gate and number for number on
everything that existed in run 1.  Five changes, each with its reason:

  A1  THE BUNDLE IS ASSEMBLED BY THE PRE-REGISTERED OUTCOME RULE, which run 1 did not
      implement.  Outcome (c) above says that when F1 fails "the elastic foundation cannot
      represent this layer over this range ... and THAT is the deliverable" -- a measurement,
      not a law.  So a segment whose F1 FAILED now gets NO fitted stiffness: its record keeps
      the layer map's own k and says the fit was refused, with the numbers.  This is the rule
      as written, applied; it is not a choice made after seeing which segment failed.
  A2  A CORRECTION TO FIXTURE S.  As written -- identity rotation, plane measured from the
      skin mesh's lowest VERTEX -- it presses vertex 780 on the rim of the cut between
      `calcn_l` and the toes, which is the 1.7 mm seam wedge docs/SOFT_BODY.md WITHDREW every
      earlier "heel" number for.  It is not P2's flat stance pose, which is the PLANT's own
      `calcn_l` transform; citing P2 for it was wrong.  Run 1 measured the layer carrying
      exactly 0 N at 0.5-2.0 mm there and 0.0909 N at 2.5 mm, which reproduces the withdrawn
      fixture's own CV5/CV6 row (0/0/0/0.094/0.096 N).  **Fixture S and its 948% cross-fixture
      figure are KEPT and RECORDED as measured, and are not used to judge the contact law**:
      one loaded point on a seam wedge is a statement about the wedge.
  A3  ADDED, FIXTURE H2: the same heel derivation at a 3 mm forefoot clearance instead of
      6 mm (12.5921 deg for `calcn_l`, 12.7154 deg for `calcn_r`), which is a DIFFERENT pose
      family that still presses the heel.  It is scored by the H fit and never trained on.
      It replaces nothing: S stays.  The reason it is added is that a cross-fixture test must
      land on the part of the segment the question is about -- the same correction
      docs/SOFT_BODY.md made when it derived the heel fixture in the first place.
  A4  ADDED, SEGMENT `calcn_r`, on its own heel fixture, by the same split and the same bars.
      One segment is one segment; a heel result that holds on the other heel is worth more
      than one that has only been seen once.
  A5  ADDED, RECORDED, NEVER GATED: the contact FOOTPRINT (major and minor extent of the
      contacting faces in the support plane) at every pose, and its minor extent against the
      segment's declared layer thickness.  The confined reading of the layer modulus is a
      statement about a load much WIDER than the layer is thick and the in-vivo source card
      does not report the ratio it was measured at, so it is measured here.
  A6  ADDED, RECORDED, NEVER GATED: a REFINEMENT ARM for a hypothesis run 1 raised.  Run 1's
      k_eff outliers are exactly the poses where the LAYER's contact patch is resolved by a
      handful of nodes (`ulna_l` steps 81/82: **4 contact nodes** against 68 foundation
      faces; `calcn_l` at 0.5 mm: 11 nodes).  The hypothesis is that those k_eff values are
      the layer's own discretisation and not a property of the contact law.  A factorial
      names a factor and not a mechanism (CLAUDE.md), so the mechanism gets its own
      instrument: `ulna_l` is rebuilt at 4 mm and 3 mm spacing and re-solved at the same
      poses, and k_eff is reported at each.  If it moves toward the deep-pose value the
      hypothesis survives; if it does not, it is excluded.  **This changes no verdict**:
      F1U is FAILED and stays FAILED.

WHAT THIS SCRIPT WRITES.  `data/derived/contact-law-fit-v1/`: `report.json`, and a
selectable `ihm.segment-contact-meshes.v1` bundle (manifest + the same meshes, sha256
verified against the source bundle) whose per-segment stiffness is the fitted one where a
segment was fitted and the layer map's own where it was not -- marked per record, never
silently mixed.  `plant_options.SEGMENT_CONTACT_BUNDLES` names it.
"""
import hashlib
import json
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.assembly import contact_law as cl                                     # noqa: E402
from ihm.assembly import soft_tissue_layer as stl                              # noqa: E402

IDENTITY = 'layer_fitted_local_confined'
PLANE_AXIS, PLANE_SIGN = 1, 1.0        # OpenSim's y-up upright floor (plant_options.SOFT_TISSUE_SUPPORT)
FLOOR_M = 0.0                          # articulated.py: the upright support plane is y = 0
DEPTHS_HEEL = tuple(np.round(np.arange(0.5, 4.01, 0.5) * 1e-3, 12))
DEPTHS_FLAT = tuple(np.round(np.arange(0.5, 2.51, 0.5) * 1e-3, 12))
POSE = {'pelvis_ty': 1.03}
DT, STEPS = 0.01, 140
BAR_F1 = 0.10                          # the layer's own convergence bar
BAR_F2 = float(np.sqrt(2.0) * BAR_F1)  # a margin over a MEASURED baseline needs both sides
OUT = ROOT / 'data/derived/contact-law-fit-v1'
WORK = ROOT / 'data/derived/contact-law-fit-runs'

RESULTS = []


def gate(name, claim, ok, detail='', recorded=False):
    print('%-5s %-74s %s%s' % (name, claim, 'RECORDED' if recorded else ('PASS' if ok else 'FAILED'),
                               '' if not detail else '   ' + detail), flush=True)
    RESULTS.append({'gate': name, 'claim': claim, 'pass': bool(ok), 'detail': detail,
                    'recorded': bool(recorded)})
    return ok


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def bundle_mesh(body):
    manifest = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_bytes())
    record = [r for r in manifest['records'] if r['body'] == body][0]
    vertices, faces = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
    return manifest, record, vertices, faces


def heel_rotation(vertices, deepest_m=0.006):
    """scripts/verify_soft_tissue.py's own heel fixture, recomputed here from the same mesh."""
    x, y = vertices[:, 0], vertices[:, 1]
    fore = x >= 0.5 * (x.min() + x.max())

    def clearance(theta):
        yp = x * np.sin(theta) + y * np.cos(theta)
        return yp[fore].min() - (yp[~fore].min() + deepest_m)

    theta = 0.0 if clearance(0.0) >= 0 else brentq(clearance, 0.0, np.pi / 3, xtol=1e-12)
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]]), float(theta)


def sample(layer, vertices, faces, rotation, translation, plane_value_m, tag):
    """One pose: the LAYER's force, and the foundation's purely geometric integral at the same pose."""
    began = time.perf_counter()
    try:
        solved = layer.solve(rotation=rotation, translation=translation, plane_axis=PLANE_AXIS,
                             plane_value_m=plane_value_m, plane_sign=PLANE_SIGN, method='fast')
    except (ValueError, RuntimeError) as error:
        return {'tag': tag, 'not_run': str(error)[:200]}
    wall = time.perf_counter() - began
    geometry = cl.foundation_geometry(vertices, faces, rotation=rotation, translation=translation,
                                      plane_axis=PLANE_AXIS, plane_value_m=plane_value_m,
                                      plane_sign=PLANE_SIGN)
    force = np.asarray(solved['segment_force_n'], float)
    return {'tag': tag,
            'rotation': np.asarray(rotation, float).tolist(),
            'translation_m': np.asarray(translation, float).tolist(),
            'plane_value_m': float(plane_value_m),
            'layer_force_n': [float(v) for v in force],
            'layer_force_magnitude_n': float(np.linalg.norm(force)),
            'layer_moment_nm': [float(v) for v in np.asarray(solved['segment_moment_nm'], float)],
            'layer_contact_nodes': int(solved['contact_nodes']),
            'layer_iterations': int(solved['iterations']),
            'layer_converged': bool(solved['converged']),
            'layer_balance_relative': float(solved['balance_force_relative']),
            'layer_wall_s': float(wall),
            'integral_m3': float(geometry['integral_m3']),
            'foundation_faces': int(geometry['contact_faces']),
            'foundation_area_m2': float(geometry['contact_area_m2']),
            'foundation_max_depth_m': float(geometry['maximum_depth_m']),
            'footprint_major_m': float(geometry['footprint_major_m']),
            'footprint_minor_m': float(geometry['footprint_minor_m']),
            'unit_force_n': [float(v) for v in geometry['unit_force_n']],
            'unit_moment_nm': [float(v) for v in geometry['unit_moment_nm']]}


def usable(rows):
    """Poses that carry information about a stiffness: both sides strictly positive."""
    keep, dropped = [], []
    for row in rows:
        if 'not_run' in row:
            dropped.append({**row, 'why': 'the layer did not solve'})
        elif row['layer_force_magnitude_n'] <= 0:
            dropped.append({**row, 'why': 'the layer carries exactly zero here'})
        elif row['integral_m3'] <= 0:
            dropped.append({**row, 'why': 'the foundation touches no face here'})
        else:
            keep.append(row)
    return keep, dropped


def fixture_calcn(layer, vertices, faces, rotation, depths, low, tag):
    rows = []
    for depth in depths:
        row = sample(layer, vertices, faces, rotation, np.zeros(3), low + depth, tag)
        row['penetration_m'] = float(depth)
        rows.append(row)
        if 'not_run' in row:
            print('      %5.2f mm  NOT RUN: %s' % (depth * 1e3, row['not_run'][:110]), flush=True)
        else:
            print('      %5.2f mm  layer %10.4f N   faces %4d   G %.6e m3   k_eff %.4e Pa/m   %.2f s'
                  % (depth * 1e3, row['layer_force_magnitude_n'], row['foundation_faces'],
                     row['integral_m3'],
                     row['layer_force_magnitude_n'] / row['integral_m3'] if row['integral_m3'] > 0
                     else float('nan'), row['layer_wall_s']), flush=True)
    return rows


def harvest_ulna_poses():
    """The plant's own `ulna_l` transforms over an uncoupled collapse.  One engine session."""
    from ihm.assembly.articulated import ArticulatedBodyPlant
    work = WORK / 'ulna_poses'
    if work.exists():
        shutil.rmtree(work)
    plant = ArticulatedBodyPlant(ROOT, work / 'plant', environment='upright',
                                 initial_pose=POSE, mechanical_fidelity=None)
    poses = []
    try:
        for step in range(STEPS):
            native = plant.native.snapshot()
            poses.append({'step': step, 'time_s': native['time_s'],
                          'transform_ground': np.asarray(
                              native['bodies']['ulna_l']['transform_ground'], float).tolist()})
            plant.advance(DT)
    finally:
        plant.close()
    return poses


def fixture_ulna(layer, vertices, faces, poses):
    rows = []
    for pose in poses:
        T = np.asarray(pose['transform_ground'], float)
        world = vertices @ T[:3, :3].T + T[:3, 3]
        gap = float(world[:, PLANE_AXIS].min()) - FLOOR_M
        if gap >= 0:
            continue
        row = sample(layer, vertices, faces, T[:3, :3], T[:3, 3], FLOOR_M, 'U')
        row['step'] = pose['step']
        row['penetration_m'] = -gap
        rows.append(row)
        if 'not_run' in row:
            print('      step %3d  %6.3f mm  NOT RUN: %s' % (pose['step'], -gap * 1e3,
                                                             row['not_run'][:100]), flush=True)
            break                       # the layer cannot carry anything deeper on this trajectory
        print('      step %3d  %6.3f mm  layer %10.4f N   faces %4d   G %.6e m3   k_eff %.4e Pa/m   %.2f s'
              % (pose['step'], -gap * 1e3, row['layer_force_magnitude_n'], row['foundation_faces'],
                 row['integral_m3'],
                 row['layer_force_magnitude_n'] / row['integral_m3'] if row['integral_m3'] > 0
                 else float('nan'), row['layer_wall_s']), flush=True)
    return rows


def score(rows, stiffness):
    predicted = np.array([stiffness * r['integral_m3'] for r in rows])
    measured = np.array([r['layer_force_magnitude_n'] for r in rows])
    return cl.summarise(predicted, measured)


def fit_report(name, train, test, declared_k, uniform_k, poisson, thickness):
    """Fit on `train`, score on `test`, against every baseline.  Never scored on its own points."""
    g = [r['integral_m3'] for r in train]
    f = [r['layer_force_magnitude_n'] for r in train]
    k_log = cl.fit_stiffness(g, f, objective='log')
    k_linear = cl.fit_stiffness(g, f, objective='linear')
    measured = np.array([r['layer_force_magnitude_n'] for r in test])
    out = {'fitted_stiffness_pa_per_m': k_log,
           'fitted_stiffness_linear_objective_pa_per_m': k_linear,
           'objective': 'log', 'objective_basis': cl.OBJECTIVE_BASIS['log'],
           'poisson_ratio': poisson, 'layer_thickness_m': thickness,
           'implied_youngs_modulus_pa': cl.material_from_stiffness(k_log, poisson, thickness),
           'declared_stiffness_pa_per_m': declared_k,
           'declared_youngs_modulus_pa': cl.material_from_stiffness(declared_k, poisson, thickness),
           'fitted_over_declared': k_log / declared_k,
           'train_penetration_mm': [round(r['penetration_m'] * 1e3, 4) for r in train],
           'test_penetration_mm': [round(r['penetration_m'] * 1e3, 4) for r in test],
           'train_self': score(train, k_log),
           'held_out': {
               'fit': score(test, k_log),
               'fit_linear_objective': score(test, k_linear),
               'B2_declared_layer_map': score(test, declared_k),
               'B3_uniform_skin_material': score(test, uniform_k)},
           'held_out_skill_vs_B2': cl.skill([k_log * r['integral_m3'] for r in test], measured,
                                            [declared_k * r['integral_m3'] for r in test]),
           'held_out_skill_vs_B3': cl.skill([k_log * r['integral_m3'] for r in test], measured,
                                            [uniform_k * r['integral_m3'] for r in test])}
    print('\n  %s' % name, flush=True)
    print('    fitted k        %.6e Pa/m  (E = %.2f Pa at the declared v %.2f, h %.4f m)'
          % (k_log, out['implied_youngs_modulus_pa'], poisson, thickness), flush=True)
    print('    linear-objective k %.6e Pa/m  (sensitivity, not the headline)' % k_linear, flush=True)
    print('    declared k      %.6e Pa/m  -> fitted/declared = %.4f' % (declared_k, k_log / declared_k),
          flush=True)
    for key, label in (('fit', 'FIT             '), ('fit_linear_objective', 'fit (linear obj)'),
                       ('B2_declared_layer_map', 'B2 declared k   '),
                       ('B3_uniform_skin_material', 'B3 uniform E/v/h')):
        s = out['held_out'][key]
        print('    held out  %s  median |rel| %7.2f%%   max %7.2f%%   max |dF| %8.4f N   '
              'skill vs zero %7.4f   vs mean %7.4f'
              % (label, 100 * s['median_absolute_relative'], 100 * s['max_absolute_relative'],
                 s['max_absolute_n'], s['skill_vs_zero'], s['skill_vs_mean']), flush=True)
    return out


def k_eff_spread(rows):
    k = np.array([r['layer_force_magnitude_n'] / r['integral_m3'] for r in rows])
    return {'per_pose_pa_per_m': [float(v) for v in k],
            'penetration_mm': [round(r['penetration_m'] * 1e3, 4) for r in rows],
            'min_pa_per_m': float(k.min()), 'max_pa_per_m': float(k.max()),
            'geometric_mean_pa_per_m': float(np.exp(np.mean(np.log(k)))),
            'max_over_min': float(k.max() / k.min()),
            'relative_spread_about_geometric_mean':
                float(np.abs(k / np.exp(np.mean(np.log(k))) - 1).max())}


def write_bundle(fits, refusals, source_manifest, report):
    """A selectable ihm.segment-contact-meshes.v1 bundle carrying the fitted stiffness."""
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / 'meshes').mkdir(parents=True)
    manifest = json.loads(json.dumps(source_manifest))     # a copy, never the source object
    for record in manifest['records']:
        source = ROOT / stl.BUNDLE / 'meshes' / record['mesh_file']
        data = source.read_bytes()
        if hashlib.sha256(data).hexdigest() != record['written_sha256']:
            raise ValueError('Source skin mesh does not match its record: ' + record['mesh_file'])
        (OUT / 'meshes' / record['mesh_file']).write_bytes(data)
        fit = fits.get(record['body'])
        if fit is None:
            refused = refusals.get(record['body'])
            record['layer']['contact_law_fit'] = {
                'fitted': False,
                'stiffness_source': 'layer map (k = E_app/h), UNFITTED',
                'basis': 'No deformable-layer measurement exists for this segment, so its '
                         'stiffness is the layer map\'s own and this bundle changes nothing '
                         'about it.' if refused is None else
                         'A layer measurement EXISTS for this segment and its fit FAILED its '
                         'own bar (F1), so by this measurement\'s pre-registered outcome (c) '
                         'the deliverable is the measurement and not a law. The layer map\'s '
                         'own k is kept and is NOT endorsed by any fit.',
                **({} if refused is None else {'fit_refused': refused})}
            continue
        record['layer']['layer_map_stiffness_pa_per_m'] = record['layer']['stiffness_pa_per_m']
        record['layer']['stiffness_pa_per_m'] = fit['fitted_stiffness_pa_per_m']
        record['layer']['contact_law_fit'] = {
            'fitted': True,
            'stiffness_source': 'fitted to ihm/assembly/soft_tissue_layer.py at ' + IDENTITY,
            'fitted_over_layer_map': fit['fitted_over_declared'],
            'implied_youngs_modulus_pa': fit['implied_youngs_modulus_pa'],
            'held_out_median_absolute_relative': fit['held_out']['fit']['median_absolute_relative'],
            'held_out_max_absolute_relative': fit['held_out']['fit']['max_absolute_relative'],
            'fixture': fit['fixture'],
            'basis': cl.STIFFNESS_BASIS}
    manifest['layer_map'] = {**manifest.get('layer_map', {}),
                             'rule': 'per-segment k FITTED to the deformable soft-tissue layer '
                                     'where a layer measurement exists; the layer map\'s own '
                                     'k = E_app/h elsewhere, marked per record',
                             'fitted_bodies': sorted(fits)}
    manifest['contact_law_fit'] = {
        'schema': 'ihm.contact-law-fit.v1',
        'source_bundle': stl.BUNDLE,
        'source_manifest_sha256': hashlib.sha256(
            (ROOT / stl.BUNDLE / 'manifest.json').read_bytes()).hexdigest(),
        'layer_identity': IDENTITY,
        'fitted_bodies': sorted(fits),
        'unfitted_bodies': sorted(r['body'] for r in manifest['records'] if r['body'] not in fits),
        'refused_bodies': sorted(refusals),
        'refusal_basis': 'A segment whose fit failed F1 carries NO fitted stiffness: the '
                         'pre-registered outcome (c) makes the measurement the deliverable.',
        'foundation_basis': cl.FOUNDATION_BASIS,
        'rate_basis': cl.RATE_BASIS,
        'measurement': 'scripts/measure_contact_law_fit.py',
        'battery': 'scripts/verify_contact_law_fit.py',
        'report': 'data/derived/contact-law-fit-v1/report.json',
        'bars': {'F1_held_out_median_relative': BAR_F1, 'F2_margin_over_measured_baseline': BAR_F2},
        'outcome': report['outcome'],
        'disclosure': cl.DISCLOSURE}
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=1, sort_keys=True) + '\n')
    (OUT / 'report.json').write_text(json.dumps(report, indent=1, sort_keys=True,
                                                default=lambda v: getattr(v, 'tolist', lambda: str(v))()) + '\n')
    return manifest


def refinement_arm(vertices, faces, poses, steps, spacings=(0.005, 0.004, 0.003)):
    """RECORDED (A6).  Is a k_eff outlier the LAYER's discretisation or the contact law?

    The poses where run 1's k_eff was an outlier are exactly the poses where the layer's own
    contact patch is carried by a handful of nodes.  That is a HYPOTHESIS about a mechanism,
    and a factorial does not test a mechanism (CLAUDE.md), so it gets its own instrument:
    rebuild the layer finer and re-solve THE SAME poses.  No gate, no bar.
    """
    print('\nA6  REFINEMENT ARM (recorded): is the k_eff outlier the LAYER\'s discretisation?',
          flush=True)
    by_step = {p['step']: np.asarray(p['transform_ground'], float) for p in poses}
    out = []
    for spacing in spacings:
        try:
            layer = stl.segment_layer(ROOT, 'ulna_l', surface='fitted', depth='local',
                                      spacing_m=spacing)
        except (MemoryError, ValueError, RuntimeError) as error:
            print('      spacing %.1f mm: NOT BUILT: %s' % (spacing * 1e3, str(error)[:110]),
                  flush=True)
            out.append({'spacing_m': spacing, 'not_run': str(error)[:200]})
            continue
        row = {'spacing_m': spacing, 'dof': int(layer.dof), 'nodes': int(len(layer.local)),
               'poses': {}}
        for step in steps:
            T = by_step[step]
            try:
                solved = layer.solve(rotation=T[:3, :3], translation=T[:3, 3],
                                     plane_axis=PLANE_AXIS, plane_value_m=FLOOR_M,
                                     plane_sign=PLANE_SIGN, method='fast')
            except (ValueError, RuntimeError) as error:
                print('      spacing %.1f mm  step %3d: NOT RUN: %s'
                      % (spacing * 1e3, step, str(error)[:100]), flush=True)
                row['poses'][str(step)] = {'not_run': str(error)[:200]}
                continue
            geometry = cl.foundation_geometry(vertices, faces, rotation=T[:3, :3],
                                              translation=T[:3, 3], plane_axis=PLANE_AXIS,
                                              plane_value_m=FLOOR_M, plane_sign=PLANE_SIGN)
            force = float(np.linalg.norm(solved['segment_force_n']))
            k_eff = force / geometry['integral_m3'] if geometry['integral_m3'] > 0 else float('nan')
            row['poses'][str(step)] = {'force_n': force, 'k_eff_pa_per_m': k_eff,
                                       'layer_contact_nodes': int(solved['contact_nodes'])}
            print('      spacing %.1f mm (%6d DOF)  step %3d  layer %9.4f N  nodes %4d  '
                  'k_eff %.4e Pa/m' % (spacing * 1e3, layer.dof, step, force,
                                       solved['contact_nodes'], k_eff), flush=True)
        out.append(row)
        del layer
    return out


def heel_fixture_rows(body, layer, vertices, faces, deepest_m, tag):
    rotation, theta = heel_rotation(vertices, deepest_m)
    low = float((vertices @ rotation.T)[:, PLANE_AXIS].min())
    print('\n%s  %s heel fixture at %.0f mm forefoot clearance: %.10f deg about z'
          % (tag, body, deepest_m * 1e3, np.degrees(theta)), flush=True)
    rows = fixture_calcn(layer, vertices, faces, rotation, DEPTHS_HEEL, low, tag)
    return rows, rotation, low, float(theta)


def main():
    began = time.perf_counter()
    print('THE CONTACT LAW THE LAYER IMPLIES: fit, held-out score, and the substitution cost',
          flush=True)
    print('layer identity %s; support: OpenSim upright floor, y = %.1f\n' % (IDENTITY, FLOOR_M),
          flush=True)

    manifest, calcn_record, calcn_v, calcn_f = bundle_mesh('calcn_l')
    _, calcn_r_record, calcn_r_v, calcn_r_f = bundle_mesh('calcn_r')
    _, ulna_record, ulna_v, ulna_f = bundle_mesh('ulna_l')
    poisson = float(manifest['skin_material']['poissons_ratio'])
    uniform_k = cl.stiffness_from_material(manifest['skin_material']['youngs_modulus_pa'],
                                           poisson, manifest['skin_material']['layer_thickness_m'])
    print('B3  the uniform skin_material triple: E %.1f Pa, v %.2f, h %.4f m  ->  k %.6e Pa/m'
          % (manifest['skin_material']['youngs_modulus_pa'], poisson,
             manifest['skin_material']['layer_thickness_m'], uniform_k), flush=True)

    report = {'schema': 'ihm.contact-law-fit-report.v1', 'layer_identity': IDENTITY,
              'bundle': stl.BUNDLE, 'uniform_stiffness_pa_per_m': uniform_k,
              'bars': {'F1': BAR_F1, 'F2': BAR_F2}, 'fixtures': {}}

    # --- U first: it needs the one engine session, and it closes before any layer is built --
    print('\nU  the plant\'s own poses: an uncoupled 140-step collapse from pelvis_ty = 1.03',
          flush=True)
    poses = harvest_ulna_poses()
    report['ulna_poses_steps'] = len(poses)

    # --- calcn_l: the heel fixture (H), the seam-wedge fixture (S), the 3 mm heel (H2) -------
    calcn = stl.segment_layer(ROOT, 'calcn_l', surface='fitted', depth='local')
    heel_rows, rotation, low_heel, theta = heel_fixture_rows('calcn_l', calcn, calcn_v, calcn_f,
                                                             0.006, 'H')
    report['heel_fixture_deg'] = {'calcn_l_6mm': np.degrees(theta)}
    low_flat = float(calcn_v[:, PLANE_AXIS].min())
    print('\nS  calcn_l at identity rotation -- the SEAM WEDGE fixture (amendment A2), recorded',
          flush=True)
    flat_rows = fixture_calcn(calcn, calcn_v, calcn_f, np.eye(3), DEPTHS_FLAT, low_flat, 'S')

    # CALL IT TWICE: the instrument must be a function of its arguments (CLAUDE.md).
    again = sample(calcn, calcn_v, calcn_f, rotation, np.zeros(3), low_heel + DEPTHS_HEEL[3], 'H')
    twice = (again['layer_force_n'] == heel_rows[3]['layer_force_n']
             and again['integral_m3'] == heel_rows[3]['integral_m3'])
    gate('T1', 'call it twice at the same pose: layer force and foundation integral bitwise', twice)
    del again

    h2_rows, _, _, theta2 = heel_fixture_rows('calcn_l', calcn, calcn_v, calcn_f, 0.003, 'H2')
    report['heel_fixture_deg']['calcn_l_3mm'] = np.degrees(theta2)
    del calcn

    # --- calcn_r: the replication (A4) --------------------------------------------------------
    calcn_r = stl.segment_layer(ROOT, 'calcn_r', surface='fitted', depth='local')
    hr_rows, _, _, theta_r = heel_fixture_rows('calcn_r', calcn_r, calcn_r_v, calcn_r_f,
                                               0.006, 'HR')
    report['heel_fixture_deg']['calcn_r_6mm'] = np.degrees(theta_r)
    del calcn_r

    ulna = stl.segment_layer(ROOT, 'ulna_l', surface='fitted', depth='local')
    print('\nU  ulna_l at those poses, solved COLD at each (never warm-started)', flush=True)
    ulna_rows = fixture_ulna(ulna, ulna_v, ulna_f, poses)
    del ulna

    fixtures = {'H': ('calcn_l heel fixture, 6 mm forefoot clearance', heel_rows, calcn_record),
                'H2': ('calcn_l heel fixture, 3 mm forefoot clearance (cross-fixture only)',
                       h2_rows, calcn_record),
                'HR': ('calcn_r heel fixture, 6 mm forefoot clearance', hr_rows, calcn_r_record),
                'S': ('calcn_l identity rotation -- the SEAM WEDGE, recorded not judged',
                      flat_rows, calcn_record),
                'U': ('ulna_l, the plant\'s own collapse poses', ulna_rows, ulna_record)}
    kept = {}
    for key in ('H', 'H2', 'HR', 'S', 'U'):
        label, rows, record = fixtures[key]
        keep, dropped = usable(rows)
        kept[key] = keep
        report['fixtures'][key] = {'label': label, 'rows': rows,
                                   'dropped': [{'penetration_mm': d.get('penetration_m', 0) * 1e3,
                                                'why': d['why'],
                                                'detail': d.get('not_run')} for d in dropped],
                                   'usable': len(keep)}
        for d in dropped:
            print('      DROPPED  %s at %.3f mm: %s' % (key, d.get('penetration_m', 0) * 1e3, d['why']),
                  flush=True)
        if keep:
            report['fixtures'][key]['k_eff'] = k_eff_spread(keep)

    # --- the fits -----------------------------------------------------------------------------
    print('\n' + '=' * 96, flush=True)
    print('THE FIT: one stiffness per segment, trained on one set of poses, scored on another',
          flush=True)
    print('=' * 96, flush=True)
    fits = {}
    for key, body in (('H', 'calcn_l'), ('HR', 'calcn_r'), ('U', 'ulna_l')):
        rows = kept[key]
        if len(rows) < 4:
            gate('F1' + key, 'fixture %s has enough usable poses to fit and score' % key, False,
                 '%d usable' % len(rows))
            continue
        train, test = rows[0::2], rows[1::2]
        record = fixtures[key][2]
        out = fit_report('%s  (%s)' % (key, fixtures[key][0]), train, test,
                         record['layer']['stiffness_pa_per_m'], uniform_k,
                         poisson, record['layer']['thickness_m'])
        out['fixture'] = key
        out['body'] = body
        report['fixtures'][key]['fit'] = out
        fits[body] = out

    # --- CROSS-FIXTURE: the heel fit on pose families it never saw ----------------------------
    report['cross_fixture'] = {}
    if 'calcn_l' in fits:
        k_fit = fits['calcn_l']['fitted_stiffness_pa_per_m']
        k_declared = calcn_record['layer']['stiffness_pa_per_m']
        for key, note in (('H2', 'the 3 mm heel fixture -- a different pose family, still the heel'),
                          ('S', 'the SEAM WEDGE fixture -- recorded, NOT used to judge the law')):
            if not kept[key]:
                continue
            entry = {'note': note, 'fit': score(kept[key], k_fit),
                     'B2_declared_layer_map': score(kept[key], k_declared),
                     'points': len(kept[key])}
            report['cross_fixture'][key] = entry
            print('\n  CROSS-FIXTURE H -> %s  (%s; %d usable points)' % (key, note, entry['points']),
                  flush=True)
            for which, label in (('fit', 'fit     '), ('B2_declared_layer_map', 'declared')):
                s = entry[which]
                print('    %s  median |rel| %8.2f%%   max %8.2f%%   max |dF| %8.4f N'
                      % (label, 100 * s['median_absolute_relative'],
                         100 * s['max_absolute_relative'], s['max_absolute_n']), flush=True)

    # --- the gates ----------------------------------------------------------------------------
    print('\n' + '=' * 96, flush=True)
    outcome, refusals = {}, {}
    for key, body in (('H', 'calcn_l'), ('HR', 'calcn_r'), ('U', 'ulna_l')):
        fit = fits.get(body)
        if fit is None:
            continue
        held = fit['held_out']['fit']['median_absolute_relative']
        declared = fit['held_out']['B2_declared_layer_map']['median_absolute_relative']
        f1 = gate('F1' + key, '%s: held-out median relative error <= %.2f' % (body, BAR_F1),
                  held <= BAR_F1, '%.4f' % held)
        f2 = gate('F2' + key, '%s: beats the UNFITTED layer-map k by more than sqrt(2)*bar'
                  % body, declared - held > BAR_F2,
                  'fit %.4f vs declared %.4f, margin %.4f, needs > %.4f'
                  % (held, declared, declared - held, BAR_F2))
        outcome[body] = {'F1': f1, 'F2': f2,
                         'outcome': 'a' if (f1 and f2) else ('b' if f1 else 'c'),
                         'fixture': key,
                         'fitted_stiffness_pa_per_m': fit['fitted_stiffness_pa_per_m'],
                         'declared_stiffness_pa_per_m': fit['declared_stiffness_pa_per_m'],
                         'held_out_median_relative': held,
                         'declared_median_relative': declared}
        if not f1:
            # PRE-REGISTERED OUTCOME (c): the deliverable is the measurement, not a law.
            refusals[body] = {'reason': 'F1 FAILED: the fit does not reproduce the layer on '
                                        'held-out poses within the layer\'s own 10% bar',
                              'held_out_median_relative': held,
                              'k_eff_max_over_min': report['fixtures'][key]['k_eff']['max_over_min'],
                              'would_have_been_pa_per_m': fit['fitted_stiffness_pa_per_m']}
            fits.pop(body)
        gate('F3' + key, '%s: k_eff spread over the fixture (what a single k CANNOT follow)' % body,
             True, 'max/min %.4f, worst %.2f%% from the geometric mean'
             % (report['fixtures'][key]['k_eff']['max_over_min'],
                100 * report['fixtures'][key]['k_eff']['relative_spread_about_geometric_mean']),
             recorded=True)
    gate('F4', 'the layer is quasi-static: no rate term is fitted, dissipation is untouched',
         True, cl.RATE_BASIS[:80], recorded=True)

    # --- A5: the footprint the confined reading assumes ---------------------------------------
    print('\nA5  the contact FOOTPRINT against the declared layer thickness (recorded)', flush=True)
    footprints = {}
    for key, body in (('H', 'calcn_l'), ('HR', 'calcn_r'), ('U', 'ulna_l')):
        record = fixtures[key][2]
        thickness = record['layer']['thickness_m']
        rows = [{'penetration_mm': r['penetration_m'] * 1e3,
                 'major_mm': r['footprint_major_m'] * 1e3,
                 'minor_mm': r['footprint_minor_m'] * 1e3,
                 'minor_over_thickness': r['footprint_minor_m'] / thickness} for r in kept[key]]
        footprints[body] = {'declared_thickness_m': thickness, 'rows': rows}
        if rows:
            print('      %-8s h %5.2f mm   footprint minor %5.1f - %5.1f mm   minor/h %.2f - %.2f'
                  % (body, thickness * 1e3, min(r['minor_mm'] for r in rows),
                     max(r['minor_mm'] for r in rows),
                     min(r['minor_over_thickness'] for r in rows),
                     max(r['minor_over_thickness'] for r in rows)), flush=True)
    report['footprints'] = footprints
    gate('A5', 'contact footprint minor extent over the declared thickness, per pose', True,
         'recorded for %d segments' % len(footprints), recorded=True)

    # --- A6: the refinement arm ----------------------------------------------------------------
    report['refinement'] = refinement_arm(ulna_v, ulna_f, poses, (81, 98, 103))
    gate('A6', 'refinement arm on the poses whose layer contact patch is a handful of nodes',
         True, 'recorded at %d spacings' % len(report['refinement']), recorded=True)

    report['outcome'] = outcome
    report['refused'] = refusals
    report['results'] = RESULTS
    report['wall_seconds'] = time.perf_counter() - began
    report['peak_rss_mb'] = peak_mb()

    write_bundle(fits, refusals, manifest, report)
    print('\nwrote %s  (%d fitted, %d refused, of %d records)'
          % (OUT, len(fits), len(refusals), len(manifest['records'])), flush=True)
    print('wall %.1f s, peak RSS %.0f MB' % (report['wall_seconds'], report['peak_rss_mb']),
          flush=True)
    failed = [r['gate'] for r in RESULTS if not r['pass'] and not r['recorded']]
    print('\n%d gates, %d FAILED%s' % (len(RESULTS), len(failed),
                                       '' if not failed else ': ' + ', '.join(failed)), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
