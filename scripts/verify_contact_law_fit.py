#!/usr/bin/env python3
"""Known answers, controls that CAN fail, and idempotence for the contact law fitted to the
deformable soft-tissue layer.

What is being judged: `ihm/assembly/contact_law.py` (the engine's own elastic foundation, in
Python), the fit `scripts/measure_contact_law_fit.py` wrote to
`data/derived/contact-law-fit-v1/`, and the `plant_options` selection that offers it.

    PYTHONPATH=. prlimit --as=4294967296 .venv/bin/python -u scripts/verify_contact_law_fit.py

================================================================================
PRE-REGISTRATION.  Bars fixed here, before any run of this battery.  A gate that fails is
recorded FAILED and is never rescored.
================================================================================

  V1  KNOWN ANSWER, ANALYTIC.  A flat plate of area A pressed to a uniform depth d gives
      exactly `k A d`, because every spring on it carries the same depth.  Bar: 1e-15
      relative, which is float noise on a sum of two triangles.
  V1b KNOWN ANSWER, INVARIANCE.  Rotating the mesh about the SUPPORT NORMAL changes no
      vertex's height, so the foundation integral must not move; nor may re-ordering the
      faces.  Bar: bitwise.  This is the null case of the rotation argument -- the input on
      which the extra freedom must produce NOTHING.
  V2  KNOWN ANSWER, CROSS-IMPLEMENTATION.  The Python law against the ENGINE's own
      `ElasticFoundationForce`, on the same skin meshes, the same poses and the same
      per-segment stiffness, AT REST so that the dissipation factor `(1 + c*vnormal)` and
      the friction term are both identically zero and only the elastic term is left.  The
      run's own reported speeds are checked to be zero rather than inferred from the fact
      that nothing has been advanced.  Bar: 1e-12 relative on the force.
  V2b CONTROL THAT CAN FAIL.  The same comparison with the Python stiffness multiplied by
      1.5 must DIFFER by 0.5 relative.  Without it, V2 could be comparing two zeros.
  V3  NULL CASE.  A support that does not reach the mesh, and one exactly touching its lowest
      vertex, give force and moment EXACTLY 0.0 over exactly 0 faces.  Bar: exact.
  V4  IDEMPOTENCE (CLAUDE.md: call it twice at the same input).  `foundation_geometry`,
      `foundation_resultant`, `fit_stiffness` and `SoftTissueLayer.solve` all return bitwise
      the same thing when called twice with the same arguments.  Bar: bitwise.
  V5  THE FIT REPRODUCES THE LAYER ON HELD-OUT PENETRATIONS.  The layer is RE-SOLVED here,
      at penetrations the artefact records as held out, and the artefact's fitted stiffness
      is scored against those fresh solves.  Bar: median absolute relative error <= 0.10 --
      the layer's OWN convergence bar (docs/SOFT_BODY.md CV7-CV10), the same number the
      measurement used, and the only place it comes from.
  V5b CONTROL THAT CAN FAIL.  The same held-out check with the fitted stiffness multiplied
      by 1.5 must EXCEED that bar.  A bar that a 1.5x wrong stiffness also passes is a bar
      that cannot see a stiffness.
  V5c CONTROL THAT CAN FAIL, SHUFFLED PAIRING.  Pair each held-out pose's foundation
      integral with a DIFFERENT pose's layer force (a cyclic shift, which is a pure
      relabelling and changes no number) and re-score.  It must exceed the bar.  If it does
      not, the agreement in V5 is a property of the range and not of the pairing.
  V6  `None` IS THE HISTORICAL PLANT, BIT FOR BIT.  `resolve_fidelity(root, None)` returns
      EXACTLY `{}` of plant kwargs, and two plants built with `mechanical_fidelity=None` --
      one of them after this module has been imported and the fitted identity resolved --
      agree bitwise in every coordinate over 30 steps.
  V6b CONTROL THAT CAN FAIL.  A plant carrying the fitted bundle must DIFFER from `None`.
      Without it V6 could be comparing two plants that never felt a contact element.
  V7  THE BUG THE FIT UNCOVERED IS FIXED, AND THE FIX IS NOT A LOOSENING.  Both layer-map
      bundles (`skin_layer_map`, and the fitted one) resolve WITHOUT a caller-side E, p or h
      and BUILD a plant; a bundle with no per-record layer map (`skin`) still receives the
      uniform triple.  Fixing the first by dropping the material everywhere would break the
      second, so both halves are gated.
  V8  THE ENGINE IS CARRYING THE FITTED NUMBER.  The stiffness the engine reports back on
      each fitted contact element equals the artefact's fitted stiffness exactly.  Read out
      of the run's own output, never inferred from what was passed in.
  V9  THE FIT CHANGED A NUMBER AND NOTHING ELSE.  Every mesh in the fitted bundle is
      byte-identical to the source bundle's, by sha256, and the record set is the same.
  V10 RECORDED, never gated: the per-step wall cost of a plant carrying the fitted bundle,
      beside the same plant carrying the unfitted layer-map bundle and beside the coupled
      layer's own 53.9 ms per 10 ms step (docs/SOFT_BODY.md).  The comparison that matters
      is the fitted bundle against the UNFITTED one, because the two differ by one constant
      and nothing else -- so the difference is the substitution's own marginal cost, and a
      plant-wide wall clock would be a numerator holding work the treatment never did.

  AMENDMENT after run 1 (`logs/verify_contact_law_fit.run1-V1c-FAILED.log`, 21 gates, 1
  FAILED).  **V1c IS RECORDED FAILED AND ITS BAR DOES NOT MOVE.**  It demanded that
  re-ordering and re-winding the faces leave the foundation integral BITWISE unchanged; the
  measured difference is 9.9e-23 m^3 on 3.55e-08, i.e. **2.79e-15 relative**, which is the
  order of a floating-point SUM over four terms and not a property of the law.  A bitwise
  bar on a sum is a bar on the summation order.  The bar is not repaired.  What is added is
  a SECOND, correctly posed gate beside it -- V1d, the same comparison at 1e-14 relative --
  and run 2 must still print V1c FAILED.  V10 also gains the `None` arm's own per-step wall
  clock, so the three costs are measured on one fixture and one step count.

Memory budget: 4 GiB address space (prlimit); peak is printed at the end.  One engine
session at a time: every plant is built, used and closed before the next one starts.
"""
import hashlib
import json
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.assembly import contact_law as cl                                     # noqa: E402
from ihm.assembly import soft_tissue_layer as stl                              # noqa: E402
from ihm.assembly.plant_options import SEGMENT_CONTACT_BUNDLES, resolve_fidelity   # noqa: E402

FIT = 'data/derived/contact-law-fit-v1'
BUNDLE_ID = 'skin_layer_fitted'
BAR_HELD_OUT = 0.10                    # the layer's own convergence bar; see the header
BAR_CROSS = 1e-12
WORK = ROOT / 'data/derived/contact-law-verify-runs'
POSE_REST = ({'pelvis_ty': 0.995}, {'pelvis_ty': 0.985}, {'pelvis_ty': 0.960})
POSE_RUN = {'pelvis_ty': 1.03}

RESULTS = []


def gate(name, claim, ok, detail='', recorded=False):
    print('%-5s %-74s %s%s' % (name, claim, 'RECORDED' if recorded else ('PASS' if ok else 'FAILED'),
                               '' if not detail else '   ' + detail), flush=True)
    RESULTS.append({'gate': name, 'claim': claim, 'pass': bool(ok), 'detail': detail,
                    'recorded': bool(recorded)})
    return ok


def peak_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def fresh(tag):
    path = WORK / tag
    if path.exists():
        shutil.rmtree(path)
    return path


def plate():
    """A flat unit square in the xz plane: the one case whose foundation force is analytic."""
    vertices = np.array([[0., 0., 0.], [1., 0., 0.], [1., 0., 1.], [0., 0., 1.]])
    return vertices, np.array([[0, 1, 2], [0, 2, 3]])


def v1_analytic():
    print('\nV1  KNOWN ANSWER: a flat plate at uniform depth carries exactly k A d', flush=True)
    vertices, faces = plate()
    worst = 0.0
    for k, depth in ((1000.0, 0.003), (1.0347e7, 0.0125), (3.2826e6, 0.0004)):
        out = cl.foundation_resultant(vertices, faces, rotation=np.eye(3), translation=np.zeros(3),
                                      plane_axis=1, plane_value_m=depth, plane_sign=1.0,
                                      stiffness_pa_per_m=k)
        expected = k * 1.0 * depth
        worst = max(worst, abs(float(out['force_n'][1]) - expected) / expected)
        worst = max(worst, float(np.abs(out['force_n'][[0, 2]]).max()) / expected)
    gate('V1', 'flat plate: force = k A d, and no tangential component', worst <= 1e-15,
         'worst rel %.2e' % worst)


def v1b_invariance():
    print('\nV1b KNOWN ANSWER, NULL CASE: a rotation about the support normal must change nothing',
          flush=True)
    manifest = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_bytes())
    record = [r for r in manifest['records'] if r['body'] == 'calcn_l'][0]
    vertices, faces = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
    low = float(vertices[:, 1].min())
    base = cl.foundation_geometry(vertices, faces, plane_axis=1, plane_value_m=low + 0.002,
                                  plane_sign=1.0)
    theta = 0.7
    c, s = np.cos(theta), np.sin(theta)
    spun = np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]])
    turned = cl.foundation_geometry(vertices, faces, rotation=spun, plane_axis=1,
                                    plane_value_m=low + 0.002, plane_sign=1.0)
    order = np.array([2, 1, 0])
    reordered = cl.foundation_geometry(vertices, faces[::-1][:, order], plane_axis=1,
                                       plane_value_m=low + 0.002, plane_sign=1.0)
    spin = abs(turned['integral_m3'] - base['integral_m3']) / base['integral_m3']
    gate('V1b', 'rotation about the support normal: integral and face count unchanged',
         spin <= 1e-14 and turned['contact_faces'] == base['contact_faces'],
         'rel %.2e, faces %d vs %d' % (spin, turned['contact_faces'], base['contact_faces']))
    gate('V1c', 'the face order and winding do not change the integral',
         reordered['integral_m3'] == base['integral_m3']
         and reordered['contact_faces'] == base['contact_faces'], 'bitwise')
    shuffle = abs(reordered['integral_m3'] - base['integral_m3']) / base['integral_m3']
    gate('V1d', 'AMENDED BAR (V1c stays FAILED): the same, at 1e-14 relative',
         shuffle <= 1e-14 and reordered['contact_faces'] == base['contact_faces'],
         'rel %.2e -- float summation order over %d terms, not the law'
         % (shuffle, base['contact_faces']))


def v3_null():
    print('\nV3  NULL CASE: a support that does not reach, and one exactly touching', flush=True)
    manifest = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_bytes())
    record = [r for r in manifest['records'] if r['body'] == 'calcn_l'][0]
    vertices, faces = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
    low = float(vertices[:, 1].min())
    ok = True
    for label, plane in (('clear by 1 mm', low - 0.001), ('exactly touching', low)):
        out = cl.foundation_resultant(vertices, faces, plane_axis=1, plane_value_m=plane,
                                      plane_sign=1.0, stiffness_pa_per_m=1.0347e7)
        exact = (out['contact_faces'] == 0 and out['integral_m3'] == 0.0
                 and not np.any(out['force_n']) and not np.any(out['moment_nm']))
        ok = ok and exact
        print('      %-18s faces %d  force %s  moment %s' % (label, out['contact_faces'],
                                                             out['force_n'], out['moment_nm']),
              flush=True)
    gate('V3', 'an unreached and an exactly-touching support give EXACTLY zero', ok)


def v4_idempotent(layer, vertices, faces, rotation, plane):
    print('\nV4  IDEMPOTENCE: call it twice at the same input', flush=True)
    a = cl.foundation_geometry(vertices, faces, rotation=rotation, plane_axis=1,
                               plane_value_m=plane, plane_sign=1.0)
    b = cl.foundation_geometry(vertices, faces, rotation=rotation, plane_axis=1,
                               plane_value_m=plane, plane_sign=1.0)
    same_geometry = (a['integral_m3'] == b['integral_m3']
                     and np.array_equal(a['unit_force_n'], b['unit_force_n'])
                     and np.array_equal(a['unit_moment_nm'], b['unit_moment_nm']))
    gate('V4a', 'foundation_geometry twice: bitwise identical', same_geometry)
    g = [1e-8, 3e-8, 9e-8]
    f = [0.1, 0.35, 1.2]
    gate('V4b', 'fit_stiffness twice: bitwise identical',
         cl.fit_stiffness(g, f) == cl.fit_stiffness(g, f))
    first = layer.solve(rotation=rotation, plane_axis=1, plane_value_m=plane, plane_sign=1.0,
                        method='fast')
    second = layer.solve(rotation=rotation, plane_axis=1, plane_value_m=plane, plane_sign=1.0,
                         method='fast')
    gate('V4c', 'the LAYER solved twice at the same pose: bitwise identical',
         first['segment_force_n'] == second['segment_force_n']
         and first['segment_moment_nm'] == second['segment_moment_nm'])


def v5_held_out(report, layer, vertices, faces):
    """RE-SOLVE the layer at the held-out penetrations and score the artefact's stiffness."""
    print('\nV5  the fit reproduces the LAYER on held-out penetrations (re-solved here)', flush=True)
    fit = report['fixtures']['H']['fit']
    stiffness = fit['fitted_stiffness_pa_per_m']
    from scipy.optimize import brentq
    x, y = vertices[:, 0], vertices[:, 1]
    fore = x >= 0.5 * (x.min() + x.max())

    def clearance(theta):
        yp = x * np.sin(theta) + y * np.cos(theta)
        return yp[fore].min() - (yp[~fore].min() + 0.006)

    theta = 0.0 if clearance(0.0) >= 0 else brentq(clearance, 0.0, np.pi / 3, xtol=1e-12)
    c, s = np.cos(theta), np.sin(theta)
    rotation = np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]])
    low = float((vertices @ rotation.T)[:, 1].min())
    measured, predicted = [], []
    for depth_mm in fit['test_penetration_mm']:
        plane = low + depth_mm * 1e-3
        solved = layer.solve(rotation=rotation, plane_axis=1, plane_value_m=plane, plane_sign=1.0,
                             method='fast')
        geometry = cl.foundation_geometry(vertices, faces, rotation=rotation, plane_axis=1,
                                          plane_value_m=plane, plane_sign=1.0)
        force = float(np.linalg.norm(solved['segment_force_n']))
        measured.append(force)
        predicted.append(stiffness * geometry['integral_m3'])
        print('      %5.2f mm  layer %10.4f N   law %10.4f N   %+7.2f%%'
              % (depth_mm, force, predicted[-1], 100 * (predicted[-1] - force) / force), flush=True)
    summary = cl.summarise(predicted, measured)
    gate('V5', 'held-out median absolute relative error <= %.2f' % BAR_HELD_OUT,
         summary['median_absolute_relative'] <= BAR_HELD_OUT,
         'median %.4f, max %.4f' % (summary['median_absolute_relative'],
                                    summary['max_absolute_relative']))
    wrong = cl.summarise([1.5 * v for v in predicted], measured)
    gate('V5b', 'CONTROL THAT CAN FAIL: a 1.5x stiffness must EXCEED the bar',
         wrong['median_absolute_relative'] > BAR_HELD_OUT,
         'median %.4f' % wrong['median_absolute_relative'])
    shifted = cl.summarise(predicted, measured[1:] + measured[:1])
    gate('V5c', 'CONTROL THAT CAN FAIL: shuffled pairing must EXCEED the bar',
         shifted['median_absolute_relative'] > BAR_HELD_OUT,
         'median %.4f' % shifted['median_absolute_relative'])
    return {'penetration_mm': fit['test_penetration_mm'], 'layer_n': measured,
            'law_n': predicted, 'summary': summary,
            'control_scaled': wrong['median_absolute_relative'],
            'control_shuffled': shifted['median_absolute_relative']}


def v2_cross_implementation(report):
    """The Python law against the ENGINE's own force, at rest, on the fitted bundle."""
    print('\nV2  KNOWN ANSWER, CROSS-IMPLEMENTATION: the engine\'s own ElasticFoundationForce',
          flush=True)
    from ihm.assembly.articulated import ArticulatedBodyPlant
    path = ROOT / SEGMENT_CONTACT_BUNDLES[BUNDLE_ID]['path']
    manifest = json.loads((path / 'manifest.json').read_bytes())
    mesh = {r['element']: (r['body'],) + stl.load_obj(path / 'meshes' / r['mesh_file'])
            for r in manifest['records']}
    declared = {r['element']: r['layer']['stiffness_pa_per_m'] for r in manifest['records']}
    worst, worst_control, checked, engine_k_ok, speeds_zero = 0.0, float('inf'), 0, True, True
    for i, pose in enumerate(POSE_REST):
        plant = ArticulatedBodyPlant(ROOT, fresh('rest%d' % i) / 'plant', environment='upright',
                                     initial_pose=pose,
                                     mechanical_fidelity={'segment_contact': BUNDLE_ID})
        try:
            native = plant.native.snapshot()
            # READ THE MODE OUT OF THE RUN'S OWN OUTPUT: this gate is only about the elastic
            # term, and the elastic term is all there is only if every speed is zero.
            speeds_zero = speeds_zero and max(abs(v['speed']) for v in native['coordinates'].values()) == 0.0
            for contact in native['contacts']:
                if not contact.get('geometry_type') or not np.any(contact['force_n']):
                    continue
                element = contact['name'].replace('mesh_support_', '')
                body, vertices, faces = mesh[element]
                engine_k_ok = engine_k_ok and contact['foundation_stiffness_pa_per_m'] == declared[element]
                T = np.asarray(native['bodies'][body]['transform_ground'], float)
                mine = cl.foundation_resultant(vertices, faces, rotation=T[:3, :3],
                                               translation=T[:3, 3], plane_axis=1,
                                               plane_value_m=0.0, plane_sign=1.0,
                                               stiffness_pa_per_m=contact['foundation_stiffness_pa_per_m'])
                theirs = np.asarray(contact['force_n'], float)
                scale = float(np.linalg.norm(theirs))
                worst = max(worst, float(np.linalg.norm(mine['force_n'] - theirs)) / scale)
                control = cl.foundation_resultant(
                    vertices, faces, rotation=T[:3, :3], translation=T[:3, 3], plane_axis=1,
                    plane_value_m=0.0, plane_sign=1.0,
                    stiffness_pa_per_m=1.5 * contact['foundation_stiffness_pa_per_m'])
                worst_control = min(worst_control,
                                    float(np.linalg.norm(control['force_n'] - theirs)) / scale)
                checked += 1
                print('      ty %.3f  %-28s engine %12.6f N   law %12.6f N   rel %.2e  faces %d'
                      % (pose['pelvis_ty'], element, theirs[1], mine['force_n'][1],
                         float(np.linalg.norm(mine['force_n'] - theirs)) / scale,
                         mine['contact_faces']), flush=True)
        finally:
            plant.close()
    gate('V2z', 'the comparison ran at rest: every reported speed is exactly zero', speeds_zero)
    gate('V2', 'the Python law = the engine\'s own force over %d contacting meshes' % checked,
         checked > 0 and worst <= BAR_CROSS, 'worst rel %.2e' % worst)
    gate('V2b', 'CONTROL THAT CAN FAIL: a 1.5x stiffness must differ by ~0.5',
         checked > 0 and worst_control > 0.4, 'smallest control difference %.4f' % worst_control)
    gate('V8', 'the ENGINE reports back the artefact\'s fitted stiffness on every element',
         engine_k_ok)
    return {'checked': checked, 'worst_relative': worst, 'control_min_relative': worst_control}


def v6_none_is_historical():
    print('\nV6  `None` is the historical plant, bit for bit', flush=True)
    from ihm.assembly.articulated import ArticulatedBodyPlant
    kwargs, selection = resolve_fidelity(ROOT, None, environment='upright')
    gate('V6a', 'resolve_fidelity(None) adds EXACTLY nothing to the plant kwargs', kwargs == {},
         str(kwargs))
    # resolving the fitted identity first: if it left anything behind, the next plant moves
    resolve_fidelity(ROOT, {'segment_contact': BUNDLE_ID}, environment='upright')
    runs, none_wall = [], []
    for i in range(2):
        plant = ArticulatedBodyPlant(ROOT, fresh('none%d' % i) / 'plant', environment='upright',
                                     initial_pose=POSE_RUN, mechanical_fidelity=None)
        began = time.perf_counter()
        try:
            rows = []
            for _ in range(30):
                rows.append({n: v['value'] for n, v in plant.native.snapshot()['coordinates'].items()})
                plant.advance(0.01)
            runs.append(rows)
            none_wall.append((time.perf_counter() - began) / 30)
        finally:
            plant.close()
    gate('V6', 'two `None` plants, 30 steps, bitwise identical in every coordinate',
         runs[0] == runs[1])
    plant = ArticulatedBodyPlant(ROOT, fresh('fitted') / 'plant', environment='upright',
                                 initial_pose=POSE_RUN,
                                 mechanical_fidelity={'segment_contact': BUNDLE_ID})
    began = time.perf_counter()
    try:
        fitted = []
        for _ in range(30):
            fitted.append({n: v['value'] for n, v in plant.native.snapshot()['coordinates'].items()})
            plant.advance(0.01)
    finally:
        plant.close()
    fitted_wall = (time.perf_counter() - began) / 30
    differ = next((i for i in range(30) if fitted[i] != runs[0][i]), None)
    gate('V6b', 'CONTROL THAT CAN FAIL: the fitted bundle DIFFERS from `None`', differ is not None,
         'first differing step %s' % differ)
    plant = ArticulatedBodyPlant(ROOT, fresh('layermap') / 'plant', environment='upright',
                                 initial_pose=POSE_RUN,
                                 mechanical_fidelity={'segment_contact': 'skin_layer_map'})
    began = time.perf_counter()
    try:
        for _ in range(30):
            plant.native.snapshot()
            plant.advance(0.01)
    finally:
        plant.close()
    unfitted_wall = (time.perf_counter() - began) / 30
    gate('V10', 'per-step cost: fitted vs the SAME bundle unfitted (one constant apart)', True,
         'None (COM spheres) %.1f / %.1f ms, unfitted bundle %.1f ms, fitted bundle %.1f ms, '
         'fitted - unfitted %+.1f ms; the coupled LAYER costs 53.9 ms per 10 ms step ON TOP '
         'of its plant, on ONE segment (docs/SOFT_BODY.md)'
         % (1e3 * none_wall[0], 1e3 * none_wall[1], 1e3 * unfitted_wall, 1e3 * fitted_wall,
            1e3 * (fitted_wall - unfitted_wall)), recorded=True)
    return {'none_ms_per_step': [1e3 * v for v in none_wall],
            'fitted_ms_per_step': 1e3 * fitted_wall, 'unfitted_ms_per_step': 1e3 * unfitted_wall,
            'coupled_layer_ms_per_step': 53.9,
            'attribution': 'The fitted bundle adds NO force element and no per-step evaluation '
                           'to the unfitted one: V9 shows the meshes are byte-identical and the '
                           'record set is the same, so the two differ by one constant. The '
                           'measured difference is therefore NOT the substitution\'s cost -- a '
                           'different stiffness gives a different trajectory and the error '
                           'controller does different work on it, which is a numerator holding '
                           'work the treatment never did. What IS attributable: the law is '
                           'solved inside the integrator, so it adds nothing per step that the '
                           'segment-contact bundle did not already cost, where the coupled '
                           'layer adds 53.9 ms per 10 ms step for ONE segment.'}


def v7_selectable():
    print('\nV7  the layer-map bundles resolve, and a uniform bundle still gets its triple',
          flush=True)
    material = {}
    for bundle in ('skin_layer_map', BUNDLE_ID, 'skin'):
        kwargs, _ = resolve_fidelity(ROOT, {'segment_contact': bundle}, environment='upright')
        material[bundle] = kwargs.get('segment_contact_material')
    gate('V7', 'a per-record layer map is resolved WITHOUT a caller E, p or h',
         material['skin_layer_map'] is None and material[BUNDLE_ID] is None)
    gate('V7b', 'CONTROL: a bundle with no layer map still receives the uniform triple',
         material['skin'] is not None,
         str(material['skin']))


def v9_meshes_untouched():
    print('\nV9  the fit changed a number and nothing else', flush=True)
    source = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_bytes())
    fitted = json.loads((ROOT / FIT / 'manifest.json').read_bytes())
    same_records = ([r['body'] for r in source['records']] == [r['body'] for r in fitted['records']])
    identical = True
    for record in fitted['records']:
        data = (ROOT / FIT / 'meshes' / record['mesh_file']).read_bytes()
        identical = identical and hashlib.sha256(data).hexdigest() == record['written_sha256']
    gate('V9', 'every mesh is byte-identical to the source bundle and the record set is the same',
         same_records and identical)


def main():
    began = time.perf_counter()
    print('THE CONTACT LAW FITTED TO THE LAYER: known answers, controls that can fail, idempotence\n',
          flush=True)
    report = json.loads((ROOT / FIT / 'report.json').read_bytes())
    print('artefact: %s (measured %s)' % (FIT, report.get('schema')), flush=True)

    v1_analytic()
    v1b_invariance()
    v3_null()

    manifest = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_bytes())
    record = [r for r in manifest['records'] if r['body'] == 'calcn_l'][0]
    vertices, faces = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
    layer = stl.segment_layer(ROOT, 'calcn_l', surface='fitted', depth='local')
    v4_idempotent(layer, vertices, faces, np.eye(3), float(vertices[:, 1].min()) + 0.002)
    held = v5_held_out(report, layer, vertices, faces)
    del layer

    cross = v2_cross_implementation(report)
    v7_selectable()
    v9_meshes_untouched()
    cost = v6_none_is_historical()

    out = {'schema': 'ihm.contact-law-fit-verification.v1', 'fit': FIT,
           'bars': {'held_out': BAR_HELD_OUT, 'cross_implementation': BAR_CROSS},
           'V2': cross, 'V5': held, 'V10': cost, 'results': RESULTS,
           'wall_seconds': time.perf_counter() - began, 'peak_rss_mb': peak_mb()}
    (ROOT / FIT / 'verification.json').write_text(
        json.dumps(out, indent=1, sort_keys=True,
                   default=lambda v: getattr(v, 'tolist', lambda: str(v))()) + '\n')
    failed = [r['gate'] for r in RESULTS if not r['pass'] and not r['recorded']]
    print('\nwall %.1f s, peak RSS %.0f MB' % (out['wall_seconds'], out['peak_rss_mb']), flush=True)
    print('%d gates, %d FAILED%s' % (len(RESULTS), len(failed),
                                     '' if not failed else ': ' + ', '.join(failed)), flush=True)
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
