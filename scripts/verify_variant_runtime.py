#!/usr/bin/env python3
"""Known answers for the live plant running the 25-body spine variant end to end.

`ArticulatedBodyPlant` was built for exactly 22 bodies and each layer refused
`articulated_spine_v1` in turn (docs/DISTRIBUTED_SURFACE_BINDING.md, "The 25-body
variant"). This file checks, on real engine sessions, that:

  A. the BASE plant is bit-identical to the committed behaviour before the variant work
     touched the skin layer -- `ihm/assembly/articulated.py` at PIN, loaded from git and
     run beside the current module; every frame field and every output file, float by
     float, over several steps. Two controls: the pinned source must really differ from
     the current one, and the comparator must catch a single 1-ulp change;
  B. a VARIANT plant constructs, steps and projects frames with its 25 segments, and
     calling the stateful projections twice at the same input gives the same bytes;
  C. its skin blend is the variant's own, the base blend is still refused for it (the
     control that shows the selection matters), and the variant blend is a pure
     refinement of the base one: merged back over the torso family it has the base
     argmax on every vertex;
  D. the skull follows `head`: in the display pose the head's named bones ride `head`
     and move under `head_extension` while nothing on torso/thorax/cervical moves; the
     SAME assertion on the base plant must fail (there the skull rides `torso`);
  E. in the force frame the skin follows too: under `head_extension` every skin vertex
     with no head weight stays exactly still, and every other vertex moves by exactly its
     head weight times the head's own motion.

Engine sessions run one at a time (each is prlimit-capped at 4 GB by the stream).

    PYTHONPATH=. .venv/bin/python scripts/verify_variant_runtime.py
"""
import gzip, hashlib, importlib.util, json, math, shutil, subprocess, sys, uuid
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
MASS = 77.6122029                      # ihm/body_constants.py: the plant's fed-weight mass
BASE_REGISTRATION = 'data/derived/mechanics/whole_body_arm26_v2/registration.json'
VARIANT_REGISTRATION = 'data/models/articulated_spine_v1/registration.json'
VARIANT_BINDING = 'data/derived/anatomy-segment-binding-articulated-spine-v1/binding.json'
PIN = '7ded91b'                        # the last commit before the skin layer was made per-plant
STEPS, DT = 5, 0.01
FAMILY = ('torso', 'thorax', 'cervical', 'head')
RECEIPT = ROOT / 'data/derived/variant-runtime-verification.json'
# Rounding bar for E: coordinates are <= ~2 m, float64 eps is 2.2e-16, and a posed vertex is
# a sum of 25 weighted multiply-adds, so rounding is ~1e-15 m, while a wrong segment or weight
# moves a vertex by millimetres. 1e-12 m sits three orders above rounding and nine below that.
RIGID_TOL_M = 1e-12


def check(name, ok, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {name}{"  " + detail if detail else ""}', flush=True)
    return bool(ok)


def diff(a, b, path='', out=None, limit=20):
    """every leaf that differs, floats compared bit for bit (sign of zero included)."""
    if out is None:
        out = []
    if len(out) >= limit:
        return out
    if type(a) is not type(b):
        out.append((path, 'type', type(a).__name__, type(b).__name__)); return out
    if isinstance(a, dict):
        if list(a) != list(b):
            out.append((path, 'keys', sorted(map(str, set(a) ^ set(b)))[:10] or 'order'))
        for k in a:
            if k in b:
                diff(a[k], b[k], f'{path}/{k}', out, limit)
    elif isinstance(a, (list, tuple)):
        if len(a) != len(b):
            out.append((path, 'len', len(a), len(b))); return out
        for i, (x, y) in enumerate(zip(a, b)):
            diff(x, y, f'{path}[{i}]', out, limit)
    elif isinstance(a, float):
        same = (a == b and math.copysign(1, a) == math.copysign(1, b)) or (math.isnan(a) and math.isnan(b))
        if not same:
            out.append((path, a, b))
    elif isinstance(a, np.ndarray):
        if a.shape != b.shape or a.dtype != b.dtype or a.tobytes() != b.tobytes():
            out.append((path, 'array differs'))
    elif a != b:
        out.append((path, repr(a)[:80], repr(b)[:80]))
    return out


def pinned_module():
    import ihm.assembly  # noqa: F401  -- the package the pinned module's relative imports resolve in
    src = subprocess.run(['git', 'show', f'{PIN}:ihm/assembly/articulated.py'], cwd=ROOT,
                         capture_output=True, check=True).stdout
    spec = importlib.util.spec_from_loader('ihm.assembly._articulated_pinned', loader=None)
    mod = importlib.util.module_from_spec(spec)
    mod.__package__ = 'ihm.assembly'
    exec(compile(src, f'<{PIN}:ihm/assembly/articulated.py>', 'exec'), mod.__dict__)
    return mod, src


def run(plant_type, registration, *, keep=False):
    """construct, STEPS steps, one re-projection; every frame and every output file."""
    out = ROOT / f'data/derived/scratch-verifyvariant-{uuid.uuid4().hex[:8]}'
    plant = plant_type(ROOT, out, environment='upright', target_mass_kg=MASS,
                       augmented_registration=registration, display_pose='opensim')
    try:
        rec = {'frames': [plant.snapshot()]}
        for _ in range(STEPS):
            rec['frames'].append(plant.advance(DT))
        rec['reproject'] = plant._project(plant.native.snapshot(), 0.)
        rec['files'] = {p.name: p.read_bytes() for p in sorted(out.iterdir()) if p.is_file()}
    except BaseException:
        plant.close(); shutil.rmtree(out, ignore_errors=True); raise
    if keep:
        return rec, plant, out
    plant.close(); shutil.rmtree(out, ignore_errors=True)
    return rec, None, None


def load_blend(relative):
    return json.loads(gzip.decompress((ROOT / relative).read_bytes()))


def bodies_at(kin, q):
    return {b: {'transform_ground': T.tolist()} for b, T in kin.forward(q).items()}


def main() -> int:
    from ihm.assembly import articulated
    from ihm.assembly.continuous_surface_binding import ContinuousSurfaceBinding, DEFAULT_ASSET
    ok = True
    receipt = {'pin': PIN, 'steps': STEPS, 'dt_s': DT, 'mass_kg': MASS}

    # ---------------------------------------------------------------- A. base bit-identity
    print(f'A. base plant against the committed behaviour ({PIN})', flush=True)
    pinned, pinned_src = pinned_module()
    current_src = (ROOT / 'ihm/assembly/articulated.py').read_bytes()
    ok &= check('control: the pinned module is a DIFFERENT source, without the per-plant skin map',
                hashlib.sha256(pinned_src).digest() != hashlib.sha256(current_src).digest()
                and not hasattr(pinned, 'SURFACE_BINDINGS') and hasattr(articulated, 'SURFACE_BINDINGS'))
    old, _, _ = run(pinned.ArticulatedBodyPlant, BASE_REGISTRATION)
    new, base_plant, base_out = run(articulated.ArticulatedBodyPlant, BASE_REGISTRATION, keep=True)
    try:
        d = diff(old, new)
        n_frames = len(new['frames']) + 1
        ok &= check('base plant: every frame field and output file bit-identical', not d,
                    f'{n_frames} frames ({STEPS} steps of {DT*1e3:.0f} ms + re-projection), '
                    f'{len(new["files"])} files, {len(d)} differences' + (f': {d[:3]}' if d else ''))
        ok &= check('base plant uses the base blend, as it always did',
                    base_plant.anatomy_plant == 'engineering_stance_v1'
                    and base_plant.surface_binding_asset == DEFAULT_ASSET
                    and len(base_plant.surface_binding.segments) == 22)
        # the comparator must be able to fail: one ulp in one centroid of one frame
        e = new['frames'][3]['entities']; k = next(iter(e)); x = e[k]['centroid_m'][1]
        e[k]['centroid_m'][1] = float(np.nextafter(x, np.inf))
        caught = diff(old, new)
        e[k]['centroid_m'][1] = x
        ok &= check('control: a single 1-ulp change is caught', len(caught) == 1, str(caught[:1]))
        receipt['base'] = {'frames_compared': n_frames, 'files_compared': sorted(new['files']),
                           'differences': len(d), 'ulp_control_caught': len(caught)}

        # the base half of D's control: the skull rides torso there
        vb = json.loads((ROOT / VARIANT_BINDING).read_text())
        skull = [i for i in vb['segment_named_bones']['head'] if i in base_plant.display_poser.index]
        bposer = base_plant.display_poser
        base_segs = {bposer.segments[bposer.segment_of[bposer.index[i]]] for i in skull}
    finally:
        base_plant.close(); shutil.rmtree(base_out, ignore_errors=True)

    # ---------------------------------------------------------------- B. the variant runs
    print('B. the variant plant, end to end', flush=True)
    vrec, plant, vout = run(articulated.ArticulatedBodyPlant, VARIANT_REGISTRATION, keep=True)
    try:
        native = plant.native.snapshot()
        segs = sorted(native['bodies'])
        ok &= check('variant plant constructs on its own anatomy plant',
                    plant.anatomy_plant == 'articulated_spine_v1' and len(segs) == 25,
                    f'{plant.anatomy_plant}, {len(segs)} native bodies')
        ok &= check('the force frame registers all 25', sorted(plant.registration.groups) == segs)
        times = [f['time_s'] for f in vrec['frames']]
        ok &= check(f'steps: time advances {STEPS} x {DT} s', np.allclose(np.diff(times), DT, rtol=0, atol=1e-12),
                    f'{times}')
        last = vrec['frames'][-1]
        ok &= check('projected frames carry 25 skin transforms and 25 display segments',
                    sorted(last['surface_transforms']) == segs and sorted(last['display_pose']['segments']) == segs)
        cents = np.array([v['centroid_m'] for v in last['entities'].values()])
        ok &= check('every entity centroid finite after stepping', np.isfinite(cents).all(), f'{len(cents)} entities')
        a, b = plant._project(native, 0.), plant._project(native, 0.)
        ok &= check('idempotent: _project twice at the same state', not diff(a, b))
        a, b = plant._display_pose(native), plant._display_pose(native)
        ok &= check('idempotent: _display_pose twice at the same state', not diff(a, b))
        receipt['variant'] = {'anatomy_plant': plant.anatomy_plant, 'bodies': segs, 'times_s': times,
                              'registration_extension': plant.registration_extension,
                              'files': sorted(vrec['files'])}

        # ------------------------------------------------------------ C. the skin blend
        print('C. the variant skin blend', flush=True)
        vasset = articulated.SURFACE_BINDINGS['articulated_spine_v1']
        ok &= check('the variant plant loads the variant blend',
                    plant.surface_binding_asset == vasset and [s['id'] for s in plant.surface_binding.segments] == segs)
        try:
            ContinuousSurfaceBinding.from_root(ROOT, plant.registration)
            refused = ''
        except ValueError as exc:
            refused = str(exc)
        ok &= check('control: the base blend is still REFUSED for the variant registration',
                    'registration supports changed' in refused, repr(refused))
        one = ContinuousSurfaceBinding.from_root(ROOT, plant.registration, asset=vasset)
        two = ContinuousSurfaceBinding.from_root(ROOT, plant.registration, asset=vasset)
        ok &= check('idempotent: the variant blend loaded twice',
                    one.identity == two.identity == plant.surface_binding.identity
                    and one.weights.tobytes() == two.weights.tobytes())
        base_blend, var_blend = load_blend(DEFAULT_ASSET), load_blend(vasset)
        bs = [s['id'] for s in base_blend['segments']]; vs = [s['id'] for s in var_blend['segments']]
        wb = np.asarray(base_blend['weights']); wv = np.asarray(var_blend['weights'])
        limbs = [s for s in bs if s != 'torso']
        same_support = all(base_blend['segments'][bs.index(s)] == var_blend['segments'][vs.index(s)] for s in limbs)
        limb_dw = max(float(np.abs(wb[:, bs.index(s)] - wv[:, vs.index(s)]).max()) for s in limbs)
        family = sum(wv[:, vs.index(s)] for s in FAMILY)
        family_dw = float(np.abs(family - wb[:, bs.index('torso')]).max())
        merged_v = np.stack([wv[:, vs.index(s)] for s in limbs] + [family], 1).argmax(1)
        merged_b = np.stack([wb[:, bs.index(s)] for s in limbs] + [wb[:, bs.index('torso')]], 1).argmax(1)
        disagree = int((merged_v != merged_b).sum())
        ok &= check('refinement: the 21 non-torso supports are the base ones', same_support)
        ok &= check('refinement: merged over the torso family, the argmax is the base argmax on every vertex',
                    disagree == 0, f'{disagree} of {len(wv)} disagree; limb columns differ by <= {limb_dw:.2e}, '
                    f'family sum vs base torso by <= {family_dw:.2e} (CG roundoff)')
        was_torso = wb.argmax(1) == bs.index('torso')
        now = {s: int(((wv.argmax(1) == vs.index(s)) & was_torso).sum()) for s in FAMILY}
        print(f'    base-torso vertices ({int(was_torso.sum())}) now: {now}')
        receipt['blend'] = {'asset': vasset, 'sha256': hashlib.sha256((ROOT / vasset).read_bytes()).hexdigest(),
                            'binding_identity': plant.surface_binding.identity, 'vertices': len(wv),
                            'limb_column_max_abs_diff': limb_dw, 'family_sum_max_abs_diff': family_dw,
                            'merged_argmax_disagreements': disagree, 'base_torso_vertices': int(was_torso.sum()),
                            'base_torso_vertices_now_argmax': now}

        # ------------------------------------------------------------ D. the skull follows head
        print('D. the display pose moves the skull with head', flush=True)
        poser = plant.display_poser
        lo, hi = poser.kin.coordinates['head_extension']['range']
        q_ref = dict(poser.reference_pose)
        q_head = dict(q_ref, head_extension=hi)      # the model's own upper limit
        M0 = np.array(plant._display_pose({'bodies': bodies_at(poser.kin, q_ref)})['segment_motion'])
        M1 = np.array(plant._display_pose({'bodies': bodies_at(poser.kin, q_head)})['segment_motion'])
        var_segs = {poser.segments[poser.segment_of[poser.index[i]]] for i in skull}
        ok &= check("the head's named bones ride head in the variant", var_segs == {'head'},
                    f'{len(skull)} bones on {sorted(var_segs)}')
        ok &= check('control: the SAME assertion fails on the base plant', base_segs != {'head'},
                    f'base: {len(skull)} bones on {sorted(base_segs)}')
        idx = np.array([poser.index[i] for i in skull])
        rest = poser.rest_centroid[idx]
        h = poser.seg_index['head']
        moved = np.linalg.norm((rest @ M1[h, :3, :3].T + M1[h, :3, 3]) - (rest @ M0[h, :3, :3].T + M0[h, :3, 3]), axis=1)
        still = [s for s in poser.segments if s != 'head']
        other = max(float(np.abs(M1[poser.seg_index[s]] - M0[poser.seg_index[s]]).max()) for s in still)
        ok &= check(f'head_extension {hi:.4f} rad moves every skull bone, and moves nothing else',
                    moved.min() > 0 and other == 0.0,
                    f'skull centroids moved {moved.min()*1e3:.1f}-{moved.max()*1e3:.1f} mm; '
                    f'max change of any other segment motion {other:.1e}')
        receipt['display'] = {'head_extension_rad': hi, 'skull_bones': len(skull),
                              'skull_moved_mm': [float(moved.min() * 1e3), float(moved.max() * 1e3)],
                              'other_segments_max_change': other, 'base_skull_segments': sorted(base_segs)}

        # ------------------------------------------------------------ E. skin in the force frame
        print('E. the skin follows head in the force frame', flush=True)
        sb = plant.surface_binding
        full = sb.sampled_binding(np.arange(len(sb.rest)))
        def skin(q):
            ent = plant.registration.project({'bodies': bodies_at(poser.kin, q)})
            return sb.project(full, sb.frame(ent))
        x0, x1 = skin(q_ref), skin(q_head)
        hcol = [s['id'] for s in sb.segments].index('head')
        w = sb.weights[:, hcol]
        off = w == 0
        off_move = float(np.abs(x1[off] - x0[off]).max())
        ok &= check('every skin vertex with no head weight stays EXACTLY still', off_move == 0.0,
                    f'{int(off.sum())} vertices, max move {off_move:.1e} m')
        # No vertex is wholly on head: the screened diffusion leaves every head vertex some
        # cervical/thorax weight above the 1e-8 tail cut. So the known answer is the blend's
        # own structure with only head's station moving: x1 - x0 == w_head * (s_head(q1) - s_head(q0)).
        def head_station(q):
            ent = plant.registration.project({'bodies': bodies_at(poser.kin, q)})
            return sb.stations(full, sb.frame(ent))[:, hcol]
        s0, s1 = head_station(q_ref), head_station(q_head)
        on = w > 0
        predicted = w[on, None] * (s1[on] - s0[on])
        residual = float(np.abs((x1[on] - x0[on]) - predicted).max())
        top = sb.weights.argmax(1) == hcol
        head_move = np.linalg.norm(x1[top] - x0[top], axis=1) if top.any() else np.zeros(1)
        rigid_move = np.linalg.norm(s1[top] - s0[top], axis=1) if top.any() else np.ones(1)
        ok &= check('every head-weighted vertex moves by exactly its head weight times the head motion',
                    residual < RIGID_TOL_M, f'{int(on.sum())} vertices, residual <= {residual:.1e} m')
        ok &= check('every vertex whose largest weight is head moves', top.any() and head_move.min() > 0,
                    f'{int(top.sum())} vertices moved {head_move.min()*1e3:.1f}-{head_move.max()*1e3:.1f} mm, '
                    f'following {np.median(head_move/rigid_move)*100:.1f}% (median) of the rigid head motion; '
                    f'head weight there {w[top].min():.3f}-{w[top].max():.6f}')
        receipt['force_frame_skin'] = {'no_head_weight_vertices': int(off.sum()), 'no_head_weight_max_move_m': off_move,
                                       'head_weighted_vertices': int(on.sum()), 'blend_residual_m': residual,
                                       'head_argmax_vertices': int(top.sum()),
                                       'head_argmax_moved_mm': [float(head_move.min() * 1e3), float(head_move.max() * 1e3)],
                                       'head_argmax_median_fraction_of_rigid': float(np.median(head_move / rigid_move)),
                                       'head_argmax_head_weight_range': [float(w[top].min()), float(w[top].max())]}
    finally:
        plant.close(); shutil.rmtree(vout, ignore_errors=True)

    receipt['all_pass'] = bool(ok)
    RECEIPT.write_text(json.dumps(receipt, indent=2, default=lambda v: v.tolist() if hasattr(v, 'tolist') else str(v)) + '\n')
    print(f'\n{"ALL PASS" if ok else "FAILED"}  receipt: {RECEIPT.relative_to(ROOT)}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
