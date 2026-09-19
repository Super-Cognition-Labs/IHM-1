#!/usr/bin/env python3
"""The soft-tissue layer IN the plant's integration loop: known answers, controls, and what
staleness costs.

docs/ACTUATION_STAGES.md's "fully present participant" mode is this coupling.  The layer
existed, deformed and pushed back on its bone, and nothing it computed reached the plant
(docs/WORKBENCH_AUTHENTICITY.md 2.1).  This battery judges the coupling that closes it.

THE COST THAT SHAPES THE DESIGN.  A loaded, warm-started solve costs 136-159 ms against the
plant's 10 ms step (docs/SOFT_BODY.md, RT1h recorded FAILED).  One solve per step is not
available.  The coupling is therefore: contact-gated every step (an unloaded segment costs one
matrix-vector product and emits nothing), sub-cycled when loaded (re-solve every
`resolve_interval_s`, hold the reaction in between), and warm-started across steps.

================================================================================
PRE-REGISTRATION.  Everything below was written before any gated run.  Bars are fixed here.
================================================================================

WHAT WAS ALREADY MEASURED WHEN THIS WAS WRITTEN, so nothing here is chosen by a result it
has already seen.  Five uncoupled probes and one coupled SMOKE run (to find crashes) preceded
it; all of them are disclosed, and none of them produced a staleness or cadence number:

  P1 At the plant's own stance the foot's SKIN sits 36.3 mm above the floor while the source
     foot contact spheres carry 616 N of the body's 761 N.  In a free drop from pelvis_ty=1.06
     the peak foot load is 770 N and the skin still clears the floor by 32.4 mm.  **The
     coupled layer on `calcn_l` therefore returns exactly 0 N in every upright trajectory the
     scaffold can produce.**  This is the skin bundle's own documented caveat
     (plant_options.SEGMENT_CONTACT_BUNDLES['skin']), now measured for the deformable layer.
  P2 At the flat stance pose `calcn_l`'s layer carries 0.119 / 0.430 / 1.278 / 2.866 / 5.574 N
     at 0.5 / 1 / 1.5 / 2 / 2.5 mm and LEAVES THE CONSTITUTIVE DOMAIN at 3.0 mm (min J 0.180).
     Its rigid core is 3.76 mm under its lowest skin point there.  Body weight is 761 N: this
     layer cannot carry a standing body, and no cadence changes that.
  P3 Over a 200-step uncoupled collapse from pelvis_ty=1.03, for every anchored segment, the
     number of steps in which its SKIN is below the floor while its own engine contact element
     is NOT.  `ulna_l` is the longest such window (69 steps), `ulna_r` 99, `tibia_r` 33,
     `calcn_l` and `calcn_r` ZERO.  A segment in that window is carried by the layer ALONE,
     with no engine contact to double-count against.  That is why the fixture below is
     `ulna_l`, and the choice is a measurement, not a preference.
  P4 pelvis_ty = 1.03 is the scaffold's own near-stance height (616 N under the left foot at
     t=0 against 761 N of weight); the model's default 0.93 starts the foot 63.7 mm inside the
     floor.
  P5 SMOKE RUN, coupled, 140 steps at a 10 ms interval, run to find crashes before this file
     existed.  It reached contact at step 81, and terminated at step 104 with
     `SoftTissueLeftDomain` (min J 0.189).  What it establishes, and all it is used for here:
     the loaded window of this fixture is about 6 steps (60 ms) long, so a cadence longer than
     that cannot re-solve inside it.  Its force values are NOT used to set any bar below.

THE FIXTURE, for every arm.  Upright plant, `initial_pose={'pelvis_ty': 1.03}`, no actuation,
no controller, 10 ms steps, at most 140 steps, coupled segment `ulna_l` and nothing else.  The
scaffold falls; its left forearm reaches the floor; the layer is the only thing on that segment
that touches.  **This is the 22-segment scaffold collapsing.  It is not the body doing
anything, and no number here is a statement about a human forearm.**

ARMS.  U uncoupled (mechanical_fidelity None); U2 the uncoupled soft-tissue identity;
R1 coupled at 10 ms (the per-step reference); A2/A5/A10/A20/A50 coupled at 20/50/100/200/500 ms;
S85 coupled at 10 ms with the reaction scaled by 0.85.

THE CADENCE RULE, as built.  Every plant step, every coupled segment: one matrix-vector product
gives the layer's minimum gap and its CORE's minimum gap against the support.  Core inside the
support -> raise (not sub-cycled, and not clamped).  Skin clear -> exactly zero, no port, no
solve.  Otherwise re-solve if `resolve_interval_s` has elapsed since this segment's last solve
or if contact has just begun, else hold the last solved reaction: the world force vector, and
the centre of pressure as a station in the SEGMENT frame, so the station follows the segment.

WHAT EACH BAR IS AND WHERE IT COMES FROM
  10%   the staleness bar (X10).  The layer's own contact force is not converged to better than
        about 10-15% (docs/SOFT_BODY.md CV7-CV10: first order, Richardson puts the finest heel
        forces 8% and 15% above their limits).  A coupling whose staleness error exceeded the
        layer's own discretisation error would make the cadence the dominant error.  10% is
        inside that, and it is the only place this number comes from.
  0.85  the force scale of the uncertainty arm (X11), the same measurement read as a factor.
  1e-6  the wrench-delivery bars (X3) and the layer balance bar (X4): both are identities that
        hold exactly in exact arithmetic, and the layer's own K5 reports 2.3e-8.
  1e-12 the canonical force round trip (X1): an orthogonal basis and an offset, applied and
        undone.  Anything above float noise here is a frame error.
  bitwise  every gate that compares two runs that MUST be the same run.

GATES
  X0  DETERMINISM.  Two identical uncoupled runs are bitwise identical in every coordinate at
      every step.  If this fails every comparison below is void and the battery stops.
  X1  The canonical force round trip: `registration.force` returns the source body, station and
      force the coupling handed it.  <= 1e-12 m and 1e-12 N.
  X2  SAME PATH.  A coupled plant whose layer is replaced by a FIXED wrench, and an uncoupled
      plant whose CALLER passes the equivalent canonical force through `forces=`, are bitwise
      identical over 20 steps.  The coupling invents no force path.
  X2b CONTROL THAT CAN FAIL.  The same comparison with the caller's force NEGATED must DIFFER.
      If it does not, X2 is comparing two plants that never felt the force at all.
  X3  The delivered point force reproduces the layer's wrench on every solve of R1:
      |p x F - M| and |M.F|/|F| are both <= 1e-6 of max(|M|, |F| r).  Frictionless contact
      against a half-space puts every nodal force along the normal, so the axial moment must be
      zero -- checked, not assumed, because a point force cannot carry it.
  X4  The layer's own balance identity on every coupled solve: balance_force_relative <= 1e-6.
  X5  OFF IS OFF, BIT FOR BIT.  U, U2 and R1 agree bitwise in every coordinate at every step up
      to the step before the layer first emits a port.  An unloaded coupled plant IS the
      historical plant.
  X6  IDEMPOTENCE.  `advance_state` called twice at the same (state, transforms, dt) returns
      bitwise the same ports and state -- once unloaded, once loaded.
  X7  Two identical R1 runs are bitwise identical.
  X8  CHECKPOINT.  From a loaded step: checkpoint, 5 steps, restore, 5 steps -> bitwise
      identical.  A coupling whose held reaction is not checkpointed fails this.
  X9  THE ERROR IS NOT SWALLOWED.  R1 terminates by raising SoftTissueBottomedOut or
      SoftTissueLeftDomain, and the plant's coordinates after the raise equal the coordinates
      before it, bitwise.  RECORDED: which, and at which step.
  X9b A pose whose rigid core is inside the support raises SoftTissueBottomedOut naming the
      segment and the depth, from the cheap per-step gate and not from a solve.
  X10 STALENESS.  On R1's own poses, the held reaction at cadence N is R1's solved reaction at
      the last step that cadence would have solved.  e_N = max_k |F_N(k) - F_1(k)| / max_k
      |F_1(k)| <= 0.10.  Peak-normalised, because a pointwise ratio at touchdown divides by
      nothing.
  X10b CONTROL THAT CAN FAIL.  e at "hold the first reaction forever" must EXCEED 0.10.  If it
      does not, the metric cannot see staleness and X10 means nothing.
  X11 CLOSED LOOP.  d_N = max over common steps of |origin_N - origin_1| for the coupled
      segment, over every arm's own run.  Bar: d_N <= d_U, where d_U is the divergence the
      layer's OWN force uncertainty produces (arm S85 against R1).  Measured, not assumed.  If
      d_U is 0 the bar is degenerate and X11 is VOID for every N.
  X12 A cadence arm must terminate within 1 step of R1, and for the same reason.  A cadence
      that changes when the tissue gives out has changed the answer, not just its resolution.
  VOID A cadence whose interval exceeds the loaded window never re-solves inside it.  Following
      this repo's own rule (scripts/verify_soft_tissue.py, third amendment), a gate whose
      solves carried no load is printed VOID and counts as FAILED.

AMENDMENT, written before any gated run of this battery (the smoke run P5 above is the only
coupled run that precedes it).  X10, X11 and X12 are stated above as gates.  They are SELECTION
CRITERIA: an interval that fails one of them is a measurement OF THAT INTERVAL, not a defect in
the coupling, and outcome (c) below is pre-registered as a result.  So each interval's three
criteria are printed PASS/FAILED and RECORDED -- they do not change the exit status -- and the
gated statements are the correctness ones (X0-X9b) plus two about the staleness instrument
itself, which CAN fail:
  X10z The replay must read EXACTLY 0.000000 at the per-step interval.  It is the reference
       replayed against itself; anything else is a bookkeeping error in the instrument.
  X10b Holding the first reaction forever must EXCEED the bar (as stated above).
The adopted cadence, and which of (a)/(b)/(c) obtained, are RECORDED.  No bar moves.

THE DECISION, fixed before the data.  The adopted cadence is the LARGEST interval passing X10,
X11 and X12.  Then, separately and without moving any bar, the amortised per-step cost of that
cadence is reported against the 10 ms plant step.  Three outcomes, all of them results:
  (a) the adopted cadence amortises below 10 ms  -> the coupling is real-time on this segment;
  (b) it does not                                -> the coupling is offered, is NOT real-time,
                                                    and the factor short is reported;
  (c) even 20 ms fails                           -> no sub-cycling is acceptable at these
                                                    costs, and that is the answer.

RECORDED, never gated: per-step wall cost unloaded and loaded, the amortised cost at each
cadence, the loading rate dF/dt of the contact, and the hold time a 10% bar would allow at that
rate (tau = 0.10 max|F| / max|dF/dt|), which is what carries this measurement to a contact that
develops at a different speed.

    PYTHONPATH=. prlimit --as=4294967296 .venv/bin/python -u scripts/verify_soft_tissue_coupled.py

Memory budget: 4 GiB address space (prlimit); peak is printed at the end.  One engine session
at a time: every arm builds its own plant, runs, and closes before the next starts.
"""
import json
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.assembly.articulated import ArticulatedBodyPlant                      # noqa: E402
from ihm.assembly import soft_tissue_layer as stl                              # noqa: E402
from ihm.assembly.soft_tissue_layer import (SoftTissueBottomedOut,             # noqa: E402
                                            SoftTissueLeftDomain, SoftTissueCoupling)

BODY = 'ulna_l'
POSE = {'pelvis_ty': 1.03}
DT = 0.01
STEPS = 140
INTERVALS = (0.010, 0.020, 0.050, 0.100, 0.200, 0.500)
STALENESS_BAR = 0.10
UNCERTAINTY_SCALE = 0.85
OUT = ROOT / 'data/derived/soft-tissue-coupled-v1'
WORK = ROOT / 'data/derived/soft-tissue-coupled-runs'

RESULTS = []


def gate(name, claim, ok, detail='', recorded=False):
    print('%-5s %-78s %s%s' % (name, claim, 'PASS' if ok else 'FAILED',
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


def fidelity(interval):
    if interval is None:
        return None
    return {'soft_tissue': {'id': 'layer_fitted_local_confined_coupled',
                            'segments': [BODY], 'resolve_interval_s': interval}}


class ScaledCoupling(SoftTissueCoupling):
    """Test double: the same coupling with its reaction multiplied by a fixed factor.

    This is how the layer's OWN force uncertainty is turned into a plant divergence, which is
    what X11's bar is made of.  It is a script-local subclass; production never scales.
    """

    scale = 1.0

    def _emit(self, body, held, rotation, origin):
        port = super()._emit(body, held, rotation, origin)
        if port is not None:
            port['force_n'] = (np.asarray(port['force_n']) * self.scale).tolist()
        return port


class FixedWrench(SoftTissueCoupling):
    """Test double: emit ONE fixed force at a fixed segment station, every step, no solve."""

    station = np.zeros(3)
    force = np.zeros(3)

    def advance_state(self, state, transforms, dt_s):
        body = self.bodies[0]
        T = np.asarray(transforms[body], float)
        point = T[:3, 3] + T[:3, :3] @ self.station
        return ({b: None for b in self.bodies},
                [{'body': body, 'point_m': point.tolist(), 'force_n': list(self.force)}],
                [{'body': body, 'in_contact': True, 'fixed_wrench': True}])


def trajectory(plant, steps, forces_at=None):
    """Coordinates, the coupled segment's origin, and the coupling receipt, step by step."""
    rows = []
    native = plant.native.snapshot()
    stop = None
    for k in range(steps):
        rows.append(sample(native, k))
        try:
            plant.advance(DT, forces=() if forces_at is None else forces_at(plant, native))
        except (SoftTissueBottomedOut, SoftTissueLeftDomain) as error:
            stop = {'step': k, 'kind': type(error).__name__, 'message': str(error)[:200]}
            break
        native = plant.native.snapshot()
        frame = plant.state.get('soft_tissue_coupling')
        rows[-1]['coupling'] = None if frame is None else frame['segments']
    return rows, stop


def sample(native, k):
    T = np.asarray(native['bodies'][BODY]['transform_ground'], float)
    return {'step': k, 'time_s': native['time_s'],
            'q': {n: v['value'] for n, v in native['coordinates'].items()},
            'origin_m': T[:3, 3].tolist(), 'coupling': None}


def same_q(a, b):
    """Bitwise equality of every coordinate over the steps both runs completed."""
    n = min(len(a), len(b))
    for i in range(n):
        if a[i]['q'] != b[i]['q']:
            return False, i
    return True, n


def divergence(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return float('nan'), 0
    d = max(float(np.linalg.norm(np.asarray(a[i]['origin_m']) - np.asarray(b[i]['origin_m'])))
            for i in range(n))
    return d, n


def reaction_series(rows):
    """The force the coupling DELIVERED at each step, and whether it solved there."""
    out = []
    for row in rows:
        entry = None if not row['coupling'] else row['coupling'][0]
        if entry is None:
            out.append((np.zeros(3), False, False))
        else:
            out.append((np.asarray(entry['force_n'], float), bool(entry['in_contact']),
                        bool(entry.get('solved', False))))
    return out


def replay(series, every):
    """What cadence `every` (in steps) would have held, on the reference run's own poses.

    A solve is due at contact onset and whenever `every` steps have passed since the last
    solve, which is exactly the rule the coupling runs.  `every = None` is "hold forever".
    """
    held = np.zeros(3)
    age = None
    out = []
    for force, contact, _ in series:
        if not contact:
            held, age = np.zeros(3), None
            out.append(held.copy())
            continue
        due = age is None or (every is not None and age >= every)
        if due:
            held, age = force.copy(), 0
        out.append(held.copy())
        age += 1
    return out


def build(tag, interval, *, steps=STEPS, double=None, forces_at=None):
    """One arm: a fresh plant, its trajectory, and its cost.  `double` swaps in a test double."""
    plant = ArticulatedBodyPlant(ROOT, fresh(tag), environment='upright', initial_pose=POSE,
                                 mechanical_fidelity=fidelity(interval))
    if double is not None:
        cls, attrs = double
        base = plant.soft_coupling
        swapped = cls(base.layers, plane_axis=base.axis, plane_sign=base.sign,
                      plane_value_m=base.plane, resolve_interval_s=base.interval_s,
                      method=base.method)
        for key, value in attrs.items():
            setattr(swapped, key, value)
        plant.soft_coupling = swapped
        plant.soft_state = swapped.new_state()
    began = time.perf_counter()
    rows, stop = trajectory(plant, steps, forces_at)
    wall = time.perf_counter() - began
    return plant, rows, stop, wall


def closed(tag, interval, **kw):
    plant, rows, stop, wall = build(tag, interval, **kw)
    plant.close()
    return rows, stop, wall


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    report = {'schema': 'ihm.soft-tissue-coupled.v1', 'body': BODY, 'pose': POSE, 'dt_s': DT,
              'steps': STEPS, 'staleness_bar': STALENESS_BAR,
              'uncertainty_scale': UNCERTAINTY_SCALE, 'intervals_s': list(INTERVALS)}
    began = time.perf_counter()

    print('\nX0  determinism of the substrate', flush=True)
    ua, _, _ = closed('u_a', None)
    ub, _, _ = closed('u_b', None)
    ok, n = same_q(ua, ub)
    if not gate('X0', 'two identical uncoupled runs agree bitwise in every coordinate', ok,
                'over %d steps' % n):
        print('\nX0 FAILED: every comparison in this battery would be void. Stopping.')
        (OUT / 'report.json').write_text(json.dumps({**report, 'gates': RESULTS}, indent=2) + '\n')
        return 1

    print('\nX1  the canonical force round trip', flush=True)
    plant, _, _, _ = build('roundtrip', 0.010, steps=0)
    native = plant.native.snapshot()
    rng = np.random.default_rng(0)
    worst_p = worst_f = 0.0
    owner_ok = True
    for _ in range(5):
        port = {'body': BODY, 'point_m': (rng.normal(size=3) * 0.1).tolist(),
                'force_n': (rng.normal(size=3) * 10).tolist()}
        canonical = plant._soft_force(port)
        back = plant.registration.force(canonical['id'], canonical['point_m'],
                                        canonical['force_n'], native)
        owner_ok &= back['body'] == BODY
        worst_p = max(worst_p, float(np.max(np.abs(np.asarray(back['point_m']) - port['point_m']))))
        worst_f = max(worst_f, float(np.max(np.abs(np.asarray(back['force_n']) - port['force_n']))))
    gate('X1', 'registration.force returns the source body, station and force it was handed',
         owner_ok and worst_p <= 1e-12 and worst_f <= 1e-12,
         'owner %s, %.2e m, %.2e N' % (owner_ok, worst_p, worst_f))
    anchor = plant.soft_force_ids[BODY]
    basis, offset = plant.registration.basis.copy(), plant.registration.global_map[:3, 3].copy()
    plant.close()

    print('\nX2  the coupling uses the caller\'s own force path, and a control that can fail',
          flush=True)
    station = np.array([0.0, 0.0, 0.05])
    force = np.array([0.0, 50.0, 0.0])

    def caller(sign):
        def at(plant, native):
            T = np.asarray(native['bodies'][BODY]['transform_ground'], float)
            point = T[:3, 3] + T[:3, :3] @ station
            return [{'id': anchor, 'point_m': (basis @ point + offset).tolist(),
                     'force_n': (basis @ (sign * force)).tolist()}]
        return at

    a2, _, _ = closed('fixed_coupled', 0.010, steps=20,
                      double=(FixedWrench, {'station': station, 'force': force}))
    b2, _, _ = closed('fixed_caller', None, steps=20, forces_at=caller(+1.0))
    c2, _, _ = closed('fixed_caller_neg', None, steps=20, forces_at=caller(-1.0))
    ok, n = same_q(a2, b2)
    gate('X2', 'a fixed wrench through the coupling == the same force through forces=', ok,
         'bitwise over %d steps' % n)
    bad, m = same_q(a2, c2)
    gate('X2b', 'CONTROL: the same comparison with the force negated must DIFFER', not bad,
         'first difference at step %d' % m)

    print('\nX5/X7  off is off, and the reference is reproducible', flush=True)
    u2, _, _ = closed('u2_uncoupled_identity', None)
    r1, stop1, wall1 = closed('r1', 0.010)
    r1b, stop1b, _ = closed('r1b', 0.010)
    series1 = reaction_series(r1)
    loaded = [i for i, (f, c, _) in enumerate(series1) if float(np.linalg.norm(f)) > 0]
    first = loaded[0] if loaded else len(r1)
    pre_ok = all(ua[i]['q'] == r1[i]['q'] and ua[i]['q'] == u2[i]['q'] for i in range(first))
    gate('X5', 'None, the uncoupled identity and the COUPLED plant agree bitwise while unloaded',
         pre_ok and first > 0, 'bitwise for %d steps, first port at step %d' % (first, first))
    ok, n = same_q(r1, r1b)
    gate('X7', 'two identical coupled runs agree bitwise', ok and stop1 == stop1b,
         'over %d steps, both stop %s' % (n, None if stop1 is None else stop1['kind']))

    print('\nX3/X4  the wrench that is delivered, and the layer\'s own balance', flush=True)
    worst_wrench = worst_axial = worst_balance = 0.0
    solves = 0
    for row in r1:
        entry = None if not row['coupling'] else row['coupling'][0]
        if entry is None or not entry.get('solved'):
            continue
        f = np.asarray(entry['force_n'], float)
        m_solved = np.asarray(entry['solved_moment_about_origin_nm'], float)
        m_delivered = np.asarray(entry['moment_about_origin_nm'], float)
        magnitude = float(np.linalg.norm(f))
        if magnitude == 0:
            continue
        solves += 1
        scale = max(float(np.linalg.norm(m_solved)), magnitude * 0.2519073143759591)
        worst_wrench = max(worst_wrench, float(np.linalg.norm(m_delivered - m_solved)) / scale)
        worst_axial = max(worst_axial, abs(float(f @ m_solved)) / magnitude / scale)
        worst_balance = max(worst_balance, float(entry['balance_force_relative']))
    gate('X3', 'the delivered point force reproduces the layer\'s wrench, axial moment included',
         solves > 0 and worst_wrench <= 1e-6 and worst_axial <= 1e-6,
         '%d solves, |pxF-M| %.2e, |M.F|/|F| %.2e' % (solves, worst_wrench, worst_axial))
    gate('X4', 'the layer\'s own force balance holds on every coupled solve',
         solves > 0 and worst_balance <= 1e-6, 'worst %.2e' % worst_balance)

    print('\nX9  the error is raised, not swallowed', flush=True)
    plant, rows9, stop9, _ = build('raise', 0.010)
    before = {n: v['value'] for n, v in plant.native.snapshot()['coordinates'].items()}
    again = None
    try:
        plant.advance(DT)
    except (SoftTissueBottomedOut, SoftTissueLeftDomain) as error:
        again = type(error).__name__
    after = {n: v['value'] for n, v in plant.native.snapshot()['coordinates'].items()}
    plant.close()
    gate('X9', 'the coupled run raises and the plant is left exactly where it was',
         stop9 is not None and again is not None and before == after,
         'stop %s at step %s; re-raised %s; plant unchanged %s'
         % (None if stop9 is None else stop9['kind'], None if stop9 is None else stop9['step'],
            again, before == after), recorded=False)

    print('\nX9b  the cheap per-step gate refuses a bottomed-out pose', flush=True)
    plant, _, _, _ = build('bottom', 0.010, steps=0)
    T = np.asarray(plant.native.snapshot()['bodies'][BODY]['transform_ground'], float)
    _, core = plant.soft_coupling.gaps(BODY, T[:3, :3], T[:3, 3])
    sunk = T.copy()
    sunk[1, 3] -= core + 0.01
    caught = None
    try:
        plant.soft_coupling.advance_state(plant.soft_state, {BODY: sunk.tolist()}, DT)
    except SoftTissueBottomedOut as error:
        caught = error
    gate('X9b', 'a pose with the rigid core inside the support raises, naming segment and depth',
         caught is not None and caught.body == BODY and abs(caught.depth_m - 0.01) < 1e-9,
         'None' if caught is None else str(caught)[:90])

    print('\nX6  call it twice at the same input', flush=True)
    clear = plant.native.snapshot()['bodies'][BODY]['transform_ground']
    pure = []
    for label, transform in (('unloaded', clear), ('loaded', _sink(clear, plant, BODY))):
        transforms = {BODY: transform}
        a = plant.soft_coupling.advance_state(plant.soft_state, transforms, DT)
        b = plant.soft_coupling.advance_state(plant.soft_state, transforms, DT)
        pure.append((label, _identical(a, b)))
    plant.close()
    gate('X6', 'advance_state is a function of its arguments, unloaded and loaded',
         all(ok for _, ok in pure), ', '.join('%s %s' % (k, v) for k, v in pure))

    print('\nX8  checkpoint and restore carry the held reaction', flush=True)
    plant = ArticulatedBodyPlant(ROOT, fresh('checkpoint'), environment='upright',
                                 initial_pose=POSE, mechanical_fidelity=fidelity(0.050))
    # stand INSIDE the sustained loaded window, so what is checkpointed is a held reaction
    mark = max((loaded[-1] if loaded else 0) - 4, 0)
    for _ in range(mark):
        plant.advance(DT)
    token = plant.checkpoint()
    branch = _two(plant)
    plant.restore(token)
    again = _two(plant)
    plant.release(token)
    plant.close()
    gate('X8', 'checkpoint from a loaded step, two steps, restore, two steps: bitwise identical',
         branch == again and len(branch) > 0,
         'standing at step %d, %d steps compared' % (mark, len(branch)))

    return _cadence(report, r1, stop1, wall1, series1, began)


def _two(plant):
    out = []
    for _ in range(2):
        try:
            out.append(_q(plant.advance(DT)))
        except (SoftTissueBottomedOut, SoftTissueLeftDomain) as error:
            out.append(type(error).__name__)
            break
    return out


def _q(frame):
    return {n: v['value'] for n, v in frame['joints'].items()}


def _sink(transform, plant, body):
    """The same segment transform, lowered until its skin is 2 mm into the support."""
    T = np.asarray(transform, float).copy()
    skin, _ = plant.soft_coupling.gaps(body, T[:3, :3], T[:3, 3])
    T[1, 3] -= skin + 0.002
    return T.tolist()


def _identical(a, b):
    sa, pa, _ = a
    sb, pb, _ = b
    if json.dumps(pa) != json.dumps(pb) or sorted(sa) != sorted(sb):
        return False
    for body in sa:
        if (sa[body] is None) != (sb[body] is None):
            return False
        if sa[body] is None:
            continue
        for key in sa[body]:
            x, y = sa[body][key], sb[body][key]
            if isinstance(x, np.ndarray):
                if not np.array_equal(x, y):
                    return False
            elif x != y:
                return False
    return True


def _cadence(report, r1, stop1, wall1, series1, began):
    """X10-X12: what staleness costs, and which cadence survives it."""
    magnitudes = np.array([float(np.linalg.norm(f)) for f, _, _ in series1])
    peak = float(magnitudes.max())
    contact_steps = int(sum(1 for _, c, _ in series1 if c))
    loaded_steps = int((magnitudes > 0).sum())
    walls = [row['coupling'][0]['wall_seconds'] for row in r1
             if row['coupling'] and row['coupling'][0].get('solved')]
    idle = [row['coupling'][0]['wall_seconds'] for row in r1
            if row['coupling'] and not row['coupling'][0]['in_contact']]

    print('\nX10  staleness on the reference run\'s own poses, and a control that can fail',
          flush=True)
    errors = {}
    for interval in INTERVALS:
        every = int(round(interval / DT))
        held = replay(series1, every)
        errors[interval] = float(max(np.linalg.norm(h - f) for h, (f, _, _) in
                                     zip(held, series1))) / peak if peak > 0 else float('nan')
    forever = replay(series1, None)
    e_forever = float(max(np.linalg.norm(h - f) for h, (f, _, _) in
                          zip(forever, series1))) / peak if peak > 0 else float('nan')
    for interval in INTERVALS:
        print('     %5.0f ms  e = %8.4f  (%d steps held)' % (interval * 1e3, errors[interval],
                                                             int(round(interval / DT))), flush=True)
    gate('X10z', 'KNOWN ANSWER: the replay reads exactly 0 at the per-step interval',
         errors[DT] == 0.0, 'e(10 ms) = %r' % errors[DT])
    gate('X10b', 'CONTROL: holding the first reaction forever must EXCEED the staleness bar',
         e_forever > STALENESS_BAR, 'e_forever %.4f against %.2f' % (e_forever, STALENESS_BAR))

    print('\nX11/X12  the closed loop: every cadence against the per-step reference', flush=True)
    arms = {0.010: (r1, stop1, wall1)}
    for interval in INTERVALS[1:]:
        arms[interval] = closed('cadence_%d' % round(interval * 1e3), interval)
    scaled = closed('uncertainty_85', 0.010,
                    double=(ScaledCoupling, {'scale': UNCERTAINTY_SCALE}))
    d_u, n_u = divergence(r1, scaled[0])
    print('     the layer\'s OWN force uncertainty (x%.2f) moves the segment %.3e m over %d steps'
          % (UNCERTAINTY_SCALE, d_u, n_u), flush=True)

    window_s = loaded_steps * DT
    rows = []
    for interval in INTERVALS:
        arm_rows, arm_stop, arm_wall = arms[interval]
        d, n = divergence(r1, arm_rows)
        arm_solves = sum(1 for row in arm_rows
                         if row['coupling'] and row['coupling'][0].get('solved'))
        void = interval > window_s
        same_stop = (stop1 is None) == (arm_stop is None) and (
            stop1 is None or (arm_stop['kind'] == stop1['kind']
                              and abs(arm_stop['step'] - stop1['step']) <= 1))
        rows.append({'interval_s': interval, 'steps_held': int(round(interval / DT)),
                     'staleness': errors[interval], 'divergence_m': d, 'common_steps': n,
                     'solves': arm_solves, 'wall_s': arm_wall,
                     'amortised_ms_per_step': 1e3 * arm_wall / max(1, len(arm_rows)),
                     'stop': arm_stop, 'same_stop': bool(same_stop), 'void': bool(void)})
        print('     %5.0f ms  e %8.4f  d %.3e m  solves %3d  stop %s' %
              (interval * 1e3, errors[interval], d, arm_solves,
               'none' if arm_stop is None else '%s@%d' % (arm_stop['kind'], arm_stop['step'])),
              flush=True)

    degenerate = not (d_u > 0)
    passing = []
    for row in rows:
        x10 = row['staleness'] <= STALENESS_BAR and not row['void']
        x11 = (not degenerate) and row['divergence_m'] <= d_u and not row['void']
        x12 = row['same_stop'] and not row['void']
        row.update(x10=x10, x11=x11, x12=x12)
        if x10 and x11 and x12:
            passing.append(row['interval_s'])
    for row in rows:
        tag = 'VOID (interval longer than the %.0f ms loaded window)' % (window_s * 1e3) \
            if row['void'] else ''
        gate('X10@%dms' % round(row['interval_s'] * 1e3),
             'held reaction within %.0f%% of the per-step reference' % (STALENESS_BAR * 100),
             row['x10'], 'e %.4f %s' % (row['staleness'], tag), recorded=True)
    for row in rows:
        gate('X11@%dms' % round(row['interval_s'] * 1e3),
             'plant divergence no worse than the layer\'s own force uncertainty',
             row['x11'], 'd %.3e m against d_u %.3e m%s'
             % (row['divergence_m'], d_u, '  VOID: d_u is 0' if degenerate else ''), recorded=True)
    for row in rows:
        gate('X12@%dms' % round(row['interval_s'] * 1e3),
             'terminates within one step of the reference, for the same reason', row['x12'],
             'stop %s' % ('none' if row['stop'] is None else
                          '%s@%d' % (row['stop']['kind'], row['stop']['step'])), recorded=True)

    adopted = max(passing) if passing else None
    solve_ms = float(np.median(walls) * 1e3) if walls else float('nan')
    amortised = None if adopted is None else \
        next(r['amortised_ms_per_step'] for r in rows if r['interval_s'] == adopted)
    rate = float(np.max(np.abs(np.diff(magnitudes))) / DT) if len(magnitudes) > 1 else 0.0
    tau = STALENESS_BAR * peak / rate if rate > 0 else float('inf')

    print('\nRECORDED (never gated)', flush=True)
    print('     loaded window            %.0f ms (%d steps in contact, %d carrying load)'
          % (window_s * 1e3, contact_steps, loaded_steps), flush=True)
    print('     peak reaction            %.4f N   (body weight is 761 N)' % peak, flush=True)
    print('     one loaded solve         %.1f ms median, %.1f ms max (%d solves)'
          % (solve_ms, max(walls) * 1e3 if walls else float('nan'), len(walls)), flush=True)
    print('     an unloaded segment      %.4f ms median over %d steps'
          % (float(np.median(idle)) * 1e3 if idle else 0.0, len(idle)), flush=True)
    print('     loading rate             %.2f N/s peak; a 10%% bar allows a hold of %.1f ms'
          % (rate, tau * 1e3), flush=True)
    print('     adopted cadence          %s' % ('NONE' if adopted is None
                                                else '%.0f ms' % (adopted * 1e3)), flush=True)
    if adopted is not None:
        print('     amortised cost there     %.1f ms per 10 ms plant step -> %s'
              % (amortised, 'REAL TIME' if amortised <= 10 else
                 '%.1fx short of real time' % (amortised / 10)), flush=True)
    if adopted is None or adopted <= DT:
        outcome = ('c', 'no sub-cycling cadence is acceptable on this contact: the reaction the '
                        'layer computes changes faster than any interval that would amortise its '
                        'cost. The coupling is offered at one solve per step and is not real time.')
    elif amortised is not None and amortised <= 1e3 * DT:
        outcome = ('a', 'the adopted cadence amortises below the plant step: real time here.')
    else:
        outcome = ('b', 'the adopted cadence is acceptable but does not amortise below the plant '
                        'step; the coupling is offered and is not real time.')
    print('     OUTCOME (%s)             %s' % (outcome[0], outcome[1]), flush=True)
    report['outcome'] = {'branch': outcome[0], 'statement': outcome[1]}

    report.update({'gates': RESULTS, 'cadence_rows': rows, 'amortised_ms_adopted': amortised, 'staleness': errors,
                   'staleness_hold_forever': e_forever, 'uncertainty_divergence_m': d_u,
                   'adopted_interval_s': adopted, 'loaded_window_s': window_s,
                   'contact_steps': contact_steps, 'loaded_steps': loaded_steps,
                   'peak_reaction_n': peak, 'loaded_solve_ms_median': solve_ms,
                   'unloaded_step_ms_median': float(np.median(idle)) * 1e3 if idle else 0.0,
                   'peak_loading_rate_n_per_s': rate, 'hold_allowed_by_bar_s': tau,
                   'reference_stop': stop1, 'plant_step_s': DT,
                   'peak_rss_mb': peak_mb(), 'wall_s': time.perf_counter() - began,
                   'basis': 'The 22-segment SCAFFOLD collapsing, with one deformable layer on '
                            'ulna_l. Not the body, and not a statement about a human forearm.'})
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'report.json').write_text(json.dumps(report, indent=2, default=str) + '\n')
    failed = [r['gate'] for r in RESULTS if not r['pass'] and not r['recorded']]
    print('\npeak RSS %.0f MB, wall %.0f s' % (peak_mb(), time.perf_counter() - began), flush=True)
    print('%d gates, %d FAILED%s' % (len(RESULTS), len(failed),
                                     '' if not failed else ': ' + ', '.join(failed)), flush=True)
    return 1 if failed else 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    finally:
        if WORK.exists():
            shutil.rmtree(WORK, ignore_errors=True)
