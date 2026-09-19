# Metabolic supply as a limit on excitation

18 Sep 2026. Closes the *mechanism* half of `WORKBENCH_AUTHENTICITY.md` §3.4 behind
an opt-in flag. **It is not a physiological fatigue model, and it has not been run
against a real BioGears interval that reports unmet energy.** The default is
unchanged: `metabolic_supply=None` still aborts on unmet energy, with the same
message, and is checked to reproduce the pre-change runtime bit for bit.

Code: `ihm/assembly/metabolic_supply.py`, wired into `EmbodiedRuntime.step`
(`ihm/assembly/embodied.py`). Select it with
`EmbodiedRuntime.from_workspace(..., metabolic_supply='energy_fraction')`.
Known answers: `scripts/verify_metabolic_supply.py`.

## 1. What "unmet" was, numerically

`embodied.py` aborted when `coupling.muscle_unmet_kcal > 1e-12`. That scalar is
written by the signed-muscle hook in the native variant the runtime loads
(`whole_body_integrity_gi_absorption` → parent `whole_body_integrity_signed_muscle_v2`,
whose `Tissue.cpp` sha256 is `8ca1eb1c…88add1d`, line 1484):

    ihm_signed::active()->unmet_muscle_kcal = tissueNeededEnergy_kcal > 0 ? tissueNeededEnergy_kcal : 0;

`tissueNeededEnergy_kcal` is the **muscle tissue's residual need in one 20 ms native
interval** after `Tissue::CalculateMetabolicConsumptionAndProduction` has drawn from,
in this order (line numbers in that `Tissue.cpp`):

| order | pool | rate / limit |
|---|---|---|
| start | need = basal share (0.8 × BMR × muscle blood-flow fraction) − 2.8% mandatory anaerobic (l.943) + 0.8 × generic-exercise demand + 0.5 × other above-basal + **the signed mechanical increment** (l.1165–1176) | |
| 1 | intracellular amino acids | 15–110 g/day by hormone factor (l.1191), O2-limited |
| 2 | intracellular TAG | rate factor 0–0.001 of the pool by hormone factor (l.1261), O2-limited |
| 3 | intracellular glucose, aerobic | 29.85 ATP/glucose (l.921), O2-limited |
| 4 | muscle glycogen, aerobic | 30.85 ATP/glucosyl (l.934), O2-limited, **no rate limit** (`TODO Glycogen should be consumed at a particular rate`, l.1376) |
| 5 | intracellular glucose, anaerobic | 2 ATP, → lactate |
| 6 | muscle glycogen, anaerobic | 3 ATP (l.935), → lactate; the mandatory anaerobic 2.8% is added back to the residual if this runs out |

ATP is not a pool. The ATP yields only convert reaction extents to kcal; there is no
ATP, PCr, Pi or intracellular H+ state (`docs/research/MUSCLE_SUPPLY_FEASIBILITY.md`
reached the same reading on 5 Sep).

Because branch 6 covers any residual while glycogen lasts, **unmet > 0 in practice
means intracellular glucose AND the whole-body muscle-glycogen scalar both ran out
inside that interval.** It is a residual of `coupling.muscle_requested_kcal`
(basal share + signed increment), not of the increment alone. The native kcal is
4184 J (`data/raw/physiology/biogears/share/etc/UCEDefs.conf:191`,
`calorie cal 4.184 J`).

Where the energy ledger lives, end to end:

1. **Plant**: cumulative Umberger/Uchida `muscle_metabolic_energy_j`,
   `signed_active_fiber_work_j`, `muscle_heat_energy_j` in the native mechanical engine.
2. **`embodied.py`**: 20 ms differences, minus the frozen zero-excitation reference
   (`metabolic_reference.json`, M0/W0/H0), through the τ = 2 s first-order lag;
   `metabolic_pending_energy_j` holds what the lag has not yet sent.
3. **Native**: `signed_step(ref, m, h, w)` → the `ihm_signed::Record`; heat goes in
   full to the Energy heat source *before* Tissue runs; the chemical increment joins
   the muscle's need above; `requested` and `unmet` come back as `coupling.*`.

Every retained run so far reports zero unmet (the 10 s supine run exchanged
0.618–23.918 W; `METABOLIC_REFERENCE.md`). Reaching unmet needs a glycogen-depleted
state — the abort is a guard on a regime nobody has yet driven the body into.

## 2. Why the law is an accounting constraint and not a fatigue model

The task was to find a published human fatigue law whose variables this engine
exposes. The engine exposes, for muscle: one whole-body glycogen scalar
(`muscle_glycogen_g`), vascular lactate, `fatigue_fraction`, `energy_deficit_w`,
and the requested/unmet pair. Candidates considered:

| law | needs | available here? |
|---|---|---|
| Xia & Frey Law 2008, *J Biomech* 41:3046 (three-compartment MA/MF/MR) | target load, time, fitted fatigue/recovery rates F, R per joint | no metabolic input at all — it is not a supply law |
| Callahan, Umberger & Kent-Braun 2016, *J Physiol* 594:3407 (in vivo fatigue from metabolites) | intracellular Pi, H+, H2PO4−, PCr/ATP kinetics | none of them exist in BioGears |
| Allen, Lamb & Westerblad 2008, *Physiol Rev* 88:287 (cellular mechanisms) | Pi, SR Ca2+ release, glycogen — mechanistic review | no quantitative law on these variables |
| Ørtenblad, Westerblad & Nielsen 2013, *J Physiol* 591:4405 (glycogen and fatigue) | local (sub-cellular) glycogen → SR Ca2+ release | association, no force-vs-glycogen law; BioGears glycogen is one whole-body scalar |
| BioGears `FatigueLevel` / `AchievedExerciseLevel` (Tissue.cpp l.1540–1559) | its own deficit / request | engine accounting, uncited in source — not published human work |

None maps. So, as the task specified for that case, what is implemented is **only the
strictly accounting cap: deliver the fraction of demand that supply covered**,
declared as an accounting constraint (`LAW` in the module, recorded in the session
manifest and every frame). BioGears' own `AchievedExerciseLevel` is the same
construct at whole-body scope; that is a coincidence of bookkeeping, not evidence.

## 3. The law, exactly

After exchange *k*, with `sent_k` the chemical increment sent (J, after the lag) and
`U_k` the native unmet (kcal):

- unsupplied increment `u_k = min(U_k·4184, max(sent_k, 0))`. Once every branch is
  exhausted the residual is linear in the need, so this is the marginal attribution.
  The remainder `U_k·4184 − u_k` is the native's **own basal muscle deficit**; it is
  reported (`native_basal_unmet_j`) and never charged to mechanics.
- covered fraction `c_k = (sent_k − u_k)/sent_k`; `c_k = 1` when `u_k = 0` or
  `U_k ≤ 1e-12 kcal` (the old abort gate, kept as the dead-band).
- interval *k+1* delivers `c_k × requested` to **every** muscle. `c_k == 1` delivers
  the request untouched (no multiplication happens).

Three details that matter:

- **Every muscle, not just the commanded ones.** The native engine holds an
  uncommanded muscle's last excitation. A cap that scaled only the command would
  leave the rest excited at zero supply (K5d makes exactly this fail). When the cap
  binds, and on the first interval after it releases, the full map is sent, built
  from the *requested* held excitation, which the limiter tracks.
- The current-feedback controllers command inside `advance_feedback_exchange`, so
  the cap is applied by a plant proxy on `advance`/`advance_observation` only
  (`CappedActuationPlant`); positional actuation is refused so nothing can bypass it.
- Supply per muscle does not exist natively; the scale is uniform.
- A `motor_blocks` muscle normally holds its previous excitation. While the cap binds
  it is commanded too, at scale × its held request — the block stops new commands,
  it does not exempt a muscle from supply.

## 4. The ledger — what "conserved" means here

Energy that supply did not cover is **not dropped**. Each exchange splits the sent
triple into `drawn` and `unsupplied` (heat and work attributed in proportion; heat
was already injected in full by the native port), and checks, every exchange:

- `raw = drawn + unsupplied + pending`, cumulatively, for chemical, heat and work —
  three independently accumulated quantities (the plant's differences, the lag's
  pending, the native's shortfall). Failure raises and aborts the run.
- `chemical = heat + work` on the drawn and unsupplied triples, with the tolerance
  `measure_resting_metabolic_reference` uses (`1e-9·(1+|m|+|w|+|h|)`).

**The ledger closing is bookkeeping, not physics.** The body conserves energy
physically only while cumulative `unsupplied` is zero; each frame says which
(`physically_conserved`). Honest limits:

- **One interval late.** The native has no restorable preview, so the interval where
  supply first falls short has already been integrated.
- **It scales excitation, not energy.** Umberger demand is not linear in excitation,
  activation has memory, and the native sees demand through the 2 s lag — so the
  next interval's energy is not guaranteed to fit. The coverage fraction is exact
  in excitation; its effect on energy is measured by the next exchange, not assumed.
- **Dynamics unmeasured.** With the 2 s lag, full release on `c = 1` and
  re-engagement on the next shortfall, a limit cycle over roughly τ is expected. It
  has not been observed on a real run.

## 5. Known answers — what `scripts/verify_metabolic_supply.py` prints

Local, 16 tests, all passed (fixtures through the real `EmbodiedRuntime.step`; the
plant holds uncommanded excitations like the native engine; the native fixture's
request is 20 W basal + increment against a scripted capacity):

| | result |
|---|---|
| K1 `None` vs pre-change `embodied.py` (git blob `542f566a…`), 25 exchanges | frames, plant commands and native demands **identical** |
| K1 capped, ample supply, 25 exchanges | plant commands **bit-identical** to uncapped; scale 1.0 every interval; coverage min 1.0; unsupplied 0.0 |
| K2 zero supply | first interval full (scale 1, coverage → 0.0); every later interval commands `{a:0.0, b:0.0, c:0.0}` (c never commanded); active work increment **0.0 W exactly** |
| K3 84 exchanges, capacity cycling ample/60/20/0/35/ample W | worst closure residual **2.8e-14 J**; 17 distinct scales; 40 capped exchanges; cumulative raw 122.31 J = drawn 24.35 + unsupplied 17.04 + pending |
| K4 idempotence | two runtimes, 40 exchanges, identical frames; `assess_exchange`, the cap and `close` each called twice at the same input identical; `close` leaves the limiter unchanged |
| K5 controls made to fail | a 1 µJ pending leak **raises**; dropping the unsupplied channel **raises**; a cap of 1−1e-6 is **detected** by K1's comparison; a commanded-only cap **leaves c at 0.2** under zero supply |
| K6 `None`, zero supply | raises the original message; runtime marked failed |
| factory | mode and `LAW` in the manifest, module in source receipts; `None` adds nothing; unknown mode rejected before output is created |

**Recorded, not rescored:** the first run of K2 FAILED on an extra assertion I had
added — that the *chemical* increment is exactly 0 at zero excitation. It was
4.26e-14 W on the first capped interval: the fixture differences cumulative float
energies, so `(h_new − h_old)/dt − 100` is round-off. The pre-registered answer
("no active work") was exactly 0.0 on that run. The assertion now reads
`|chemical| ≤ 1e-12 W`; active work remains asserted exactly zero.

Existing regressions still pass: `scripts.verify_embodied_runtime` (16),
`scripts.verify_candidate_embodied_factory` + `scripts.verify_regional_embodied_factory` (7).

**Native mechanical plant (`--native`), 18 Sep 2026, passed.** One real
`NativeMechanicalStream` (supine, 77.6122029 kg, 80 muscles), 5 × 20 ms intervals per
arm from one checkpoint. 40 muscles commanded to 0.3; the other 40 left at the
engine's held default excitation (0.01 for every muscle):

| check | result |
|---|---|
| scale 1 vs the uncapped command: energies, all 80 excitations, every body transform | **bit-identical** |
| scale 0 vs an explicit all-zero command | **bit-identical** |
| scale 0: all 80 engine-reported excitations | **exactly 0.0** — including the 40 uncommanded muscles, which a commanded-only cap would leave at 0.01 |
| scale 0.5: engine excitation of every muscle vs 0.5 × requested held map | equal (to 1e-15 absolute) |

This shows that the cap's command reaches the real engine for every muscle and that
full coverage costs nothing. It does NOT show the cap's dynamics against a real
BioGears shortfall; nothing here drove the native physiology.

The first `--native` launch failed before starting an engine:
`NativeMechanicalStream` needs an output directory that does not yet exist, and the
script passed it the directory `mkdtemp` had just created. The script now passes a
fresh child of it. Report: `data/derived/metabolic-supply-verification/report.json`
(not committed; `data/derived` is gitignored).

## 6. What would make this a physiological result

- A run from a glycogen-depleted state, so the native reports unmet and the cap's
  real dynamics (the expected limit cycle, the size of the first-interval
  shortfall) are measured rather than argued.
- A restorable native preview (the contract in `ihm/coupling/muscle_supply.py`),
  which would remove the one-interval lateness.
- An intracellular metabolite state (Pi, H+, PCr) with one native owner; only then
  does a published force law such as Callahan et al. 2016 have its inputs.
