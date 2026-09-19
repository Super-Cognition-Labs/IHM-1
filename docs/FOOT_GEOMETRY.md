# the foot the scaling step left behind: one joint offset that was never scaled

`docs/FOOT_JOINTS.md` closed the toe-GeometryPath route as *"not viable as is"*
because installing it produced an impossible sign pattern — the dorsal extensor
and the plantar flexors sharing a sign — and named the cause: the toe muscles'
last calcaneal via points were scaled and the mtp joint offset was not. This
document repairs that geometry in a NEW variant, `data/models/corrected_foot_v1`,
and re-tests the toe joint on it.

`data/models/engineering_stance_v1` is **read and never written**. It is the
identified plant every recorded result stands on, and `FOOT_JOINTS.md`'s
recommendation is explicitly to leave it alone.

Everything below is the SCAFFOLD's — the 22-body plant. Nothing here is a
statement about a human foot.

    build the variant   OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_corrected_foot --check
    weld mtp in spine   .venv/bin/python -m scripts.build_mtp_welded_spine --check
    run                 cd <repo> && OPENBLAS_NUM_THREADS=1 nohup nice -n 10 \
                            .venv/bin/python -u -m scripts.run_corrected_foot_arms --part all \
                            > logs/corrected_foot_arms.log 2>&1 &
    report              data/models/corrected_foot_v1/corrected_foot_report.json
                        per-step traces under data/derived/corrected-foot-arms/ (gitignored)

## 1. What the scaling step did, and the one thing it did not

`subject_walk_scaled.osim` is `Rajagopal2016.osim` put through OpenSim's
ScaleTool. ScaleTool records the per-body factors it used in the model it emits,
as `<Mesh><scale_factors>`, and multiplies every `PhysicalOffsetFrame`
translation by **its own parent body's** factors.

Audited over the 44 joint offset frames the two files share
(`scripts/build_corrected_foot.py`, check C3 — it fails if this stops being true):

| class | n | |
|---|---:|---|
| translation is the zero vector | 20 | nothing to scale |
| on `ground`, no factors recorded | 1 | `ground_pelvis/ground_offset` |
| **generic × its own body's factors, to 1e-12** | **19** | the rule |
| re-placed by the knee's own spline rebuild | 2 | `walker_knee_{l,r}/femur_{l,r}_offset` |
| **byte-identical to the UNSCALED generic value** | **2** | **`mtp_{l,r}/calcn_{l,r}_offset`** |

So the mtp calcaneal offset — `0.1788 −0.002 ±0.00108` m — is the **one joint
offset in the file the scaling step missed.** The ankle's tibia offset was
scaled (`−0.01 −0.4 0` → `−0.011635086131 −0.46517425728 0`, tibia_r's
`1.1635086131 1.1629356432 1.1352204367`). The subtalar's talus offset was scaled
(talus_r's `1.2281617902 1.3418506218 1.2869909462`). The toe muscles' calcaneal
via points were scaled, by calcn_r's `1.155608663 1.1170848239 1.3052527858`. The
axis those points are meant to straddle was not moved with them.

**And the same two joints are the only ones whose `orientation` was not carried
over.** Rajagopal declares `-3.14159 ±0.619901 0` on *both* frames of `mtp_{l,r}`
— an oblique metatarsal break — and `subject_walk_scaled.osim` carries `0 0 0`,
which makes it a pure calcn-z hinge. Every other joint's orientation in the
scaled file is its parent's, to the byte (`ankle_r` `0.175895 −0.105208
0.0186622`, `subtalar_r` `−1.76819 0.906223 1.8196`, and so on).

`Rajagopal2016.osim` and `RajagopalLaiUhlrich2023.osim` declare the mtp frames
**identically** (check C2), so "Rajagopal's value" is not a choice between two
sources.

**Why nobody upstream would have noticed.** Rajagopal ships `mtp_angle` with
`<locked>true</locked>`. A locked hinge has no moment arm to get wrong and no axis
to point in the wrong direction. `subject_walk_scaled.osim` unlocks it
(`<locked>false</locked>`), and this repository's engine runs it free, so the
defect became reachable here and not there.

## 2. The corrected offset, and its derivation

    corrected translation = Rajagopal2016's own translation
                            × this model's own calcn <Mesh><scale_factors>,
                            componentwise
    corrected orientation = Rajagopal2016's own orientation, on BOTH frames of
                            the joint (angles are dimensionless and ScaleTool
                            carried every other joint's orientation over unchanged)

| | shipped | corrected | Δ |
|---|---|---|---:|
| `mtp_r` calcn offset x | 0.1788 | **0.20662282894439998** | **+27.8228 mm** |
| y | −0.002 | −0.0022341696478 | −0.2342 mm |
| z | +0.00108 | +0.001409673008664 | +0.3297 mm |
| `mtp_r` orientation | `0 0 0` | `−3.14159 +0.619901 0` | axis (0,0,1) → (0.580955, 2e−6, −0.813936) in calcn |
| `mtp_l` orientation | `0 0 0` | `−3.14159 −0.619901 0` | mirrored |

Both frames of the joint take the orientation, so the calcn→toes transform at
q = 0 is unchanged apart from that translation: only the hinge **axis** moves.
`scripts/build_corrected_foot.py` re-derives both halves from the two source files
at build time; nothing in it is typed, and `--check` re-derives them again.

**A consequence that is not a moment arm.** The mtp axis carries the
`toes_{l,r}` BODY with it, so the toes, their two source contact spheres
(`lateralToe`, `medialToe`) and the toe muscles' toes-frame points all move
27.82 mm distally. **The corrected foot is a longer foot, with its forefoot
contact 27.8 mm further forward.** That is a change to the contact model, and
`docs/FOOT_JOINTS.md`'s own Q1 result — where a 48.5 mm change in support-plane
position moved every ankle by more than a radian — is the reason it cannot be
waved through. It is why the arms below are a 2×2 and not a pair.

## 3. PRE-REGISTRATION — written and committed before any arm below was run

Nothing in this section is edited after the runs. Results go in their own
section, in a separate commit.

### The arms

| arm | model | path set | G | P |
|---|---|---|:-:|:-:|
| `base` | `engineering_stance_v1/model.osim` | shipped fitted (80) | 0 | 0 |
| `base_toe_paths` | the same, unchanged | shipped minus the 8 toe muscles | 0 | 1 |
| `corrected_fitted` | `corrected_foot_v1/model.osim` | shipped fitted (80) | 1 | 0 |
| `corrected_toe_paths` | the same | shipped minus the 8 toe muscles | 1 | 1 |

**G** = the geometry correction (which also moves the toes body and its two
contact spheres 27.8 mm distally). **P** = edl/ehl/fdl/fhl, both sides, on their
own `GeometryPath`s instead of the fitted polynomial — the route `1e29ee3` took
for the subtalar. `base_toe_paths` is `FOOT_JOINTS.md`'s `base_toe_geometry_paths`
rebuilt here, so the whole 2×2 runs without touching
`data/derived/foot-joint-arms`, which another worker is using.

`corrected_fitted` is the arm that carries the geometry and contact change with
**no new mtp moment arm**, because a fitted polynomial has no argument for a joint
it was not fitted over, whatever the geometry does. It is what separates the two
factors.

### A — the sign pattern, and the gate on it

Engine `moment_arms` at the rest pose, `free` environment, **called twice at the
same input and asserted equal** before any value is read.

**GATE S1**, on `corrected_toe_paths`, on **both** sides:

    sign(edl) == sign(ehl)            the two toe extensors agree
    sign(fdl) == sign(fhl)            the two toe flexors agree
    sign(edl) != sign(fdl)            extensors oppose flexors
    min |arm| >= 1.0 mm               no arm small enough to make the test vacuous

The floor is there because a sign test on an arm of 0.05 mm is a test of rounding.
1.0 mm is `docs/UPPER_BODY_ACTUATION.md`'s own "effectively unactuated" floor,
borrowed; it is **not measured here**, and it is declared as borrowed.

**S1 PASS** = the pattern is anatomically possible and the fix has installed toe
flexion and extension. **S1 FAIL** = *the fix is wrong*, and this document will
say so in those words. An impossible pattern after the correction is not a
smaller problem than an impossible pattern before it.

S1 says nothing about magnitudes being *right*. Two extensors and two flexors
with opposite signs is the weakest claim that can be false, which is the point.

### A prediction, made before the engine ran

The corrected arms were computed **offline from the XML alone** on 18 Sep 2026,
before any engine session touched the corrected model, as

    arm = a · ((p − T) × u)

with `a = R(orientation)·ẑ` the hinge axis in calcn, `T` the corrected calcn
offset, `p` the last calcaneal path point and `u` the unit vector from `p` to the
first toes-frame point mapped into calcn at q = 0. The same routine on the
**shipped** geometry reproduces the engine's recorded values
(−4.7804 / +24.1416 / +8.3729 / +10.3578 mm for edl/ehl/fdl/fhl, from
`foot_joints_report.json`) **to 0.01 mm**, up to one global sign that is the
engine's axis convention. That agreement is what licenses a bar on the prediction.

| muscle | predicted arm about `mtp_angle`, both sides |
|---|---:|
| edl | **−6.052 mm** |
| ehl | **−7.551 mm** |
| fdl | **+6.834 mm** |
| fhl | **+7.129 mm** |

### Instrument checks — declared before the run

**A failure in any of I1–I6, I8, I9 VOIDS the comparison with prior numbers.**
Internal comparisons within this run would still stand, and would be reported as
such.

| | check | bar | provenance of the bar |
|---|---|---|---|
| I1 | `moment_arms` called twice, bit-equal, every arm | exact | IHM-1 `CLAUDE.md`, *call it twice* |
| I2 | `base` reads **exactly 0** about mtp for all 8 toe muscles | exact 0 | `WORKBENCH_AUTHENTICITY.md` §0.2a, a known answer |
| I3 | `base_toe_paths` reproduces the recorded shipped-geometry arms | 1e−6 m | `foot_joints_report.json`, full precision |
| I4 | `soleus_r` about `ankle_angle_r` = −0.0497 m, every arm | 5e−5 m | `FOOT_JOINTS.md`'s own control |
| I5 | `corrected_fitted` **still** reads exactly 0 about mtp | exact 0 | a fitted polynomial cannot see a joint it was not fitted over |
| I6 | `corrected_toe_paths` reproduces the offline prediction | 0.05 mm | 5× the 0.01 mm agreement measured on the known answer |
| I7 | lumbar arms reproduce `UPPER_BODY_ACTUATION.md:110` (42.69 / −52.82 / −62.79 mm) | 0.05 mm | **NON-VOIDING** — see below |
| I8 | `base` reproduces 0.122880 / 0.117805 / 1.449794 / 0.219587 rad | 5e−5 rad | `foot_joints_report.json`, full precision |
| I9 | `base_toe_paths` reproduces ΔA +0.0784 / +0.0106 / −0.0574 / −0.0137 | 5e−4 rad | the same report |

**I7 is declared non-voiding in advance, and this is not a hedge chosen after the
fact.** `UPPER_BODY_ACTUATION.md` does not record the pose at which that table was
measured, and a moment arm is pose-dependent — the same document says so about
`pro_sup`. So a disagreement there is evidence about a pose, not about this plant.
It is recorded FAILED if it fails, and diagnosed.

### B — what the toe hinge costs now

`docs/FOOT_JOINTS.md`'s Q2 protocols, **imported from its own runner** rather than
restated, so that "the same protocol" is a fact about the call graph:

- **supine tonic** — supine, 0.02 excitation on every muscle, 200 × 10 ms, worst
  excursion past each declared range — **unstopped** and **stopped** (crawl.py's
  stop set: every declared range but the ±10 rad sentinels, 30 N·m/rad, 85.94
  N·m·s/rad, 0.35 rad);
- **crawl seed** — `crawl.SEED`, prone, upright environment, 3 s — **unstopped**
  and **stopped**. The full horizon runs; where `crawl.diverged()` would have
  ended a search the first such event is recorded, not acted on.

Measured on each arm: `mtp_angle` extrema and excursion past ±0.5236; ankle
excursion; worst excursion over every other declared coordinate; for the crawl,
pelvis forward travel, first divergence event, peak fibre velocity; wall clock.

**Baselines, the numbers to beat or match** (base plant, `foot_joints_report.json`):

| protocol | base `ankle_angle_r` |
|---|---:|
| supine unstopped | 0.1229 |
| supine stopped | **0.1178** |
| crawl 3 s unstopped | **1.4498** |
| crawl 3 s stopped | **0.2196** |

**Reading.** No gate is set on mtp travel — it is the measurement asked for.
*"Changes the ankle"* is reported when |ΔA| against `base` exceeds **0.022 rad**,
`FOOT_JOINTS.md`'s reporting line, carried over unchanged and still *a reporting
convention with no measured provenance*. The 2×2 is read as four differences on
A = max(ankle_r, ankle_l), all four reported:

    G at P=0   corrected_fitted     - base                 the geometry/contact change alone
    G at P=1   corrected_toe_paths  - base_toe_paths
    P at G=0   base_toe_paths       - base                 the toe arms alone, wrong geometry
    P at G=1   corrected_toe_paths  - corrected_fitted     the toe arms alone, right geometry

If G alone moves A by more than the reporting line, then **the 27.8 mm of extra
forefoot is doing the work and the moment arms are not**, and it will be said that
way. The simulation is deterministic, so there is no sampling error to clear; it
is also one nonlinear trajectory per arm, so a difference between two arms is a
property of those two trajectories and not of a population.

**Not measured, and not claimed:** upright stance, walking push-off, the
identified linearization. `linearization.npz` was solved on the base plant with
mtp free at its shipped offset, so neither arm here is the identified plant.

### C — mtp welded in the articulated spine variant

`FOOT_JOINTS.md`'s standing recommendation is to weld mtp in any NEW plant. The
one new plant this repository has is `articulated_spine_v1`, so:

    data/models/articulated_spine_v1/model_mtp_welded.osim
    data/models/articulated_spine_v1/passive_mtp_welded.xml
    data/models/articulated_spine_v1/registration_mtp_welded.json

built by `scripts/build_mtp_welded_spine.py`, beside the `*_muscled` pair already
there. `mtp_{l,r}` becomes a `WeldJoint` on the joint's own frames (rest pose =
the joint at q = 0) and the two `mtp_angle` passive terms (−25·q − 2·q̇) are
removed through a source override, because they name a coordinate that no longer
exists.

**`articulated_spine_v1/model.osim` is not written** (check W2 fails if it is).

**GATE F2 IS NOT RESCORED.** F2 FAILED on `registration.json` —
`thoracic_extension` 0.3940 against a 0.2202 bar — and it stays FAILED. The
welded model is a **different plant**; every number measured on it is a new
measurement of that plant. Whatever it prints, F2's verdict does not move. This is
written here, before the run, so that it cannot be read as a response to the
numbers.

Protocol: F2's, unchanged — supine, 0.02 tonic, 2 s, worst excursion past each
declared range, bar measured on the **base** model in the same run. Arms: `base`,
`spine` (the committed `registration.json`), `spine_mtp_welded`. Instrument
checks: the bar reproduces 0.2202 and `spine` reproduces `thoracic_extension`
0.3940, both to 5e−5.

Reported: every coordinate's worst excursion on all three, the per-coordinate
difference the weld makes, and which coordinates are over the bar on each.
