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

## Results — run 18 Sep 2026, after the pre-registration commit (`fb9b90d`)

Raw: `data/models/corrected_foot_v1/corrected_foot_report.json` (every arm's
model sha256, the sha256 of the protocol module, and the scored result); per-step
traces under `data/derived/corrected-foot-arms/` (gitignored, regenerable).
Every number below is the scaffold's.

### Instrument checks: all nine pass

| | check | result |
|---|---|---|
| I1 | `moment_arms` called twice, bit-equal, every arm | pass (asserted in the runner) |
| I2 | `base` reads exactly 0 about mtp, all 8 | pass — the §0.2a known answer |
| I3 | `base_toe_paths` reproduces the recorded shipped-geometry arms | pass, ≤1e−6 m |
| I4 | `soleus_r` −0.0497 m about `ankle_angle_r`, every arm | pass, −0.049708 on all four |
| I5 | `corrected_fitted` still reads exactly 0 about mtp | pass |
| I6 | the corrected arms reproduce the offline prediction | pass, worst disagreement **0.0004 mm** |
| I7 | lumbar 42.69 / −52.82 / −62.79 mm (non-voiding) | pass: 42.6898 / −52.8238 / −62.7912 |
| I8 | `base` reproduces 0.122880 / 0.117805 / 1.449794 / 0.219587 | pass, all four |
| I9 | `base_toe_paths` reproduces ΔA +0.0784 / +0.0106 / −0.0574 / −0.0137 | pass, all four |

Not pre-registered, and reported because it happened: part A was run twice, in
two separate processes hours apart, and the two `moment_arms` dictionaries are
**identical for all four arms** — not just the repeat call inside one session.
(The second run happened because `--check`'s own digest had to be narrowed to the
files the builder writes; it was scanning the whole output directory and so began
failing once the runner's report landed there. The model bytes never changed —
`model.osim`'s sha256 is the same before and after — but the builder's own sha256
is recorded in each registration's `sources`, so the report was refreshed rather
than left quoting a builder that no longer exists.)

I6 is the load-bearing one. The engine and a fifteen-line piece of axis algebra
that never loads OpenSim agree on all eight corrected arms **to four ten-thousandths
of a millimetre**, and the prediction was written into
`scripts/run_corrected_foot_arms.py` and committed before the engine saw the model.

### A — the sign pattern is anatomically possible. GATE S1 PASSES.

Engine moment arms about `mtp_angle` at the rest pose, mm, right side (the left is
equal to 1e−14 on every arm):

| muscle | side | `base` | `base_toe_paths` | `corrected_fitted` | **`corrected_toe_paths`** | predicted |
|---|---|---:|---:|---:|---:|---:|
| edl | dorsal, extensor | 0 | **−4.780** | 0 | **−6.052** | −6.052 |
| ehl | dorsal, extensor | 0 | **+24.142** | 0 | **−7.551** | −7.551 |
| fdl | plantar, flexor | 0 | **+8.373** | 0 | **+6.834** | +6.834 |
| fhl | plantar, flexor | 0 | **+10.358** | 0 | **+7.129** | +7.129 |

| arm | extensors agree | flexors agree | extensors oppose flexors | min \|arm\| | verdict |
|---|---|---|---|---:|---|
| `base_toe_paths` | **no** (− vs +) | yes | — | 4.78 mm | **IMPOSSIBLE** |
| `corrected_toe_paths` | yes (both −) | yes (both +) | **yes** | 6.05 mm | **POSSIBLE** |

`base_toe_paths` reproduces `FOOT_JOINTS.md`'s impossible pattern to 1e−6 m, on a
rebuilt arm, which is the control that the two runs are measuring the same thing.
The correction turns it into two extensors together and two flexors together,
opposing, every arm well clear of the 1 mm floor. **The toe GeometryPaths now
encode toe flexion and extension.**

The magnitudes are Rajagopal's own, scaled: the generic model's arms at the same
construction are +5.64 / +6.59 / −6.47 / −5.79 mm and the corrected ones are
+6.05 / +7.55 / −6.83 / −7.13 in the opposite sign convention, i.e. 1.07–1.23× the
generic, which is the calcn scaling in those directions.

**It was the offset, not the axis, that carried the sign error.** The same algebra
on the two half-fixes: scaling the offset and leaving the axis pure-z gives
−5.27 / −14.49 / +8.81 / +4.51 mm, which is already a possible pattern; restoring
the oblique axis and leaving the offset unscaled gives +0.08 / +10.75 / +7.07 /
−2.33, which is still impossible and puts edl at 0.08 mm. Both halves are
restored here because both are equally derived from the parent, but the
attribution is the offset's.

**A consequence that could not have gone the other way, and the engine agrees.**
`corrected_toe_paths`'s ankle arms are **bit-identical** to `base_toe_paths`'s
(+39.52996182446783 mm for edl_r on both, and the same for the other seven).
Everything distal of the ankle rotates rigidly with it, so only the tibia→calcn
crossing segment sets the ankle arm, and moving the mtp axis cannot touch it.
`corrected_fitted`'s ankle arms are likewise bit-identical to `base`'s. So the
geometry correction moves the mtp arm and **nothing else in the moment-arm table.**

### B — what the toe hinge costs now

A = max(ankle_r, ankle_l) excursion past the declared range, rad. Reporting line
0.022 rad, declared before the data.

| protocol | `base` | `base_toe_paths` | `corrected_fitted` | **`corrected_toe_paths`** |
|---|---:|---:|---:|---:|
| supine unstopped | **0.1229** | 0.2012 (+0.0784 **over**) | 0.1260 (+0.0031) | **0.1094 (−0.0134)** |
| supine stopped | **0.1178** | 0.1284 (+0.0106) | 0.1182 (+0.0004) | **0.1136 (−0.0043)** |
| crawl 3 s unstopped | **1.4498** | 1.3924 (−0.0574 **over**) | 1.4435 (−0.0063) | **1.3738 (−0.0760 over)** |
| crawl 3 s stopped | **0.2196** | 0.2059 (−0.0137) | 0.2261 (+0.0065) | **0.2004 (−0.0192)** |

The 2×2, on A:

| protocol | G at P=0 | G at P=1 | P at G=0 | P at G=1 |
|---|---:|---:|---:|---:|
| supine unstopped | +0.0031 | −0.0918 | **+0.0784** | −0.0166 |
| supine stopped | +0.0004 | −0.0149 | +0.0106 | −0.0046 |
| crawl unstopped | −0.0063 | −0.0186 | **−0.0574** | **−0.0697** |
| crawl stopped | +0.0065 | −0.0055 | −0.0137 | **−0.0257** |

**The 27.8 mm of extra forefoot does almost nothing to the ankle.** G alone
(`corrected_fitted` − `base`) is +0.0031 / +0.0004 / −0.0063 / +0.0065 rad — every
one inside the reporting line, in a plant where `FOOT_JOINTS.md` Q1 found a 48.5 mm
plane shift worth more than a radian. The confound the pre-registration was built
around is real in principle and small in fact, and it is small **because it was
measured**, not because it was assumed. The supine plane sits at −0.403499 m on
every arm here, unchanged: the toes' contact spheres moved forward, not down, and
the plane is set by the torso ball.

**The toe arms' effect reverses sign and shrinks.** Installing toe GeometryPaths
on the *wrong* geometry pushed the supine ankle **+0.0784 rad**, over the line —
`FOOT_JOINTS.md`'s reason for listing option 2 as shifting the ankle "by up to
0.078 rad through changed extensor tension". On the corrected geometry the same
installation moves it **−0.0166 rad**, 4.7× smaller and the other way. The one
arm still over the line is the unstopped crawl, where the toe extensors' passive
tension reduces the ankle excursion from 1.4498 to 1.3738 rad (−5.2%). It remains
a 1.37 rad excursion past a declared range: **the corrected foot does not repair
the unstopped crawl's ankle, and nothing here claims it does.**

**mtp travel** (declared ±0.5236 rad):

| protocol | `base` mtp_r | **`corrected_toe_paths`** mtp_r | closest approach, corrected |
|---|---|---|---:|
| supine unstopped | −0.005 … +0.008 | **−0.029 … −0.004** | 0.494 rad |
| supine stopped | −0.004 … +0.008 | −0.034 … −0.004 | 0.490 |
| crawl unstopped | −0.027 … +0.092 | **−0.109 … +0.000** | 0.415 |
| crawl stopped | −0.144 … +0.076 | −0.073 … +0.040 | 0.449 |

The toes are now *driven*: in supine the whole range sits on the flexed side
instead of straddling zero, which is what a muscle with a real arm does to a joint
held otherwise by a −25·q spring. The travel is still nowhere near the bound in
any protocol, and the largest excursion of any arm is still a fifth of a radian
inside it.

**Crawl travel and divergence.** Pelvis forward travel over 3 s, stopped:
0.0729 (base) → 0.0762 (`corrected_fitted`) → **0.0781 m** (`corrected_toe_paths`),
+5.2 mm, +7.2%; unstopped 0.0678 → 0.0690 → 0.0702 m. Most of that is G, not P.
One trajectory per arm and no bar, so this is reported and not claimed. The first
divergence event is the same in every unstopped arm (`hip_rotation_l` past 0.35–0.36
rad), and the stopped crawl never diverges in any arm. Peak normalised fibre
velocity in the unstopped crawl falls 13.45 → 7.08 ofl/s, both under the 15 guard.

**Gates.** F2's bar, remeasured on each arm in the same protocol, is 0.2202
(`base`), 0.2191, 0.2199, 0.2199 — it moves by at most 0.0011 rad. **Neither F2
(`thoracic_extension` 0.3940) nor G-S (ankles 1.1399 / 0.8709) flips on any arm.**
Wall clock ranged 0.08–0.12 s per advance under shifting machine load; no
wall-clock difference is claimed.

### C — welding mtp in the articulated spine variant costs nothing measurable

`registration_mtp_welded.json`, F2's protocol, bar measured on the base model in
the same run: **0.220193** on `knee_angle_r`, reproducing F2's 0.2202, and `spine`
reproduces F2's `thoracic_extension` 0.3940.

| coordinate | `spine` | `spine_mtp_welded` | Δ |
|---|---:|---:|---:|
| `ankle_angle_r` | 1.3563 | 1.3562 | −0.00003 |
| `ankle_angle_l` | 1.3522 | 1.3513 | −0.00095 |
| `thoracic_extension` | 0.3940 | 0.3940 | −0.00001 |
| `pro_sup_l` | 0.2535 | 0.2535 | +0.00003 |
| `knee_angle_{l,r}` | 0.2249 / 0.2244 | 0.2249 / 0.2244 | <1e−5 |
| `subtalar_angle_r` | 0.1263 | 0.1284 | +0.00203 |
| `hip_rotation_r` | — | — | **−0.00606** (the largest of all 46) |

The same seven coordinates are over the bar on both plants: both ankles, both
knees, both `pro_sup`, `thoracic_extension`. The largest change the weld makes to
any coordinate is **0.0061 rad**, on `hip_rotation_r` — 2.8% of the bar and a
quarter of the reporting line.

Before the weld, `mtp_angle_{l,r}` in that plant came no closer than **0.448 rad**
to its ±0.5236 bound, held by the same passive spring as in the base plant and by
no muscle at all. Welding it removes two coordinates nothing could drive and costs
nothing measurable in this protocol, which is what `FOOT_JOINTS.md` predicted for a
new plant.

**GATE F2 IS NOT RESCORED.** F2 FAILED on
`data/models/articulated_spine_v1/registration.json` — `thoracic_extension` 0.3940
against a 0.2202 bar — and it is still FAILED. `model.osim` there is byte-identical
to the file F2 ran on (check W2). The welded model is a separate plant, its numbers
above are a new measurement of that plant, and they are not a new score for
anything.

## What this changes, and what it does not

* **`FOOT_JOINTS.md`'s option 2 — "give the toe muscles GeometryPaths" — is no
  longer "not viable as is".** Its two stated blockers are gone: the sign pattern
  is possible (S1), and the ankle shift it caused (+0.078 rad supine) is down to
  −0.017. It is viable **on the corrected geometry**, in `corrected_foot_v1`, and
  it remains not viable on `engineering_stance_v1`'s geometry, where nothing here
  has changed.
* **`FOOT_JOINTS.md`'s recommendation for `engineering_stance_v1` stands
  unchanged.** Leave it alone. Its `linearization.npz` was solved with mtp free at
  its shipped offset, and the corrected model has a foot 27.8 mm longer with
  different forefoot contact: adopting it there is a re-identification, not an
  edit. Nothing in this document is a reason to touch that plant.
* **The defect is upstream's, and it is still in `subject_walk_scaled.osim`.**
  Everything derived from that file — every plant in `data/models/` that is not
  `corrected_foot_v1` — carries a foot whose mtp axis is 27.8 mm proximal of where
  its own scaling put the muscles, and a metatarsal break that is a pure calcn-z
  hinge instead of Rajagopal's oblique one.
* **Not measured, so none of this reaches it:** upright stance and walking
  push-off, which is where toe loading is largest and where `WORKBENCH_AUTHENTICITY.md`
  §0.2a says a wrong toe would matter most. Also not measured: whether the longer
  forefoot is *better* against any external reference. It is what Rajagopal's
  proportions times this subject's own scale factors give, which is an argument
  about internal consistency and not about this subject's foot.

## Needs correcting elsewhere, and not edited here

* **`docs/WORKBENCH_AUTHENTICITY.md` §0.2a** says the geometry-path route is *"not
  viable as is"* and *"fix the foot geometry … before giving the toes geometry
  paths"*. The fix exists now, and §0.2a should point at this document. That file
  is outside this work's territory: **flagged, not edited.**
* **`docs/ARTICULATED_SPINE.md`** lists the variant's registrations
  (`registration_foot_paths.json`, `registration_muscled.json`) and does not yet
  know about `registration_mtp_welded.json`. Same: flagged, not edited.
* **`data/models/engineering_stance_v1` is deliberately NOT corrected**, and this
  is not an oversight to be tidied later. See the section above.
