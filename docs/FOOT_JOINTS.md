# the foot's joints: why the ankle collapses on the spine variant, and what the toe hinges cost

Two open questions, both about the foot, both carried by the SCAFFOLD (the
22/25-body plant), never by the body. Nothing below is a statement about a human
foot.

1. **The ankle.** On `articulated_spine_v1` the ankle leaves its declared range
   by 1.36 rad (kinematic variant) and still by 1.14 / 0.87 rad once the
   subtalar-crossing muscles have real arms (gate G-S, FAILED,
   `docs/ARTICULATED_SPINE.md`). The base plant leaves it by 0.12. So the
   missing subtalar arm is at most part of the cause, and the separating run —
   re-weld the subtalar alone — has not been made.
2. **The toes.** `engineering_stance_v1` keeps `mtp_{l,r}` as free PinJoints
   while its muscle paths were fitted with them welded, so edl/ehl/fdl/fhl have
   exactly 0 arm about `mtp_angle` (`docs/WORKBENCH_AUTHENTICITY.md` §0.2a).
   What does that cost, and what would each fix change?

    build arms   OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_foot_joint_arms --check
    run          OPENBLAS_NUM_THREADS=1 nohup nice -n 10 .venv/bin/python -u \
                     -m scripts.run_foot_joint_arms --part all > logs/foot_joint_arms.log 2>&1 &
    arms         data/derived/foot-joint-arms/<arm>/{model.osim,registration.json}  (gitignored,
                 deterministic; `--check` builds twice and compares bytes)

## PRE-REGISTRATION — written and committed before any arm below was integrated

Nothing in this section is edited after the runs. Results go in their own
section, in a separate commit.

### Q1 — the ankle

The variant differs from the base in three separable ways. Each is a factor:

| factor | what it is |
|---|---|
| **T** | torso repartitioned into torso / thorax / cervical / head, with the thoracic, neck and atlanto-occipital joints (9 coordinates, each with its stop and damping). In supine this ALSO changes the contact model: the engine gives every body a posterior ball inscribed in its inertia ellipsoid and puts the plane under the lowest ball, so T brings three new balls and a torso ball of 0.309 m instead of 0.258 m. |
| **W** | `radius_hand_{l,r}` WeldJoint → UniversalJoint (+ stops, damping) |
| **S** | `subtalar_{l,r}` WeldJoint → PinJoint (+ stop, damping) |

**Arms.** The full 2×2×2 factorial, named `t{T}w{W}s{S}`, plus two:

| arm | plant | role |
|---|---|---|
| `base`, `base_repeat` | `engineering_stance_v1/registration.json`, run twice | the bar; and *call it twice*: the two must agree exactly |
| `t0w0s0` | base → the variant builder's `unweld(S,W)` → `reweld(S,W)` | round-trip control: must equal `base` exactly, in the XML (checked statically) and in the engine |
| `t0w0s1` | base + subtalar only | is subtalar freedom *sufficient* on the base trunk? |
| `t0w1s0` | base + wrists only | are the wrists sufficient? (expected: no) |
| `t0w1s1` | base + wrists + subtalar | |
| `t1w0s0` | variant with wrists and subtalar re-welded | is the trunk partition sufficient? |
| `t1w0s1` | variant with wrists re-welded | |
| **`t1w1s0`** | **variant with ONLY the subtalar re-welded** | **the separating control `ARTICULATED_SPINE.md` names** |
| `t1w1s1` | `articulated_spine_v1/registration.json`, as committed | must reproduce F2's 1.3563 / 1.3522 |
| `tweld` | variant with all seven new joints re-welded | T's bodies and masses (and so its balls and plane) with none of its articulation: splits T into *mass/contact redistribution* vs *spine articulation* |
| `foot_paths` | `registration_foot_paths.json`, as committed | the G-S plant; must reproduce 1.1399 / 0.8709; used for the timing question |

T=0 arms are built from the base with the variant builder's own functions, so
every freed joint, stop and damping element is the variant's **element for
element** (checked). T=1 arms are the committed variant with joints re-welded
using the joint's own frames — the rest pose is the joint at q = 0 — and that
coordinate's stop and damping removed (checked: no reference to a welded
coordinate survives in the file). All arms run the shipped fitted path set, i.e.
**no muscle has a subtalar arm in any factorial arm** — the kinematic variant's
condition, on which the doc's regression was measured.

**Protocol**: F2 / G-S's, unchanged. Supine, 0.02 tonic excitation on every
muscle, 200 × 10 ms, worst excursion past each coordinate's declared range.
**Bar**: the base model's own worst excursion under the identical protocol,
measured in the same run (expected 0.2202 rad, `knee_angle_r`).

**Instrument checks — the run is void for comparison with prior numbers if any
fails** (internal comparisons within the run would still stand, and would be
reported as such):

- `base` == `base_repeat`, every worst excursion bit-equal;
- `t0w0s0` == `base`, bit-equal;
- the bar reproduces 0.2202, and `t1w1s1` and `foot_paths` reproduce F2's and
  G-S's printed ankle values, within 5e-5 (the rounding of a 4-dp print —
  not a physical tolerance).

**Primary reading.** An arm's ankle **recovers** iff both ankle excursions are ≤
the bar; otherwise it **collapses**. A = max(ankle_r, ankle_l) excursion.

| outcome | meaning |
|---|---|
| `t1w1s0` recovers | the collapse is caused by the subtalar's **freedom itself** (in this plant, with these contact proxies), not by the trunk or wrists. The G-S result then says muscle arms about the subtalar are not enough to control that freedom. |
| `t1w1s0` collapses | the cause is **elsewhere in the variant**; the subtalar is excluded as necessary. |
| `t0w0s1` collapses | subtalar freedom is **sufficient** on the base trunk. |
| `t1w0s0` collapses | the trunk partition is **sufficient**. Then: `tweld` collapses → it is the mass / contact redistribution (balls, plane), not the spine's articulation; `tweld` recovers → it is the spine's articulation. |
| `t0w1s0` collapses | the wrists are sufficient (not expected; reported if so). |
| collapse only where two factors are both 1 | an **interaction**; neither factor alone is the cause. |

Main effect of each factor = mean over its four pairs (arms differing only in
that factor) of ΔA; all eight A values are reported. The simulation is
deterministic, so there is no sampling error to clear — but it is also a single
nonlinear trajectory per arm, so a difference between two arms is a property of
those two trajectories and not of a population. That is why the reading is on
recover/collapse against the bar and on the *pattern* across eight arms, not on
the size of one ΔA.

**Secondary: does the collapse track the subtalar in time?** (the other
candidate the doc names — "the plantarflexors load a hinge"). For every arm with
a free subtalar and for `foot_paths`, per side: *t_A* = first time the ankle
excursion exceeds the bar; *t_S* = first time the subtalar reaches its own
declared bound; Pearson r between ankle and subtalar angle over the 2 s; the
direction of the ankle's extreme; plantarflexor and dorsiflexor tendon force and
foot contact force at t_A.

| outcome | meaning |
|---|---|
| t_S < t_A | the subtalar hit its bound before the ankle crossed the bar: **consistent with** a loaded hinge (not proof) |
| t_S ≥ t_A, or the subtalar never reaches its bound | the ankle collapses first: the hinge is **not the onset driver** |
| the collapse is in dorsiflexion | "plantarflexors" cannot be the driver of its direction |

This timing reading is subordinate to the factorial: if `t1w1s0` collapses, the
hinge is excluded as necessary whatever the timing says.

### Q2 — the toe hinges (base plant only; `engineering_stance_v1` is read, never written)

| arm | plant |
|---|---|
| `base` | `engineering_stance_v1` as identified |
| `base_mtp_welded` | `mtp_{l,r}` welded at q = 0, as upstream did before fitting (`exampleMocoInverse.cpp:51`). The two `mtp_angle` passive terms (−25·q − 2·qdot) are removed through a source override because the coordinate no longer exists. |
| `base_toe_geometry_paths` | model unchanged; shipped path set with edl/ehl/fdl/fhl (both sides) removed, so those 8 run on their own GeometryPaths — the route `1e29ee3` took for the subtalar. The other 72 paths are byte-equal in content (checked). |

Protocols (`docs/NATIVE_JOINT_LIMITS.md`'s two):

- **supine tonic**, as Q1, **unstopped** and **stopped** (crawl.py's stop set:
  every declared range except the ±10 rad sentinels, 30 N·m/rad, 85.94
  N·m·s/rad, 0.35 rad; the mtp-welded arm has no mtp stop because it has no mtp);
- **crawl seed**, `crawl.SEED`, prone, upright environment, 3 s — the stop sweep's
  horizon — **unstopped** and **stopped**. Runs the full horizon; where
  `crawl.diverged()` would have ended a search, the first such event is
  recorded, not acted on.

Measured on each: `mtp_angle` extrema and excursion past its declared ±0.5236;
ankle excursion; worst excursion over every other declared coordinate; for the
crawl, pelvis forward travel, first divergence event, fibre-velocity peak; wall
clock. Plus engine moment arms of the 8 toe muscles about ankle and mtp on
`base` and `base_toe_geometry_paths`, each called twice (must be equal):
`base` must read exactly 0 about mtp (the §0.2a finding, a known answer), and
the geometry arm must read nonzero (the control that the fix is installed).

**Reading.** No gate is set on mtp travel — it is the measurement asked for.
"Changes the ankle" is reported when |ΔA| against `base` exceeds **0.022 rad**:
10% of F2's bar, a *reporting convention with no measured provenance*, stated
here so it cannot be chosen after the fact. "Changes a gate" means: measured on
that arm, would F2's bar move enough to flip F2 (thoracic_extension 0.3940) or
G-S (ankles 1.1399 / 0.8709), and in the crawl, does the first-divergence event
appear, disappear, or change coordinate.

**Not measured, and not claimed**: upright stance, walking, the identified
linearization. `linearization.npz` was solved with mtp free, so a plant with mtp
welded is not the identified plant, whatever these numbers say.

## Results — run 18 Sep 2026, after the pre-registration commit (`e00844d`)

Raw: `data/models/articulated_spine_v1/foot_joints_report.json` (scored results,
every arm's model sha256); per-step traces under
`data/derived/foot-joint-arms/results/` (gitignored, regenerable). Every number
below is the scaffold's.

### Instrument checks: all seven pass

| check | result |
|---|---|
| `base` run twice, every worst excursion bit-equal | pass |
| round trip `t0w0s0` == `base`, bit-equal (and element-equal in the XML) | pass |
| bar reproduces F2's 0.2202 | 0.220193, `knee_angle_r` |
| `t1w1s1` reproduces 1.3563 / 1.3522 | pass |
| `foot_paths` reproduces G-S's 1.1399 / 0.8709 | pass |

Not pre-registered, but also reproduced: the crawl seed on the base plant, 3 s,
gives `ankle_angle_r` 1.4498 unstopped and 0.2196 stopped. The stop sweep in
`NATIVE_JOINT_LIMITS.md` recorded 1.450 and 0.220.

### Q1 — the ankle collapse is caused by the trunk partition, not the subtalar

| arm | T | W | S | ankle r | ankle l | verdict | supine plane x (m) |
|---|:-:|:-:|:-:|---:|---:|---|---:|
| base | – | – | – | 0.1229 | 0.1207 | recovered | −0.4035 |
| t0w0s0 | 0 | 0 | 0 | 0.1229 | 0.1207 | recovered | −0.4035 |
| t0w0s1 | 0 | 0 | 1 | 0.1129 | 0.1080 | recovered | −0.4035 |
| t0w1s0 | 0 | 1 | 0 | 0.1204 | 0.1182 | recovered | −0.4035 |
| t0w1s1 | 0 | 1 | 1 | 0.1102 | 0.1052 | recovered | −0.4035 |
| t1w0s0 | 1 | 0 | 0 | **1.3851** | **1.2663** | **COLLAPSES** | −0.4520 |
| t1w0s1 | 1 | 0 | 1 | **1.3563** | **1.3522** | **COLLAPSES** | −0.4520 |
| **t1w1s0** | 1 | 1 | 0 | **1.3847** | **1.2665** | **COLLAPSES** | −0.4520 |
| t1w1s1 | 1 | 1 | 1 | **1.3563** | **1.3522** | **COLLAPSES** | −0.4520 |
| tweld | welded | 0 | 0 | **1.4014** | **1.2876** | **COLLAPSES** | −0.4520 |
| foot_paths | 1 | 1 | 1 + arms | 1.1399 | 0.8709 | COLLAPSES | −0.4520 |

Bar 0.2202. Main effects on A = max(ankle_r, ankle_l): **T +1.254 rad**, W −0.001,
S −0.019.

The pre-registered readings:

* **The separating control `t1w1s0` does NOT recover.** Re-welding the subtalar
  alone leaves the ankle at 1.385 / 1.267. **The subtalar is excluded as a
  necessary cause.**
* **Subtalar freedom is not sufficient either.** `t0w0s1` stays at 0.113 / 0.108,
  inside the bar and slightly *better* than the base.
* **The trunk partition is sufficient.** `t1w0s0` collapses. T is present in all
  four collapsed factorial arms, and W and S are not needed for any of them.
* **`tweld` collapses** (1.401 / 1.288). It has all of T's bodies and masses and
  none of its articulation, so the cause is **T's mass/contact redistribution,
  not the spine joints**. The wrists do nothing.
* **Timing.** In the S=1 collapsed arms of the factorial, the ankle crosses the
  bar at 0.34 s and the subtalar reaches its bound at 0.35/0.36 s: **the ankle
  goes first, so the hinge is not the onset driver.** In the `foot_paths` plant
  the subtalar reaches its bound (now in eversion, +0.80) at 0.28 s, before the
  ankle crosses at 0.30 s. That is "consistent with a loaded hinge", but the
  pre-registration makes timing subordinate to the factorial, and the factorial
  has excluded the hinge as necessary.
* **The collapse is not the plantarflexors.** At the moment of crossing, in
  every collapsed arm, the plantarflexors carry **25–31 N** and the dorsiflexors
  **468–476 N**. The ankle is being driven into plantarflexion *against* its own
  dorsiflexors. In the base at the same time it is 171 N against 221 N.

**What T changes that could do this, measured.** The supine environment gives
every body a posterior ball inscribed in its inertia ellipsoid and puts the plane
under the lowest ball (`native_mechanical_stream.cpp:206-217`). Repartitioning the
torso took the head and rib cage out of it, and its ball grew from 0.2577 to
0.3090 m. That lowered the plane by **48.5 mm**. The arithmetic closes: +51.3 mm
of radius against a 2.8 mm COM shift, and the engine reports −0.4035 → −0.4520.
In every T arm the legs therefore fall 48.5 mm further before the talus ball
lands. The trace shows the landing:

| | talus contact peak (at 0.20 s) | ankle at 0.40 s | toes on the plane | final ankle |
|---|---:|---:|---|---:|
| base | 363 N | −0.888 | never | −0.640 |
| t1w0s0 | **1,142 N** | −1.364 | from 0.45 s | −2.258 |
| tweld | **1,026 N** | −1.222 | from 0.48 s | −2.243 |

A 3.1× larger heel impact comes 0.14 s before the ankle crosses the bar. After
that the foot folds until the toes rest on the plane, and it stays there.

**Not separated, and not claimed:** whether the plane drop (contact) or the mass
redistribution itself drives the fold. `tweld` changes both, and the plane
position is computed inside the engine, which this work may not edit. The run that
would separate them needs an engine option for the supine plane offset: `t0w0s0`
with the plane lowered 48.5 mm, and `tweld` with it raised 48.5 mm.

**What this does to earlier claims.**

* `ARTICULATED_SPINE.md`'s "un-welding the subtalar makes the ankle far worse"
  is **withdrawn**. The subtalar is the one change that makes no difference. The
  ankle regression arrived with the torso partition. The supine contact proxy is
  the leading candidate mechanism, but it is not yet separated from the mass
  change.
* **G-S's premise is withdrawn; its verdict stands.** G-S asked whether subtalar
  arms repair the ankle, and it FAILED. That verdict is kept, but the collapse it
  tested was never a subtalar collapse. The partial improvement it saw (16% and
  36%) is a property of the foot-paths plant and is not explained here.
* **F2's bar compared plants whose support planes differ by 48.5 mm.** The base
  lies on a plane set by a 0.258 m torso ball and the variant on one set by a
  0.309 m ball. The confound F2 names for `thoracic_extension` (new balls on the
  plane) is real, and it reaches the feet as well. F2 stays FAILED. It is not
  rescored.
* **This is a scaffold artefact about ball placement, not a finding about any
  foot.** Supine results on any plant that repartitions a body have the same
  exposure, because a body's ball, and therefore the plane, follows its inertia.

### Q2 — the toe hinges cost almost nothing measurable in these protocols

**Moment arms at rest** (engine, each call made twice, bit-equal):

| muscle (r; l identical) | ankle, base | ankle, geometry | **mtp, base** | **mtp, geometry** |
|---|---:|---:|---:|---:|
| edl | +39.3 mm | +39.5 | **0** | −4.8 |
| ehl | +42.6 | +42.5 | **0** | +24.1 |
| fdl | −11.3 | −11.5 | **0** | +8.4 |
| fhl | −18.0 | −18.2 | **0** | +10.4 |

The §0.2a known answer reproduces: exactly 0 about mtp on the base. The geometry
arms have the **wrong pattern**. ehl, which runs on the dorsal side, and fdl/fhl,
which run on the plantar side, share one sign, and edl is opposite to ehl. That
cannot be right for two extensors and two flexors under any sign convention. The
cause is in the file: in this plant **and in upstream `subject_walk_scaled.osim`**,
the toe muscles' last calcaneal via points sit **8–22 mm distal to the mtp axis**
they are meant to straddle. In `RajagopalLaiUhlrich2023.osim` they sit 5–17 mm
proximal to it. The scaling multiplied the calcn path points by 1.156 (and the
subtalar offset by 1.228), but left the mtp offset at Rajagopal's unscaled
0.1788 m. The model's mtp axis is also pure calcn-z, where Rajagopal's is oblique.
**So the toe GeometryPaths in this plant do not encode toe flexion and extension.**
Installing them "gives the toes muscles" in name only.

**mtp travel** (declared ±0.5236 rad):

| protocol | base mtp_r | base mtp_l | closest approach to the bound |
|---|---|---|---:|
| supine tonic, unstopped | −0.005 … +0.008 | −0.005 … +0.008 | 0.515 rad |
| supine tonic, stopped | −0.004 … +0.008 | −0.004 … +0.008 | 0.515 |
| crawl seed 3 s, unstopped | −0.027 … +0.092 | −0.005 … +0.092 | 0.432 |
| crawl seed 3 s, stopped | **−0.144** … +0.076 | −0.104 … +0.077 | 0.379 |

The passive term (−25·q − 2·qdot) holds the toes within 0.15 rad of neutral in
every protocol measured.

**Consequences of each option** (ΔA against the base in the same protocol;
reporting line 0.022 rad, declared before the data):

| protocol | weld mtp: ΔA | weld: travel | toe geometry paths: ΔA | geometry: travel |
|---|---:|---:|---:|---:|
| supine unstopped | −0.0001 | – | **+0.0784** (over the line) | – |
| supine stopped | 0.0000 | – | +0.0106 | – |
| crawl unstopped | +0.0008 | 67.8 → 68.0 mm | **−0.0574** (over the line) | 67.8 → 67.8 mm |
| crawl stopped | +0.0016 | 72.9 → 73.7 mm | −0.0137 | 72.9 → 73.7 mm |

Gates: welding moves F2's bar by 2e-5 rad and the geometry paths move it by
−0.0011. **Neither flips F2 or G-S.** In the unstopped crawl the first-divergence
event is the same in every arm (`hip_rotation_l` past 0.35 rad at 0.42 / 0.42 /
0.43 s). The stopped crawl never diverges in any arm.

The geometry paths change the ankle through **passive tendon force, not through
mtp**. At 0.5 s of the supine run, with the ankle near −0.8 rad, edl carries 32 N
and ehl 11 N against the fitted set's 64 N and 42 N. The toe extensors resist
plantarflexion less, so the ankle goes further. Why the geometric and fitted
lengths diverge at that angle (plausibly the fit's domain, which is walking) was
not isolated.

Wall clock: the machine was under shifting load, and the same plant ranged from
0.10 to 1.03 s per advance between runs. No wall-clock difference is claimed.

**The options, with what they were measured to cost:**

1. **Weld mtp as upstream did.** No measured consequence to the ankle, to any
   gate, or to crawl travel (+0.2 / +0.8 mm) in these four protocols. It matches
   how the paths were fitted, and it removes two coordinates no muscle can drive.
   *Unmeasured cost:* the identified stance has `mtp_angle` = 0.113 rad in
   `initial_pose.json`, and a weld at 0 changes that posture. `linearization.npz`
   was solved with mtp free. **In `engineering_stance_v1` it would therefore be a
   re-identification, not an edit.**
2. **Give the toe muscles GeometryPaths (the 1e29ee3 route).** **Not viable as
   is.** The arms it installs have the wrong pattern, because the via points
   were scaled and the mtp offset was not. It also shifts the ankle by up to
   0.078 rad through changed extensor tension. It needs the foot geometry
   reconciled first: the mtp offset scaled with the calcn, and ideally
   Rajagopal's oblique mtp axis.
3. **Leave it.** Measured cost in these protocols: the toes stay within 0.15 rad
   of neutral, held by a passive spring. They behave as nearly welded, and
   welding them changes nothing measured.

**Recommendation from this evidence:** leave `engineering_stance_v1` as it is,
because its cost is below everything measured here. Weld mtp in any NEW plant,
where it costs nothing and removes a dishonest degree of freedom. Do not install
toe GeometryPaths until the mtp axis and offset are fixed.

*Not measured, so this recommendation does not reach it:* upright stance and
walking push-off. That is where toe loading is largest, and where "toes that no
muscle could move" (§0.2a) would matter most.
