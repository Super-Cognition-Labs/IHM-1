# the foot's joints: why the ankle collapses on the spine variant, and what the toe hinges cost

Three questions, all about the foot, all carried by the SCAFFOLD (the 22/25-body
plant), never by the body. Nothing below is a statement about a human foot.

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
3. **The plane.** Q1 below found the ankle collapse is carried by the TRUNK
   REPARTITION, and named the supine contact proxy's plane as the leading
   candidate mechanism — but could not separate the plane's 48.5 mm drop from the
   mass change that caused it, because the engine computed the plane. The engine
   now takes an option for it. Q3 is the separating run.

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

## PRE-REGISTRATION — Q3, the plane against the mass (written and committed before any arm below ran)

Q1 left one thing unseparated and said so: *"whether the plane drop (contact) or
the mass redistribution itself drives the fold. `tweld` changes both, and the
plane position is computed inside the engine, which this work may not edit."*

**The engine option now exists.** `NativeMechanicalStream(...,
support_plane_source_x_m=<float>)` PINS the supine plane instead of hanging it
under the lowest inertia-inscribed proxy sphere; it is refused in the upright
environment, pinning it at the default reproduces the default exactly, and
`execution.json` records `support_plane_override_x_m` and a basis string
(`scripts/verify_engine_options.py`). Nothing else about any plant changes.

Nothing in this section is edited after the runs. Results go in their own
section, in a separate commit.

    run   OPENBLAS_NUM_THREADS=1 nohup nice -n 10 .venv/bin/python -u \
              -m scripts.run_plane_arms > logs/plane_arms.log 2>&1 &

**Protocol, bar and scoring are Q1's, unchanged**: supine, 0.02 tonic excitation
on every muscle, 200 x 10 ms, worst excursion past each coordinate's declared
range; the bar is the base model's own worst under the identical protocol
measured in the same run (expected 0.2202 rad, `knee_angle_r`); an arm
**recovers** iff BOTH ankle excursions are <= the bar, otherwise it
**collapses**; A = max(ankle_r, ankle_l).

`P0` is the base plant's OWN derived plane and `PT` the repartitioned trunk's,
each taken at full precision from the engine in this run, not from the 5-dp
prints (expected -0.40350 and -0.45201; if they do not reproduce, the run is
void against the prior numbers).

### Arms

Plane DERIVED, exactly as every Q1 arm ran — these carry the bar and the
instrument checks, and they are re-run because **the engine has been rebuilt
since Q1** (`build-cx6s8y89` -> `build-8crm1_k9`) and a rebuilt engine is a
changed instrument until it reprints the same numbers:

| arm | plant | role |
|---|---|---|
| `base`, `base_repeat` | `engineering_stance_v1` | the bar; and *call it twice* |
| `t0w0s0` | base -> unweld(S,W) -> reweld(S,W) | round-trip control, must equal `base` |
| `t1w1s1` | the committed variant | must reproduce 1.3563 / 1.3522 |
| `foot_paths` | `registration_foot_paths.json` | must reproduce 1.1399 / 0.8709 |
| `tweld` | T's bodies and masses, all seven new joints welded | the repartition with no articulation; must reproduce 1.4014 / 1.2876 |
| `t1w0s0` | T alone in the factorial | must reproduce 1.3851 / 1.2663 |

Plane PINNED — the new arms:

| arm | plant | plane | role |
|---|---|---|---|
| `base@P0` | base | pinned at its own `P0` | **control**: pinning at the default must be a no-op, bit-equal to `base` |
| `tweld@PT` | `tweld` | pinned at its own `PT` | **control**: bit-equal to `tweld` |
| **`tweld@P0`** | `tweld` | **raised 48.5 mm to the base's** | **primary A** — repartitioned trunk on the base plant's floor: does the ankle recover? |
| **`base@PT`** | base | **lowered 48.5 mm to the variant's** | **primary B** — base trunk, no mass change at all: does the ankle collapse? |
| `t1w0s0@P0` | `t1w0s0` | raised to `P0` | secondary: is A a property of `tweld` alone, or of the repartition? |
| `base@P0-d` | base | `P0` lowered by d = 12.125, 24.25, 36.375 mm | the **ladder**: does the fold track plane depth continuously? d = 0 is `base@P0` and d = 48.5 is `base@PT`, so the ladder is five points |

### A hazard declared before the run, and it is asymmetric

Pinning `tweld`'s plane UP to `P0` puts its 0.309 m torso proxy sphere **48.5 mm
inside the plane at t = 0**. Lowering the base plant's plane to `PT` starts it
48.5 mm **above** the floor, which is a clean free fall. So **primary B is clean
and primary A carries an initial-penetration artefact**, and there is no arm that
raises a floor without one, because the plane is by construction the lowest
sphere's tangent.

Therefore: every arm reports the signed gap from every proxy sphere to the plane
at t = 0, reconstructed from the engine's own emitted mass properties
(`radius^2 = 5(I1 + I2 - I0)/2m`), and the contact force on every traced body at
t = 0 and at its peak. If `tweld@P0`'s torso contact at t = 0 is large, its
verdict is **weaker evidence than B's**, and the LADDER — every point of which is
a clean drop with no penetration — carries the dose reading. This is said now so
it cannot be said afterwards only if it is convenient.

### Instrument checks — the run is void against Q1's numbers if any of 1-7 fails

Internal comparisons within this run would still stand and would be reported as
such.

1. `base` == `base_repeat`, every worst excursion bit-equal
2. `t0w0s0` == `base`, bit-equal
3. the bar reproduces F2's 0.2202 (within 5e-5, the rounding of a 4-dp print)
4. / 5. `t1w1s1` reproduces 1.3563 / 1.3522
6. / 7. `foot_paths` reproduces G-S's 1.1399 / 0.8709

Three more, new, and **each can fail for the reason it exists**:

8. `base@P0` bit-equal to `base` AND `tweld@PT` bit-equal to `tweld` — pinning at
   the default is a no-op. (A pin that silently did nothing at all would also pass
   this, which is why 9 is separate.)
9. the engine's reported `support_plane_source_x_m` equals the pinned value
   **exactly** in every pinned arm, and the arm's own `execution.json` carries
   `support_plane_override_x_m` and the pinned basis string — read back out of the
   run's own record, never inferred from the flag passed. A pinned arm whose plane
   differs from an unpinned one by 48.5 mm proves the pin bites.
10. the plane rule reconstructed from the engine's own snapshot: for every
    DERIVED-plane arm, `min over bodies of (com_ground_x - r) - plane == 0` to
    1e-12, the lowest body is `torso`, and its radius reads 0.2577 m on the base
    and 0.3090 m under the repartition. This is the check that the reconstruction
    used to report gaps is the engine's own rule and not a story about it.

Reported separately, not among the seven: `tweld` and `t1w0s0` reproducing
1.4014 / 1.2876 and 1.3851 / 1.2663.

### What each outcome would mean — all four, including the one where neither explains it

| `tweld@P0` (repartitioned mass, floor HELD up) | `base@PT` (base mass, floor LOWERED) | reading |
|---|---|---|
| **recovers** | **collapses** | **THE PLANE.** The 48.5 mm drop is sufficient on an unchanged mass and necessary under the repartition. The collapse is the contact proxy's plane placement, and `FOOT_JOINTS.md`'s leading candidate is confirmed. |
| **collapses** | **recovers** | **THE MASS.** The plane drop is neither sufficient nor necessary. The leading candidate is **withdrawn**: something in the repartitioned inertia does it, and the next question is what. |
| **recovers** | **recovers** | **NEITHER ALONE.** The collapse needs both terms together; each is necessary and neither is sufficient. The single-factor story is wrong and the candidate is downgraded to "one of two jointly necessary terms", not confirmed. |
| **collapses** | **collapses** | **OVER-DETERMINED, and this experiment cannot assign it.** The plane drop alone collapses the base plant AND holding the floor up does not rescue the repartition. Reported as two independently sufficient routes, or as a third cause common to both; not as a separation. |

And the fifth outcome, which is not in the table because it voids it: **if any of
checks 8-10 fails the option is not doing what it says, and nothing here is
reported as a separation.** Recorded as VOID.

**The ladder's reading**, A against plane depth on the base plant at
d = 0, 12.125, 24.25, 36.375, 48.5 mm: monotone and graded -> the fold tracks fall
height continuously; flat then a jump -> a threshold, and where; non-monotone ->
reported as non-monotone, with nothing fitted to it.

### What is NOT separated by this run, whatever it returns

* The plane's **position** is separated from the mass. The proxy **spheres** are
  not: `tweld@P0` still carries the 0.309 m torso ball and the three new balls,
  just against a higher floor. "The plane" below always means its position.
* A pinned plane is an **instrument**, not a better model of a bed. It is no more
  anatomical than a derived one.
* G-S and F2 stay **FAILED** and are not rescored, whatever this finds.
* Nothing here is a statement about a human foot, ankle or trunk. It is the
  scaffold throughout.

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

> **SEPARATED, Q3 below (18 Sep 2026): it is NOT the plane.** Holding the floor up
> under the repartitioned trunk leaves the ankle at 1.3961 / 1.2820; dropping the
> floor 48.5 mm under the base trunk leaves it at 0.1532 / 0.1504, inside the bar.
> The plane's position is worth ≤0.03 rad in either direction. **The "leading
> candidate mechanism" named in this section and in the paragraph above it is
> withdrawn**; the factorial it sits under is unaffected. The talus-impact table
> below reproduces, but `base@PT` reproduces that impact (992 N) without the fold,
> so it is a co-symptom.

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
2. **Give the toe muscles GeometryPaths (the 1e29ee3 route).** ~~**Not viable as
   is.**~~ The arms it installs have the wrong pattern, because the via points
   were scaled and the mtp offset was not. It also shifts the ankle by up to
   0.078 rad through changed extensor tension. It needs the foot geometry
   reconciled first: the mtp offset scaled with the calcn, and ideally
   Rajagopal's oblique mtp axis. **DONE, 18 Sep 2026 —
   `docs/FOOT_GEOMETRY.md`.** The mtp calcaneal offset is the ONE joint offset in
   `subject_walk_scaled.osim` the scaling step missed (19 of 19 comparable
   non-zero offsets obey the rule; this one is the unscaled generic value), and
   the same two joints are the only ones whose orientation was dropped. Scaled to
   0.20662282894439998 m and given Rajagopal's oblique axis back, in the new
   variant `data/models/corrected_foot_v1`, the pattern becomes **possible** —
   edl −6.05, ehl −7.55, fdl +6.83, fhl +7.13 mm, extensors together against
   flexors together — and the supine ankle shift drops from +0.078 to −0.017 rad.
   **The route is viable on the corrected geometry**, and still not viable on this
   plant's, which is unchanged.
3. **Leave it.** Measured cost in these protocols: the toes stay within 0.15 rad
   of neutral, held by a passive spring. They behave as nearly welded, and
   welding them changes nothing measured.

**Recommendation from this evidence:** leave `engineering_stance_v1` as it is,
because its cost is below everything measured here. Weld mtp in any NEW plant,
where it costs nothing and removes a dishonest degree of freedom. Do not install
toe GeometryPaths until the mtp axis and offset are fixed. *(Both amended 18 Sep
2026, `docs/FOOT_GEOMETRY.md`: mtp is now welded in `articulated_spine_v1` as a
separate registration — largest change to any of its 46 coordinates 0.0061 rad,
F2 not rescored and still FAILED — and the mtp axis and offset are fixed in
`corrected_foot_v1`. The recommendation to leave `engineering_stance_v1` alone is
unchanged and is reinforced: the corrected foot is 27.8 mm longer at the
forefoot, so adopting it there is a re-identification, not an edit.)*

*Not measured, so this recommendation does not reach it:* upright stance and
walking push-off. That is where toe loading is largest, and where "toes that no
muscle could move" (§0.2a) would matter most.

## Q3 Results — run 18 Sep 2026, after the pre-registration commit (`e8d25d8`)

Raw: `data/models/articulated_spine_v1/plane_arms_report.json` (every arm's plane,
override, basis string, worst excursions and reconstructed proxy spheres);
scored copy and per-step traces under `data/derived/foot-joint-arms/results/`
(gitignored, regenerable). Engine build `build-8crm1_k9`; tree at `8e8f8c6`,
which has the pre-registration `e8d25d8` as an ancestor and differs from it only
in another worker's bed-catalogue fix. 15 arms, 268 s. Every number is the
scaffold's.

### Instrument checks: all sixteen pass

Q1's seven, re-run because the engine was rebuilt after Q1:

| check | result |
|---|---|
| `base` run twice, every worst excursion bit-equal | pass |
| round trip `t0w0s0` == `base`, bit-equal | pass |
| bar reproduces F2's 0.2202 | 0.220193, `knee_angle_r` |
| `t1w1s1` reproduces 1.3563 / 1.3522 | pass |
| `foot_paths` reproduces G-S's 1.1399 / 0.8709 | pass |

**So the rebuilt engine reprints every Q1 number**, and `tweld` and `t1w0s0`
reprint 1.4014 / 1.2876 and 1.3851 / 1.2663 as well (reported separately, as
pre-registered). The two planes reproduce at full precision:
**P0 = −0.403498987616 m**, **PT = −0.452012769824 m**, drop **48.5138 mm**.

The three new ones, each of which could have failed for the reason it exists:

| check | result |
|---|---|
| 8a `base@P0` bit-equal to `base`; 8b `tweld@PT` bit-equal to `tweld` | pass — pinning at the default is a no-op |
| 9a-c every pinned arm reports the pinned plane exactly, records `support_plane_override_x_m` and the PINNED basis in its own `execution.json`, and no derived arm records one | pass |
| 9d the pin bites: `tweld@P0`'s plane is 48.5138 mm from `tweld`'s | pass |
| 10a-c the plane rule reconstructed from the engine's own emitted mass properties puts the plane exactly under the lowest sphere (residual **0.0**, all seven derived arms), the lowest body is `torso` every time, and the torso radius reads **0.2577 m** (base) and **0.3090 m** (repartitioned) | pass |

Nothing is void. The separation stands.

### The answer: the MASS, not the plane

| arm | plant | plane (m) | ankle r | ankle l | verdict |
|---|---|---:|---:|---:|---|
| `base` | base | −0.403499 | 0.1229 | 0.1207 | recovered |
| `base@P0` | base | −0.403499 *pinned* | 0.1229 | 0.1207 | recovered (control) |
| **`base@PT`** | **base** | **−0.452013 *pinned*** | **0.1532** | **0.1504** | **RECOVERED** |
| `tweld` | repartitioned | −0.452013 | 1.4014 | 1.2876 | COLLAPSES |
| `tweld@PT` | repartitioned | −0.452013 *pinned* | 1.4014 | 1.2876 | COLLAPSES (control) |
| **`tweld@P0`** | **repartitioned** | **−0.403499 *pinned*** | **1.3961** | **1.2820** | **COLLAPSES** |
| `t1w0s0` | T alone | −0.452013 | 1.3851 | 1.2663 | COLLAPSES |
| `t1w0s0@P0` | T alone | −0.403499 *pinned* | 1.3846 | 1.2610 | COLLAPSES |

Bar 0.2202, as in Q1.

The pre-registered quadrant is **(A collapses, B recovers) → THE MASS**:

* **Holding the floor up does not rescue the repartition.** `tweld@P0` gives
  **1.3961 / 1.2820** against `tweld`'s 1.4014 / 1.2876 — moving the floor 48.5 mm
  is worth **−0.0054 rad** on a 1.40 rad collapse. `t1w0s0@P0` repeats it at
  **−0.0005 rad**.
* **Dropping the floor does not collapse the base plant.** `base@PT` gives
  **0.1532 / 0.1504**, +0.0303 rad, comfortably inside the bar.
* **The trunk repartition is worth the same at either floor**: +1.2732 rad
  measured at P0 and +1.2482 rad at PT, against Q1's main effect of +1.254 rad
  measured with the floor free to move. The factor and the floor are not
  interacting.

**`docs/FOOT_JOINTS.md`'s leading candidate — "the plane drop is the mechanism" —
is WITHDRAWN.** The plane's position is worth at most 0.03 rad in either plant
and in either direction. The collapse is carried by what the repartition does
besides moving the floor.

### The ladder: graded, monotone, and two orders of magnitude too small

Base plant, floor lowered in equal steps, nothing else changed. Every point is a
clean drop with no penetration.

| drop below P0 | plane (m) | ankle r | ankle l | A | verdict |
|---:|---:|---:|---:|---:|---|
| 0 mm | −0.403499 | 0.1229 | 0.1207 | 0.1229 | recovered |
| 12.125 | −0.415624 | 0.1306 | 0.1285 | 0.1306 | recovered |
| 24.25 | −0.427749 | 0.1392 | 0.1367 | 0.1392 | recovered |
| 36.375 | −0.439874 | 0.1461 | 0.1434 | 0.1461 | recovered |
| 48.5 | −0.452013 | 0.1532 | 0.1504 | 0.1532 | recovered |

**Monotone increasing, no threshold, no collapse at any depth.** The fold tracks
plane depth continuously at **6.26e-4 rad/mm**, with a span of 0.0303 rad over the
whole 48.5 mm. At that slope the floor would have to drop **156 mm** to put the
ankle on F2's bar and about **2.0 m** to reach the 1.40 rad the repartition
produces — linear extrapolation that far outside the measured range is not
evidence, but the order is the point.

### Q1's proposed mechanism is falsified too, by the impact it named

Q1 wrote: *"the legs therefore fall 48.5 mm further before the talus ball lands
… a 3.1× larger heel impact comes 0.14 s before the ankle crosses the bar."*
The pinned arms reproduce that impact **without** the fold, and reproduce the
fold **without** the extra fall. The initial gap from each foot proxy to the
plane is identical between the arms that share a plane — `base` and `tweld@P0`
both start with the talus ball 173.70 mm above the floor, `base@PT` and `tweld`
both at 222.21 mm — so the fall geometry is matched exactly:

| arm | talus ball starts | talus contact peak | ankle A |
|---|---:|---:|---:|
| `base` | 173.70 mm up | 363 / 374 N | 0.1229 |
| `tweld@P0` | **173.70 mm up (same fall)** | 815 / 779 N | **1.3961** |
| `base@PT` | 222.21 mm up (48.5 mm further) | **992 / 925 N** | **0.1532** |
| `tweld` | 222.21 mm up | 1,026 / 1,066 N | 1.4014 |

`base@PT` lands its heels **2.7× harder than `base`** — as hard as the collapsed
arms — and its ankle moves 0.03 rad. `tweld@P0` lands them *softer* than `tweld`
and folds anyway. **Heel impact magnitude is not the driver.** Neither is the
fall distance. The ankle's direction and the muscle picture are unchanged from
Q1's: plantarflexion to −2.27 rad, and at the crossing 33 N of plantarflexor
against 468 N of dorsiflexor in `tweld@P0`, against 171 N / 221 N in the base.

### The declared hazard, measured

It behaved as declared and it did not decide anything. `tweld@P0` and
`t1w0s0@P0` start with the torso proxy **48.51 mm inside the floor** and
**2,809 N** of torso contact at t = 0 — by far the largest force anywhere in the
run. Their ankles still land within **0.0054 rad** of their own unpinned
versions. So the penetration transient is a large force on the trunk that does
essentially nothing to the ankle, and the arm that carries no penetration at all
(`base@PT`, and every rung of the ladder) gives the same verdict from the other
side. Both directions agree; the asymmetry did not matter.

*Also recorded, because it is over a bar:* on `base@PT` the **worst** coordinate
is `knee_angle_l` at **0.2204**, 0.0002 rad above F2's 0.2202 — and the ladder's
worst coordinate is **not** monotone (0.2202, 0.2194, 0.2143, 0.2226, 0.2204)
even though the ankle is. The ankle verdict is defined on the ankles and is
unaffected. **F2 is not rescored and stays FAILED**; this is noted because a bar
measured on a plant whose floor can move is a bar with a term in it, which is
what F2's own confound note already says.

### What this leaves, and it is still the contact proxy — just not the plane

**Exonerated:** the plane's *position*. It is a graded, small term.

**Not exonerated, and not separated by this run:** the proxy **sphere** itself.
Pinning the plane does not pin the ball. The repartitioned torso's ball is
0.3090 m against the base's 0.2577 m at every arm, so wherever the floor is, that
trunk's contact point sits **51.3 mm further from its own mass centre** — the
trunk is propped that much higher off the surface relative to the hips and heels.
That geometric difference is present in `tweld@P0` exactly as it is in `tweld`,
which is consistent with the two scoring the same, and it is a property of the
*inertia-inscribed sphere*, not of any anatomy.

So the surviving candidates are two, and this run does not choose between them:

1. the **proxy sphere radius** — a contact artefact, `r² = 5(I₁+I₂−I₀)/2m` on a
   body whose inertia the repartition changed;
2. the **repartitioned inertia itself**, acting through the dynamics rather than
   through contact.

**What would separate them:** an engine option for the per-body proxy radius, the
exact analogue of the plane option used here — run `tweld` with the torso ball
pinned at the base's 0.2577 m. That is the next arm, and it is named here so it
is not re-derived. Until it is run, "the contact proxy is the leading candidate"
remains true of the *sphere*; it is false of the *plane*, and the plane is what
Q1 named.

### What this does to earlier claims

* **`docs/FOOT_JOINTS.md` Q1's "the supine contact proxy's plane drop is the
  leading candidate mechanism" is WITHDRAWN.** The plane drop is worth ≤0.03 rad;
  the collapse is 1.25 rad. Q1's *factorial* stands unchanged — the collapse is
  the trunk repartition — and so does every number in it. Only the mechanism
  attributed to it is withdrawn, for the second time on this line: first the
  subtalar, now the plane.
* **The talus-impact table in Q1 is not withdrawn but no longer supports what it
  was cited for.** Its measurements reproduce; the impact is now known to be a
  co-symptom, not the cause.
* **`docs/ARTICULATED_SPINE.md` and `docs/WORKBENCH_AUTHENTICITY.md` §0.2 both
  carry the withdrawn sentence** and need the same correction. ARTICULATED_SPINE
  is corrected in this commit; WORKBENCH_AUTHENTICITY is outside this work's
  territory and is flagged, not edited.
* **G-S and F2 stay FAILED.** Neither is rescored, here or anywhere.
* **This is still the scaffold.** Nothing here is a statement about a human foot,
  ankle or trunk, and a pinned plane is an instrument, not a better bed.

---

## Q4 — pre-registration: is it the proxy RADIUS or the repartitioned INERTIA?

Written and committed **before any Q4 arm runs**. Q3 excluded the support plane and
left two candidates standing, both of which travel with the torso repartition:

* the **proxy sphere radius**, 0.2577 m on the base plant and 0.3090 m on the
  repartitioned one, unchanged by pinning the floor — a pure contact artefact;
* the **repartitioned inertia itself**, which changes the trunk's dynamics whatever
  the contact geometry.

The engine option added for this (`proxy_radius_m={'torso': <m>}`, commit to follow)
pins a named body's proxy radius instead of inscribing it in that body's inertia
ellipsoid. It is refused in the upright environment and refuses a body the model does
not have.

**Arms.** `supine_tonic`, TONIC 0.02, 200 steps of 10 ms, bar **0.2202 rad** — F2's
protocol line for line, as Q1 and Q3 used.

| arm | plant | torso proxy radius |
|---|---|---|
| `base` | base | derived (0.2577 m) |
| `tweld` | repartitioned, all new joints welded | derived (0.3090 m) |
| **`tweld@R0`** | repartitioned | **pinned 0.2577 m** — the base plant's ball |
| **`base@RT`** | base | **pinned 0.3090 m** — the repartitioned ball |

**What each outcome means**, fixed now:

| `tweld@R0` | `base@RT` | reading |
|---|---|---|
| recovers | collapses | the **radius** is the mechanism: a contact artefact, and the repartitioned inertia is exonerated |
| still collapses | stays fine | the **inertia** is the mechanism: not a contact artefact at all, and the sphere is exonerated |
| recovers | stays fine | the two are **not separable by this arm**; each is necessary and neither sufficient |
| still collapses | collapses | **both** are independently sufficient — report both, claim neither alone |

A fifth outcome is possible and will be reported as such: **neither** explains it, in
which case the factor is the repartition through some third route and Q1's factor
stands with no mechanism attached. That is the outcome this pre-registration exists to
keep available, because the mechanism has already been guessed wrong twice here.

**Instrument checks that must pass before any arm is read.** Pinning a body's radius at
its OWN derived value must reproduce the unpinned run bit-identically; `base` and
`tweld` must reprint Q1/Q3's 0.1229/0.1207 and 1.4014/1.2876; the pinned radius must be
read back out of each run's own `execution.json`, never inferred from the flag passed;
and the support plane must move exactly as the pinned radius dictates, since the plane
still hangs under the lowest ball.
