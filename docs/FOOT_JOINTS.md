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

## Results

Not yet run.
