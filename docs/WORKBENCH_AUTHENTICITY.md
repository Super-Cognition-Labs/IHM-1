# what in the workbench is not the real thing, and what it takes to make it real

A complete register of every shortcut, proxy, fallback and unreachable capability
between the live body and what the workbench shows. Read `ACTUATION_STAGES.md`
first: it says what the body is *for*, and this file says where the body is not
that yet.

Measured by reading the code on 18 September 2026, not surveyed. Every row names
the file and the line that carries the shortcut.

**The headline.** Nothing in the workbench fabricates data — no placeholder
geometry, no synthetic trajectories, no mocked responses; the app fails loudly
instead (`docs/APP.md`, and the client's only "standing in" string is a refusal).
The gap is not dishonesty. It is that **the body being integrated is a
22-segment scaffold with legs, and the body being displayed is 4,000 anatomical
entities carried along kinematically.** Everything below is a consequence of
that, or of a capability that exists in the plant and cannot be reached from the
app.

---

## Tier 0 — the plant is not a human body yet

These are not app problems. They are the model.

### 0.1 ~~There are no muscles above the pelvis~~ — WITHDRAWN 18 Sep 2026

**This entry was wrong, and the way it was wrong is the one this programme keeps
recording: a quantity measured correctly, on the wrong object.** It counted
`Millard2012EquilibriumMuscle` elements in the source file
`.../example3DWalking/subject_walk_scaled.osim` — 80, all attaching at the pelvis
or below — and reported that as the body. **Nothing on the live path loads that
file unaugmented.** The registration the runtime passes by default inserts more
muscles at load, so the file and the plant are different objects, and a static
parse cannot see it.

Measured by constructing each plant and counting `snapshot()['muscles']`:

| plant | muscles | above the pelvis | who loads it |
|---|---:|---:|---|
| bare source file | 80 | 0 | nothing on the live path |
| `whole_body_arm26_v2` | **92** | **12** | the `EmbodiedRuntime` default (`embodied.py:365`) |
| `engineering_stance_v1` | **98** | **18** | forced for every stance controller |

The 18 are `arm26_{BIClong,BICshort,BRA,TRIlat,TRIlong,TRImed}` and
`gait2392_{ercspn,extobl,intobl}`, both sides. The elbow, three shoulder axes and
three lumbar axes are bidirectionally muscled, and the upper body moves under
muscle alone with the torque ports at zero — right elbow +105.98° against +0.10°
contralateral in a paired null run.

The 13 `CoordinateActuator` torque motors are also not what that entry implied:
`coordinate_actuation` is passed by four offline scripts and by nothing in
`ihm/app/` or `ihm/assembly/`, so on the live path all 13 sit at zero and the
upper body is moved by muscle or by nothing.

**The girdle now exists, from a licence-clean donor** (`294f3a9`…`3a354ce`,
`docs/SHOULDER_GIRDLE.md`). `data/models/shoulder_girdle_v1`: **29 bodies, 60
coordinates, 158 muscles**, mass conserved to 1e-9, built from Seth 2019's
thoracoscapular model — CC BY 4.0 on SimTK *and* Apache-2.0 through opensim-core's
CLA, which **converge on commercial use with attribution**, the exact question on
which MoBL-ARMS's statements diverge. Clavicle and scapula per side, a
sternoclavicular and a 4-DOF scapulothoracic joint, the acromioclavicular constraint,
`acromial` re-parented onto the scapula, 30 donor muscles per side. The arm does not
move: <1e-9 m against the base plant. The torso debit reconstructs to **0.0 kg** and
7.1e-16 kg·m². 16 of 16 checks pass, re-run here.

**The donor's own inertia was not physically realisable, and was replaced rather than
clipped.** Verified here independently: its clavicle and scapula violate the triangle
inequality and its **radius carries a negative principal moment** (−4.5e-06 kg·m²).
They were rebuilt from this body's own anatomy.

**Can it press itself up? The girdle is not the limit — the elbow is.** Against a body
weight of 761.1 N, arms straight gives **3,265 N (4.29 body weights)**, bounded by the
scapula; at the worst-case lever it is **345 N (0.45 body weights)**, bounded by the
**elbow's 43.7 N·m**, which this build does not touch. Shoulder torques rose 4–31×
and are bidirectionally balanced for the first time. Before the girdle the question
had **no answer at all**: the humerus was jointed straight to `torso`, so a press-up
load reached the trunk through joint reactions with no muscle in the path. **This is
capacity at one pose, not demonstrated behaviour** — nothing integrates, and a
per-axis sum of `Fmax·|r|` is not a tension-feasible torque cone.

**Still missing:** forearm rotation (the Arm26 donor fuses it — `pro_sup` has a
2.16 mm moment arm), the wrist, the hand. Recorded from the build: on the left,
`scapula_elevation_l` and `scapula_upward_rot_l` run **opposite** to their right-side
namesakes, because the mobilizer's convention is not derivable from its declaration —
and assembly success was no test of it (16 of 48 candidate conventions assembled,
wrong ones included).

Two cautions measured alongside: the lumbar drive ran **30.1° outside its declared
±90° range** with nothing clamping it (see Tier 1 — joint stops are now
reachable), and the shoulder coordinates declare `±10 rad`, which is no range at
all to clamp to.

`docs/UPPER_BODY_ACTUATION.md` carries the full census, the moment-arm controls
and the acquisition position: **MoBL-ARMS 4.1 is already on disk** (50 muscles,
thorax/clavicle/scapula) and the source record marking it `unacquired` is stale.
The blockers are girdle bodies, a one-time torso mass partition, a coordinate map,
and an unresolved licence — not acquisition.

### 0.2 The skeleton is missing most of its joints

    WeldJoint  radius_hand_{l,r}     no wrist
    WeldJoint  subtalar_{l,r}        no inversion/eversion
    CustomJoint back                 ONE joint, 3 coordinates, for the whole spine

There is **no head body at all** — no skull, no neck, no cervical or thoracic
spine, no scapula, no clavicle, no fingers. `torso` is one rigid body from the
sacrum to the crown. The body cannot nod, turn to look, shrug, arch, grasp, or
roll its foot.

Consequence in the display: the skull, every vertebra, every rib, both scapulae
and both clavicles ride `torso` as one rigid lump, because there is nothing else
for them to ride.

*To close:* a cervical/thoracic chain, a shoulder girdle and a wrist in the
scaffold — or the scaffold's retirement per stage 3. Note the ordering trap: each
new joint must arrive with a passive stop (0.3 below), or it becomes one more
coordinate a search can exploit.

**Partly done, in a variant, 18 September 2026 — and the result is that joints
are the easy half.** `data/models/articulated_spine_v1` (`docs/ARTICULATED_SPINE.md`)
adds `head`, `cervical` and `thorax` bodies, a thoracic and two neck joints, both
wrists and both subtalars: 22 bodies → 25, 33 coordinates → 48, four `WeldJoint`s
→ zero, total mass conserved to 1e-9 because the torso was **repartitioned**
using the inertial records this repo had already measured and left unused
(`data/research/cervical_inertia/v2`, `data/research/thoracic_mechanism/native_composition_v1`).
Every new coordinate carries a declared range from a donor's own declaration and
a `CoordinateLimitForce` at that range, inside the model file. It loads and
integrates at +17% wall clock.

Two measurements made it a partial close. The first has since been taken
apart, and the explanation it gave is **withdrawn**:

* ~~No muscle has a moment arm about any of the fifteen new coordinates, because
  the source paths are fitted polynomials, so the missing subtalar arm is what
  collapsed the ankle.~~ The zeros were **two different problems**, separated by
  measuring each coordinate two independent ways — OpenSim on the model's own
  muscle geometry, and which bodies each path touches read from the file — which
  agree on all 21 coordinates checked:
  * **representation:** 11 muscles per side genuinely cross the subtalar (arms up
    to 32.7 mm) and read 0 only because the fitted polynomials take no subtalar
    input;
  * **absence:** nothing crosses the thoracic, neck or wrist joints at all. The real
    geometry reads 0 too, and no refit can help — they need muscles.

  **What is drivable now** (`1e29ee3`, `registration_muscled.json`: 25 bodies,
  **148 muscles**, 48 coordinates). The 22 subtalar-crossing foot muscles run on
  the model's own geometry, as the 18 arm and trunk muscles already did, with the
  other 58 paths bit-for-bit unchanged. Two refits were tried first, each
  pre-registered, and each FAILED its do-no-harm gate on the right knee
  (`gasmed_r` 1.66 mm, then `gaslat_r` 4.50 mm, against 1.36 mm). The fitter was
  also observed discarding 706 and 1,170 of 1,485 samples as NaN on two of four
  fits, for a reason not found, so that route was closed rather than redrawn. The
  neck gets **50 muscles transferred verbatim** from the MASI cervical donor —
  sternocleidomastoid, splenius, semispinalis, longus colli and capitis, scalenes —
  every one's rest length matching the donor's own to 1.5e-11 m, 160 of 168
  muscle–joint pairs agreeing in sign with it. Excluded: 10 anchored to a clavicle
  or scapula this body lacks, and 18 lying wholly inside the lumped cervical body.
  At full activation the neck extensors reach about 53 N·m against ~9.5 N·m to hold
  the head with the neck horizontal — **capacity, not demonstrated behaviour** — but
  flexion at the skull joint reaches only 1.12 N·m against the 2.69 N·m needed to
  lift the head from supine there.

* **Pre-registered gate G-S — does giving the subtalar its arms repair the
  ankle? FAILED** (1.140 / 0.871 rad against 0.220) — and its premise is now
  **withdrawn**. A full factorial over the variant's three changes — trunk
  repartition T, wrists W, subtalar S — pre-registered in `e00844d` before any arm ran,
  results in `a910288`, `docs/FOOT_JOINTS.md`:

  | arm | ankle r / l (rad) |
  |---|---|
  | base with only the subtalar freed | 0.113 / 0.108 — fine |
  | base with only the trunk repartitioned | 1.385 / 1.266 — **collapses** |
  | variant with only the subtalar re-welded | 1.385 / 1.267 — does not recover |
  | repartitioned trunk, **every new joint welded** | 1.401 / 1.288 — **collapses** |

  Effect on the worse ankle: T +1.254 rad, W −0.001, S −0.019. **The collapse is the
  trunk repartition, and it is not the spine's new articulation** — it persists with
  all seven new joints welded. "Plantarflexors loading an uncontrollable hinge" is
  excluded too: at the crossing they carry 25–31 N against 468–476 N in the
  dorsiflexors, and the ankle crosses before the subtalar reaches its bound.

  ~~**Leading candidate: the contact proxy's PLANE.**~~ **WITHDRAWN 18 Sep 2026** —
  the second mechanism withdrawn on this question, after the subtalar. The engine
  option added in `3a6e8f0` made the separating run possible, and it was
  pre-registered (`e8d25d8`) before running (`293b913`, `docs/FOOT_JOINTS.md` Q3):

  | arm | floor | ankle r / l |
  |---|---|---|
  | repartitioned trunk, floor **pinned UP** at the base plant's −0.403499 m | held | **1.3961 / 1.2820 — still collapses** |
  | base trunk, floor **dropped** to the repartitioned −0.452013 m | held | **0.1532 / 0.1504 — still fine** |

  Moving the floor 48.5 mm is worth **−0.0054 rad** on the repartitioned plant and
  **+0.0303** on the base one; the repartition is worth **+1.2732 rad with the floor
  pinned**, against +1.254 with it free — factor and floor do not interact. A 5-rung
  ladder is monotone and graded with no threshold (0.1229 → 0.1532, slope 6.26e-4
  rad/mm); at that slope the floor would need ~2.0 m to reach 1.40.

  **The impact story is falsified by its own witness.** `base@PT` lands its heels at
  992/925 N — 2.7× the base arm's, as hard as the collapsed arms — and its ankle moves
  0.03 rad; `tweld@P0` lands *softer* than `tweld` and folds anyway. Heel impact is a
  co-symptom.

  **What survives.** Pinning the plane does not pin the BALL: the repartitioned torso's
  proxy sphere is **0.3090 m against 0.2577 m** at every arm, wherever the floor is, so
  that trunk's contact point sits 51.3 mm further from its own mass centre. Two
  candidates remain and this run does not choose between them — the **proxy radius**
  (still a contact artefact) or the **repartitioned inertia itself**. The separating arm
  needs an engine option for the per-body proxy radius: run `tweld` with the torso ball
  pinned at 0.2577 m.

* **F2's bar compares plants lying on planes 48.5 mm apart.** Recorded; F2 stays
  FAILED and is not rescored.
* The original gate — no new coordinate may leave its range by more than the base
  model already leaves its own — still **FAILED** on `thoracic_extension` (0.394
  rad against 0.220). Recorded, not rescored.

**Still blocked:** the thoracic joint has no donor muscles (the only candidates'
attachments sit 32 mm above a joint centre that is itself uncertain to 62 mm); the
wrists are blocked on MoBL-ARMS 4.1 and its unresolved licence, though NOT on the
girdle, since the wrist muscles need no clavicle or scapula; the 10 girdle-anchored
neck muscles wait on the girdle.

Still absent in that variant: scapula, clavicle, fingers, per-vertebra thoracic
bodies, and any actuator on the new coordinates. The shoulder girdle is absent
for a stated reason — the thoracic material partition is rib cage only, so there
is no mass in it to give a scapula, and both `acromial` joints still sit on
`torso`. That is the neck-first half of the double torso partition
`docs/UPPER_BODY_ACTUATION.md` warns about, and it composes: the debit is done
with `ihm.assembly.cervical_inertia.partition_body`, exact to 1e-15, so the
girdle's 0.86 kg per side can come out of the 17.938 kg residual core afterwards.

In the display nothing changes yet: `data/derived/anatomy-segment-binding/`
assigns every entity to one of the original 22 segments, so the skull still rides
`torso` until that binding knows the three new bodies exist.

### 0.2a The identified plant has two toe hinges no muscle can drive

Found while diagnosing the spine variant, and it is **not** a property of the
variant — it is in `engineering_stance_v1`, the plant every existing result was
measured on. Upstream fitted the muscle-path polynomials **with the toes welded**
(`exampleMocoInverse.cpp:51`: `ModOpReplaceJointsWithWelds({"mtp_r", "mtp_l"})`);
this repo's engine keeps `mtp` as a free `PinJoint`. So the toe muscles cross a
joint their paths cannot see. Measured on the base plant through the engine:

| muscle | about `ankle_angle_r` | about `mtp_angle_r` |
|---|---:|---:|
| `edl_r` | +39.3 mm | **0** |
| `ehl_r` | +42.6 mm | **0** |
| `fdl_r` | −11.3 mm | **0** |
| `fhl_r` | −18.0 mm | **0** |
| `soleus_r` (control) | −49.7 mm | 0 |

These are the toe extensors and flexors. Both toe joints are therefore unactuated
hinges held only by passive terms.

**What it costs, measured** (`docs/FOOT_JOINTS.md`, base plant only, untouched):
**almost nothing in the protocols this repo runs.** The toe joint stays within
±0.009 rad supine and never comes closer than 0.38 rad to its ±0.524 bound in the
crawl, held by its passive spring. Welding it, as upstream did, moves the ankle by at
most 0.0016 rad and crawl travel by under 1 mm; F2's bar moves 2e-5 rad and no gate
flips. The same runs reproduced three numbers recorded weeks earlier — the unstopped
crawl's 1.450 rad, the stopped crawl's 0.220 rad, F2's 0.2202 — which also confirms
on the rebuilt engine that the §1.3 damping fix moved no plant.

**The defect underneath it is now diagnosed exactly, and fixed in a variant**
(`fb9b90d` pre-registered, `118096e` results, `docs/FOOT_GEOMETRY.md`).
`subject_walk_scaled.osim` is `Rajagopal2016.osim` through OpenSim's ScaleTool, which
multiplies every joint offset by its own parent body's factors. Audited over the
shared offset frames — and verified here independently against both source files —
**exactly two nonzero offsets are byte-identical to the unscaled original, and they
are the left and right toe joints.** The ankle's went −0.400 → −0.465 m, the
subtalar's −0.0488 → −0.0599; `mtp`'s stayed at 0.1788. The same two frames are also
the only ones whose `orientation` was dropped, from Rajagopal's oblique
`−3.14159 ±0.619901 0` to a plain `0 0 0` hinge. Rajagopal ships `mtp_angle` **locked**,
which is why it never showed upstream.

The correction is the rule the other offsets obey — Rajagopal's own translation times
this model's own calcn scale factors, with its orientation restored — giving
**0.1788 → 0.20663 m (+27.82 mm)**. Nothing is typed; the builder re-derives it and
refuses if the audit stops holding.

**The sign pattern is now anatomically possible**, and the number was predicted before
the engine saw the model and reproduced to **0.0004 mm**: extensors `edl` −6.05 and
`ehl` −7.55 mm, flexors `fdl` +6.83 and `fhl` +7.13 mm. Previously `ehl` read **+24.1**,
sharing a sign with the flexors. It was the **offset**, not the axis, that carried the
sign error. Toe paths on the corrected geometry now *reduce* ankle excursion
(−0.013 supine unstopped, −0.076 in the crawl) where the shipped geometry *raised* it
by +0.078.

**It does not repair the unstopped crawl** (1.37 rad), F2's bar moves ≤0.0011 rad, and
neither F2 nor G-S flips. `mtp` is also welded in a separate spine-variant plant
(`model_mtp_welded.osim`); the largest change to any of its 46 coordinates is
0.0061 rad. The recommendation to leave `engineering_stance_v1` alone is **reinforced**:
the corrected foot is 27.8 mm longer, so adopting it means re-identifying that plant.

**Recommendation from that evidence:** leave `engineering_stance_v1` as it is (its
stance linearization was solved with mtp free at 0.113 rad, so changing it means
re-identifying the plant); weld mtp in any new plant; fix the foot geometry — scale
the mtp offset, ideally adopt Rajagopal's oblique mtp axis — before giving the toes
geometry paths. **Not measured:** upright stance and push-off, where toe loading
matters most.

### 0.3 Two bodies of different stature and different mass

| | mechanical | anatomical |
|---|---|---|
| stature | 1.7973 m | 1.7195 m |
| mass | 77.6122029 kg | 70.7713 kg |

**Audited 18 September 2026** (`docs/BODY_PARAMETERS.md`, "The two-body constants,
counted"). Two things changed and one did not.

**The mechanical mass has a derivation, and it is not a mechanical one.** It was
recorded here, in `CLAUDE.md` and in `ihm/body_parameters.py` as an unattributed
literal. It is the **BioGears StandardMale physiology patient's weight at t = 0**:

    170 lb x 0.45359237          = 77.1107029 kg   the declared patient
    + 0.5 L water + 500 mg Ca + 1 g Na =  0.5015 kg   StomachContents at t=0
                                   -----------------
                                   = 77.6122029 kg

exact in IEEE-754 double arithmetic, and present verbatim as
`<Weight unit="kg" value="77.6122029"/>` in 8 of the 19 BioGears `*@0s` state
files this repository already ships
(`data/raw/physiology/biogears/share/data/states/StandardMale@0s.xml:8`).
**So the plant is scaled to the fed weight of a physiology reference patient** —
the same 170 lb this repository's own composition ledger rejected as unreachable
in the acquired envelope. That is a modelling question, not a tidy-up.

**The count was wrong, and the instrument is why.** "67 places across 58 files"
was a stale 2026-09-10 figure — but the repository cannot have *one* such number,
because the shell `grep` here is a wrapper passing `--ignore-files` and silently
skips everything in `.gitignore`. `grep -rn` from the prompt returns **124** lines
in 0.031 s; `/usr/bin/grep -rnI` over the same tree returns **67,932** in 42,554
files, and `git grep` (tracked only) returns **109 in 82**. Of the 67,932,
**67,871 are under `data/`** — recorded history, and not to be changed.
**Zero** lines under `ihm/` now carry it outside `ihm/body_constants.py`; 22 live
literals in 20 `scripts/` files still do, four of them added by other work the
same day.

**The seam is now declared in one place.** `ihm/body_constants.py` holds each
constant once with its provenance; `ihm/body_parameters.py` imports and
re-exports them; `scripts/verify_body_constants.py` fails if a raw literal
reappears under `ihm/`, re-derives the mass two independent ways, and checks the
two declared mirrors (`garment_wardrobe.py`, `anatomy_pose.py`) for exact
equality. It has a `--self-test` that makes each check fail on purpose.

**What did not change:** nerve routes still scale by the anatomical stature and
the displayed anatomy still scales by the mechanical one
(`scripts/materialize_body_variant.py:193`, declared `known_seam`). The seam is
declared, not closed.

*To close:* one body, one stature, one mass. The literal was the easier half and
is done; the scale seam is the real work, and it now has exactly one place to be
changed from.

---

## Tier 1 — capabilities the plant has that the workbench cannot ask for

**This is the cheapest and most embarrassing tier.** `NativeMechanicalStream`
accepts these; `ArticulatedBodyPlant.__init__` (`ihm/assembly/articulated.py:119`)
simply does not forward them, and `EmbodiedRuntime.from_workspace` never offers
them. Every one is reachable only from an offline `scripts/*.py`.

| capability | what it is | who can use it today |
|---|---|---|
| `coordinate_limits` | the opt-in joint stops | `scripts/crawl.py`, `measure_tissue_mechanics.py`, `scene_object_contact.py` |
| `segment_contact_meshes` | **real segment surfaces as contact geometry, including the `skin` layer** | `scripts/measure_segment_contact_meshes.py` |
| `segment_contact_material` | its E, ν, layer thickness | same |
| `segment_contact_replaces_source_feet` | | same |
| `tissue_ligaments` + 3 companions | the **117 ligament/capsule force elements** | `scripts/crawl.py`, `measure_tissue_mechanics.py` |
| `scene_objects` / `scene_contact_material` | native sphere colliders | `scripts/scene_object_contact.py` |

So the live body in the workbench right now:

* **has no joint limits** — and the model declares `<clamped>true</clamped>` on
  every rotational coordinate while holding **zero** `CoordinateLimitForce`
  (`docs/NATIVE_JOINT_LIMITS.md`). A prone search once bought 973 mm of travel
  with the ankles folded to 145°, and that result was withdrawn.
* **stands on 28 spheres inscribed in inertia ellipsoids**, one per segment at
  its centre of mass (`docs/SEGMENT_CONTACT_SURFACES.md`: "a femur represented
  as a ball"), when the code to stand on real concave segment surfaces — and on
  the *skin* layer specifically — exists and is measured.
* **carries none of the 117 tissue force elements** that `docs/TISSUE_MECHANICS.md`
  built and evaluated.

### 1.1 CLOSED, 18 Sep 2026 — the keywords are forwarded

`ihm/assembly/plant_options.py` resolves server-owned identities into plant
keywords, in the same shape as `controller_selection.py`: a client names
`{"joint_stops": true, "segment_contact": "skin", "tissue_ligaments": "admissible"}`
and never a path, a mesh or a material. `ArticulatedBodyPlant` forwards them,
`EmbodiedRuntime.from_workspace` accepts them, and `EmbodiedSessions.create`
validates them **before a resource slot is taken** — an unknown bundle or a skin
selection in a supine environment is a 400 on the create call, not a body that
dies in initialization and leaves the caller reading an engine log.

**`None` is the historical plant to the float.** Nothing in the resolver runs
unless it is named, so every existing measurement is unaffected and the
workbench's controls send absence rather than an explicit off.

What a live body can now be asked for, measured by constructing each plant:

| selection | what the plant does |
|---|---|
| nothing (default) | 28 `fall_support_*` / `contact*` spheres — the historical body |
| `segment_contact: skin` | 20 `mesh_support_skin_*` elements, **source feet replaced** |
| `segment_contact: skin_carried` | the same 20, **source feet kept** |
| `segment_contact: skin_layer_map` | the per-patch measured-depth variant |
| `segment_contact: bone_all` / `bone_proxy` | real bone surfaces, declared as a collider layer |
| `joint_stops: true` | 22 coordinate limits, 8 coordinates left free |
| `tissue_ligaments: admissible` | **66 of 117 elements**, the kinematically admissible subset |

**What this does NOT establish, and an earlier draft of this entry wrongly claimed:
that the body stands on its skin.** Installing 20 skin contact elements is not the
same as those elements carrying the weight.
`scripts/measure_segment_contact_meshes.py` records — and **excludes the arm from
its own default run** for it — that *the skin never reaches the floor in the stance
pose*. The foot's skin surface sits above the source foot spheres' effective plane,
so replacing those spheres leaves a standing body with nothing under it, and it
falls. That is a property of the pose, not of the skin. Both arms are therefore
offered, neither is called the truth, and the selection carries the caveat: `skin`
replaces the feet (the intended end state, which collapses from this pose) and
`skin_carried` keeps them (the cost of carrying the geometry, separated from the
cost of a plant that is collapsing — two things one number would mix).

The bundle's **own declared skin material** now travels with it — E = 3000 Pa,
ν = 0.45, thickness 6.6 mm, from this body's canonical skin-layer entities — so the
arm measures what this body's declared skin does rather than what a default
contact stiffness does. The first version of the resolver dropped it.

Two things the resolver refuses to hide. The **8 unranged coordinates** — the six
shoulder angles and the two knee couplers — declare ±10 rad, which is the absence
of a range, so they are left free and listed by name on every resolution rather
than stopped at a limit the model never gave. And the **unfiltered** tissue set is
offered only under the id `all`, carrying the sentence that
`docs/TISSUE_MECHANICS.md` measured it **worse than no tissue on every drop**.

The joint-stop table is not new arithmetic: `ihm/assembly/plant_options.py`
reproduces `scripts/crawl.py`'s 22 rows **identically**, checked as a known answer,
so a stopped live body and a stopped offline crawl are the same plant.

**One unit defect inside that stop, found 18 Sep 2026 and not fixed here.**
`scripts/native_mechanical_stream.cpp:293` converts the stop's limits, its
stiffness and its transition into `CoordinateLimitForce`'s degree convention and
**does not convert `damping`**, which that class reads in N·m/(degree/s). The
declared `damping_nm_s_per_rad: 1.5` is therefore applied as 85.9 N·m/(rad/s),
57.3× the stated value, while `execution.json` records the stated one. It sits
inside a term that is already an explicit engineering constant, so nothing
published rests on it — but the number in the receipt is not the number the plant
used, and a sweep over damping would have been reading the wrong axis. The engine
is a shared binary and changing it changes the plant for every prior run.
Detail: `docs/ARTICULATED_SPINE.md`.

**The frame stopped lying about its own contact.** `body_environment.scope` was a
constant asserting "inertia-derived spheres are engineering proxies, not anatomical
surface". With a skin bundle loaded that sentence is false, so it is now derived
from what is actually loaded, and the per-frame `limitations` list gains an entry
when stops or tissue are absent rather than only when they are present.

**Reachable from the workbench**: a Contact select, Joint stops and Tissue forces
checkboxes, and an Anatomy pose select, in `app/src/scene-interaction.js`. Their
meaning is a value a test can read (`app/src/mechanical-fidelity.js`,
`app/test/mechanical-fidelity.test.js`), following `transport.js`'s pattern. The
note under them states what is **not** true of the body rather than congratulating
the reader for switching something on.

Verification: `scripts/verify_plant_fidelity.py`.

### 1.2 Still open

`scene_objects` / `scene_contact_material` are still not forwarded — native sphere
colliders remain reachable only from `scripts/scene_object_contact.py`. The
workbench's draggable objects are display geometry to the embodied path.

### 1.3 ~~A unit defect in the shared engine~~ — FIXED 18 Sep 2026, and no plant moved

`scripts/native_mechanical_stream.cpp` converted a joint stop's limits, stiffness
and transition from radians to degrees for `CoordinateLimitForce`, and passed
`damping` straight through. OpenSim's header declares that property
**`Nm/(degree/s)`** for a rotational coordinate, so a declared 1.5 was applied as
**1.5 × 180/π = 85.94 N·m·s/rad — 57.3× the number**, on every stopped run.

**The fix, and why it changes no result.** The engine now divides damping by 180/π
like the stiffness. The constants in `scripts/crawl.py` and
`ihm/assembly/plant_options.py` are re-declared at **85.94 N·m·s/rad — the value the
plant has actually been running** — so the engine receives exactly what it always
did. `(1.5·180/π)/(180/π) == 1.5` holds exactly in IEEE-754, so this is not
"approximately the same plant", it is the same bits. Measured:

| check | result |
|---|---|
| same pre-fix build run twice (is identity even testable?) | max \|Δq\| **0.0** — deterministic |
| pre-fix build, 1.5 declared vs post-fix build, 85.94 declared, 50 steps × 33 coordinates | max \|Δq\| **0.0 — bit-identical** |
| **control that must fail:** post-fix build fed the OLD raw 1.5 (57.3× less damping) | max \|Δq\| **0.444 rad** (`hip_flexion_r`) — the damping engages, so the identity is informative |
| `plant_options.joint_stops()` vs `crawl.joint_stops()` | identical |
| intake-mass variant engine, rebuilt | loads and steps |

So no stopped run has to be redone: every recorded number was measured at 85.94 and
the plant still delivers 85.94. The build is `data/runtime/mechanical-stream/build-cx6s8y89`
(variant `build-o25sa609`); the pre-fix build is stale by design and reproducible
from git history.

**What the corrected number says.** Against a critical damping of roughly
2√(kI) ≈ 6 N·m·s/rad for a limb segment at k = 30 N·m/rad, 85.94 is about 14×
overdamped. That is likely part of why a stopped plant integrates *faster* than an
unstopped one (`crawl.py`'s sweep). It is an engineering constant, declared as
such, and it has not been retuned — retuning would move every result.

**A second provenance defect, in this register's own resolver, fixed in the same
change.** `plant_options.py`'s `firm` and `stiff` profiles carried damping 2.0/3.0
and transitions 0.25/0.20 that no measurement supported — they had been typed. The
only sweep in the repo varied stiffness alone at the damping and transition above,
so both profiles now use those, and a resolution says which profile was adopted from
the sweep and which was only swept.

The spine variant's 15 in-model stops (`data/models/articulated_spine_v1`) write
damping 1.5 directly into the `.osim`, in OpenSim's own Nm/(degree/s) — also
85.94 N·m·s/rad. Unaffected by the engine change, and consistent with the base plant.

---

## Tier 2 — contact, and the fact that the skin never touches anything

`ACTUATION_STAGES.md` is explicit: *"Contact with the world is never bone against
world… it is the skin that meets the floor."* Today, in the upright environment
the workbench actually runs, contact is:

* 12 source foot contact spheres, plus
* one `SmoothSphereHalfSpaceForce` per non-foot body, radius inscribed in its
  inertia ellipsoid at its centre of mass.

The floor does carry the whole body (761.5 N against 761.4 N of weight), so it is
not wrong — it is a body-shaped collection of balls.

**Supine is better and is the exception.** A bed selection from the environment
catalogue sets `surface_contact_manifest`, which is a genuine
`ihm.supine-skin-foundation.v1` skin foundation with a measured bed material
(`ihm/native/mechanical_stream.py:98`). So the workbench's most authentic contact
is lying down, and its least authentic is standing up — the opposite of what the
programme needs.

Further limits on the skin path that is there:

* ~~at most **128** skin sensor points per session~~ — **FIXED 18 Sep 2026.** The cap
  was a magic number with no stated basis, and it made a whole body's cutaneous
  afference 128 numbers. Both sides now bound the selection by the quadrature's OWN
  size: **21,382 points**. Measured: 500 sensors bind and all 500 are emitted in a
  frame (`surface_foundation.sensor_points`), where 129 used to be refused; asking for
  more points than exist is refused by name. A plant that selects no sensors is
  bit-identical across the change (`scripts/verify_engine_options.py`).
* no self-contact, no skin-on-skin, no sliding, no volume preservation.

**And the supine skin foundation could not be loaded at all.** The environment
catalogue's `bed-support-skin-quadrature` selection pins
`data/derived/supine-surface-contact-exmzq9pq/manifest.json`, whose `source_files`
hashes no longer match: `canonical/mechanics.json` (last written 8 Sep) and
`scripts/build_supine_surface_contact.py` (changed 9 Sep) both moved, so
`NativeMechanicalStream` refuses it with `Surface foundation source identity
mismatch`. **This predates all of today's work by ten days** — the single most
authentic contact in the workbench, the bed foundation, has been unselectable. The
artefact regenerates in **3.1 s** (`scripts/build_supine_surface_contact.py`).

**Fixed the same day.** Of the five supine artefacts on disk, exactly two validate
against current sources, and the regenerated one is byte-identical to the existing
valid one in the data that matters (`quadrature.npz` and the foundation file; only
manifest metadata differs), so no new artefact was kept. The catalogue builder's
three pinned paths now point at the valid `supine-surface-contact-5jqy1juo` — the one
the nociception work already binds to — and the catalogue was rebuilt (90 artefacts,
its own checks passing). Verified end to end: `resolve_selection(... 'bed-support-skin-quadrature')`
returns the valid manifest, and a plant built from it **loads and steps**, on the
foundation's own support plane (−0.269 m) with 20 contact elements, where it
previously refused to start.

*To close:* Tier 1's wiring gives upright the skin layer. Beyond that: a
deformable skin with self-contact, and an afferent allocation that is not capped
at 128 points.

### 2.1 A deforming soft-tissue layer exists — and is NOW in the plant loop (18 Sep 2026)

`docs/SOFT_BODY.md`; built in `a1b2812`, made convergent, depth-correct and fast in
`17ecc4e`. Per segment, tissue from the skin down to the body's own measured depth
is meshed and solved with the same constitutive law as the shipped supine
foundation; its reactions on the rigid core are the force and moment on the bone.
Selectable as `soft_tissue: layer_fitted_local_confined`. **It is the first
implementation of any part of ACTUATION_STAGES' "participant" mode.**

**COUPLED, later the same day** (`docs/SOFT_BODY.md`, first section;
`scripts/verify_soft_tissue_coupled.py`, 32 gates, 0 FAILED over two runs). A
`*_coupled` identity makes `ArticulatedBodyPlant.advance` pose each coupled segment's
layer at that segment's own transform every step and apply what it transmits through
the same `forces=` port a caller uses. The sentence that used to stand here — *nothing
it computes flows back into the integrator* — is withdrawn. What the measurement says,
and none of it is comfortable:

* **Off is off, bit for bit.** An unloaded coupled plant is the historical plant in
  every coordinate for 81 steps, because a segment whose skin is clear of the support
  emits no force port at all. An unloaded segment costs below the timer's resolution.
* **Sub-cycling was built, measured, and NOT adopted.** Against a per-step reference,
  holding the reaction for 2 steps already costs 26% of the peak reaction and moves the
  plant 1.45x further than the layer's own 15% force uncertainty does. The reason is one
  number: the reaction grew at **1,473 N/s**, so a 10% bar allows a **1.0 ms** hold — a
  tenth of the plant step. A cadence is set by how fast the contact develops, not by
  what the solve costs. The coupling runs at one solve per step and is **5.4x short of
  real time** (53.9 ms of coupling per 10 ms step, on ONE 14,874-DOF segment).
* **Not at the foot.** At the scaffold's own stance the foot's skin sits **36.3 mm above
  the floor** while the source foot spheres carry 616 N, and it still clears by 32.4 mm
  at the bottom of a 770 N landing — so a coupled `calcn_l` returns **exactly 0 N** in
  every upright trajectory the scaffold can produce. And if it touched: 5.57 N at 2.5 mm,
  leaving the constitutive domain at 3.0 mm, against 761 N of weight.
* **It does not replace the engine's own contact.** No option here removes one segment's
  engine contact, so a coupled segment that has one carries both. The double count is
  reported in every frame beside the layer's own force rather than left to be inferred.
  The measurement was run on `ulna_l`, chosen because it is the segment with the longest
  window in which its skin is below the floor and its engine element is not — 69 steps,
  during which the engine's force on it peaked at 1.7e-04 N against the layer's 14.7 N.
* **When the rigid core would enter the support the step RAISES** and the plant is left
  exactly where it was. Gated, with the plant's coordinates compared bitwise across the
  raise.

**The first heel result is withdrawn.** Every earlier "heel" figure — the
1.113 / 0.457 / 0.657 N non-converging sequence, and the "4.4–17× softer than the
column it replaces" — measured `calcn_l`'s lowest skin vertex, which is on the rim of
the cut between the heel and toe segments, where tissue is at most 1.7 mm thick. That
was a seam wedge, not the heel pad. A real heel pose is now derived from the mesh (the
smallest rotation keeping the forefoot off the floor, 17.26°), not typed.

**Depth.** The "rigid core 30 mm above the sole" was measured between two points
55 mm apart. Under the heel the core sits at the declared depth. What was wrong is
the declared depth itself: one median per segment (18.61 mm), while the measured
depth map reads 6.4–14.1 mm under the heel and 3.9–8.9 mm under the forefoot, so the
median left the forefoot with no core at all. A local rule now reads the map point
by point — 6.9 mm heel, 5.9 mm forefoot, the latter inside the published fat pad
(5.8–7.1 mm). `radius_l` fails that rule on one vertex two skin partitions disagree
about and is refused rather than passed.

**Convergence — monotone, first order, not converged.** The mesh is now cut to the
skin and the core surface, with nodes exactly on both. Heel force by cell size is
monotone under both depth rules, but successive differences do not shrink from the
first step (CV8, CV10 **FAILED**, recorded). Observed order ≈ 1; Richardson
extrapolation puts the finest forces ~8% and ~15% high. **No contact force from this
layer is converged to better than about 10%.**

**Speed, and how far the plant loop is.** A cold heel solve went from 2.0–3.3 s to
0.17–0.44 s at 190 MB (fill-reducing ordering, numeric factorisation reused across
iterations and steps, an exact Hessian, a line search that no longer stalls). A
loaded, warm-started heel step costs **136–159 ms median against the plant's 10 ms**
— 9–16× short of real time (RT1h **FAILED**). Those absolute times were measured on a
heavily loaded machine; the ratios are fair, the absolutes are an upper bound. A
small-strain condensation reaches 1–3 ms but reads 2–27% low, so it was measured and
not offered.

**The amendments were audited before this was accepted.** The pre-registration
(`e140edf`, 18:40) was amended three times; every gated run — the first after it is
run 5 at 19:33 — postdates all three, and the census was pre-registered (19:39)
before its first run (20:05). A depth-map *build* at 18:59 preceded the 19:17
amendment; it is geometry, and the amendment concerns geometry.

*Still to close:* an IMPLICIT coupling (an explicit reaction is a step behind by
construction, and on a contact this fast even one step is 10x outside the bar),
real-time cost, a converged mesh, a plant option that removes ONE segment's engine
contact so the layer can replace it instead of adding to it, bone as the rigid core
(skin and bone are not co-registered in the bundle), muscle as a volumetric
activation-coupled solid, and sliding between soft surfaces for fascia and bursae.

---

## Tier 3 — the control path

### 3.1 The control path: a recruitment pool, and hip/knee reflexes from the source

**Two corrections to how this entry used to read.** The reflex source is **Geyer &
Herr 2010** (IEEE TNSRE 18(3):263), not "Song & Geyer": the pinned file's sha256
matches `SOURCE_SHA256` and its title page says so. The wrong name came from the
download URL's filename, `song.pdf`, and was carried from here into two work
orders. And the headcount was not "80 muscles" — the live plant has 92 by default
and 98 for stance (§0.1).

What the entry used to describe is still the default, deliberately: eight ankle
muscles `{tibant,soleus,gasmed,gaslat}_{l,r}` with Geyer & Herr's Table I constants,
and every other muscle driven by `baseline + gain·reflex + drive` with gain a
linear function of a regional firing rate (`cortical_gain_per_hz = .005`). That is
a dial, not recruitment, and `None`/`None` still selects it — reproduced field by
field against the original file at `ab3d2b8` over 30 of 30 frames
(`scripts/verify_recruitment.py`, all pass).

**What exists now beside it** (`docs/MOTOR_RECRUITMENT.md`):

- **A size-principle recruitment pool** (`ihm/assembly/recruitment.py`) — the
  Fuglevand, Winter & Patla 1993 motor-unit pool, every number from the
  open-access restatement by Potvin & Fuglevand 2017, pinned by sha256. 120 units,
  recruited small-first, firing from 8 imp/s to a 35→25 imp/s ceiling; the last unit
  recruits at 74.6% of maximum drive, as the paper states. It produces excitation
  only; native mechanics still owns activation, fibre dynamics and force. It is
  **one pool shape for every muscle** — no source read supports ordering units
  across muscles by Fmax or fibre length, so it does not.
- **Knee and hip reflexes** (VAS, HAM, GLU, HFL, stance and swing) transcribed from
  Geyer & Herr's Table I and Appendix I, not by analogy. Two defects in the source
  are recorded rather than fixed silently: the printed stance HFL law can only
  lower stimulation, contradicting the paper's own p. 265 (the p. 265 reading is
  used); and `k_bw = 1.2` lies outside its own printed range of 1.3–5.0. Declared as
  unmeasured choices: which plant muscles map to each source group, how one source
  muscle is split across several, and that all reflexes run at 20 ms where the
  source uses 20/10/5 ms for ankle/knee/hip.

**Measured against EMG, pre-registered before scoring** (commits `313f7db` then
`02b3861`), one right-leg walking trial on the model the plant descends from:

| gate | result |
|---|---|
| P1 — knee/hip reflexes beat a foot-contact indicator | **FAILED**: skill +0.023 vs +0.040, 95% CI of the difference −0.124 to +0.126 |
| P2 — pool vs rate decoder | **null, as predicted**: with no descending drive the pool is idle and the arms differ by a dial moving 1.005–1.103; CI −0.053 to +0.019 |
| P3 — context, not a gate | the existing ankle reflexes carry signal beyond gait phase: 0.54 vs 0.09, soleus r = 0.85 |

Medial hamstrings is predicted **backwards** — its +0.065 skill exists only because
the fitted map has a negative slope, and its raw correlation is −0.25. So skill
after a fitted map can reward a law that is wrong in sign; read the raw correlation
beside it. The comparison ran on a **replay**, not the native plant, because the
engine cannot prescribe motion and its static evaluator returns no tendon force.

**What still stops it:** the pool can only be tested on a task with a measured
descending command, and none is retained. Everything above is one subject and
about one stride.

### 3.2 The cortico-cortical path carries 0.03–0.07% and severing it changes nothing

`ihm/app/brain_graph.py:338`, and the app draws this honestly:

> `'Severing the kernel changed the motor command by nothing, because nothing was
> getting through to sever.'`

### 3.3 The one load-bearing cortical path does not depend on what the cortex learned

The stance policy owns all 98 muscle commands and is genuinely load-bearing —
intact holds at 0.17 mm, severed falls at 2.90 s. But a kernel with its **site
rows permuted** recovers the same push (1.7452 vs 1.7844 mm) and falls at the
same 2.90 s. The E/I dynamics are the motor owner; the *trained content* is not
doing the work.

### 3.4 There is no metabolic feedback onto force

`ihm/assembly/embodied.py:646` raises `'Native muscle energy demand is unmet;
mechanical supply feedback is not yet supported'`. A body that cannot afford its
own contraction aborts the run rather than fatiguing.

**Partly closed, opt-in, 18 Sep 2026** (`228b154`, `765e649`, `docs/METABOLIC_SUPPLY.md`).
The default is unchanged and still aborts with the original message.
`metabolic_supply='energy_fraction'` instead scales the next interval's excitation, for
every muscle including uncommanded ones, by the fraction of the movement's extra
energy demand that native fuel actually covered, and books what supply did not
cover in a separate "unsupplied" account rather than dropping it.

**It is declared as an accounting rule, not physiology, because no published
fatigue law runs on what this engine exposes.** Xia & Frey Law 2008 takes no
metabolic input; Callahan, Umberger & Kent-Braun 2016 needs phosphate, pH and PCr,
which BioGears does not model; Allen/Lamb/Westerblad 2008 and Ørtenblad et al. 2013
describe mechanisms, not a law. "Unmet" means the muscle's glucose and its single
whole-body glycogen store both ran out within one 20 ms interval, after amino acids,
fat, and aerobic and anaerobic glucose and glycogen were drawn in that order
(`Tissue.cpp` 1165–1484); every run on record reports zero unmet.

Verified here, not taken from the report: `scripts/verify_metabolic_supply.py` 16 of
16, including the default reproducing the pre-change `embodied.py` exactly, the energy
ledger closing to 2.8e-14 J over 84 exchanges, and four fail-on-purpose controls that
fire. **Two limits stand:** the cap acts one interval late (the engine cannot preview
a step), and it scales excitation, not energy, so the next interval is not
guaranteed to fit the supply. **Not measured:** behaviour when BioGears genuinely runs
short, which needs a glycogen-depleted starting state; an oscillation near the 2 s
exchange lag is expected and unmeasured.

*Remaining, in order:* a fatigue law with a mechanism the engine can supply
(phosphate and pH are not modelled) → a depleted-state run to see the cap act.

---

## Tier 4 — what the renderer does to the truth

### 4.1 The verified anatomy poser is not called

`ihm/assembly/anatomy_pose.py` exists, and passes an unusually good battery in
`scripts/verify_anatomy_pose.py`: FK matching Simbody to 7.8e-16 over 12 native
frames, exact rest pose to 5e-16, distal-only motion with 0 violations, bit
identical idempotence, frame guards that each fire, and 11 mirrored-pair
overrides.

**Nothing calls it.** `ArticulatedBodyPlant` builds `CanonicalRegistration`
instead (`articulated.py:127`), whose own manifest describes itself as *"one
shared global rigid fit; COM and bounding-box center are not measured homologous
joint landmarks"* and whose per-frame limitations include *"Canonical anatomical
joint locations and surface continuity remain uncalibrated"*.

So the workbench renders the true anatomy through a demonstrably coarser
registration than the one the repo has already verified.

**CLOSED, 18 Sep 2026 — and it was worse than "coarser".** Measured before
wiring anything, at the plant's own first frame with `pelvis_ty = 0.93`:

| | translation of every entity |
|---|---|
| shipped `registration.project()` | **1.28e-16 m** — zero to machine precision |
| verified `AnatomyPoser` | **109.79 mm median, 223.97 mm max** |

`project()` is relative to the plant's own reference bodies, which is correct for
a **material embedding** — forces, the surface binding and the world exchange all
apply loads in a frame that must not move when a caller asks for an initial pose.
As a **picture** it means the anatomy is drawn where the atlas left it while the
mechanical body sits 88.5 mm lower and walks away underneath. The 109.79 / 223.97
pair is independently recorded in `docs/DISCONNECTS.md` as the model's default
pose against the registered rest pose, which is the known answer that says the
poser is the one telling the truth.

So the poser is carried **beside** `entities`, never instead of it — swapping it
in would move the frame forces are applied in. `ArticulatedBodyPlant(display_pose=
'opensim'|'anatomical')` adds `frame['display_pose']`, and the session accepts
`display_pose`.

**The frame ships 22 segment motions, not 3,995 entity poses.** A rigid binding
gives every entity on a segment the same motion, so the per-entity form was 3,995
copies of 18 distinct matrices — **1,159,032 bytes a frame**, which is not a
per-frame payload. The compact form is **8,111 bytes** (8,077 before 7ded91b added the `plant` field), a 143× reduction, and the
static entity→segment map goes once as `display_pose_map.json`.

`scripts/verify_display_pose.py`, all passing: reconstruction from the 22 motions
reproduces the per-entity pose with **max |difference| = 0.000e+00 m** (exact, not
approximate); the same state twice gives bit-identical motions (the CLAUDE.md
call-it-twice check); every motion is a proper rotation to 1.11e-15; the 5 unbound
entities are listed with a reason on every frame.

*Still true, and unchanged by this:* one rigid segment per entity, no soft-tissue
deformation, no volume preservation, no sliding, and the 220 entities that straddle
a joint still tear at it.

### 4.2 Everything the poser cannot do

From `docs/DISCONNECTS.md` §1, and unchanged by calling it:

* one rigid segment per entity; no soft-tissue deformation, no volume
  preservation, no sliding. **220 entities straddle a joint and tear at it.**
* the skin is a linear blend, not a skin model.
* the lymphatic network graph is posed by nothing — no single rigid segment
  carries it.
* **nothing flows back**: the anatomy feels no contact and applies no force.

### 4.3 The default view is a recording, and it barely moves

`/api/body/trajectory` with no run serves `data/derived/canonical/trajectory.json`.
That trajectory moves a maximum of **6.35 mm over 30 s** — it is breathing and
perfusing, not moving. The transport now labels recording vs live correctly
(`app/src/transport.js`), which is a real fix; the remaining issue is that the
thing a visitor sees first is the least dynamic artefact in the repo.

---

## Tier 5 — physiology

BioGears runs as the unmodified upstream engine at a pinned revision, which is
the right call. Its own declared limitations (`ihm/native/__init__.py`) are:

* **"No explicit posture action found in upstream source."** The supine body and
  the standing body are thermally and hydrostatically the same body.
* radiation uses an upstream **standing** area factor even supine.
* `EvaporativeHeatLoss` repeats convective heat loss because of an upstream
  telemetry defect, and is **excluded from energy accounting**.
* `'experimental_validation': False`, uncertainty
  `'unquantified; deterministic source simulation, not calibrated human evidence'`.
* a stopped intervention initiates recovery but does not guarantee return to
  baseline.

The mechanics→physiology exchange is a **first-order lag with
`METABOLIC_EXCHANGE_TAU_S = 2.` seconds**, declared in the source as *"an
explicit engineering choice, not an identified substrate kinetic"*. Physiology→
mechanics is 3.4 above: absent.

*To close:* a posture term is upstream work or a fork. The evaporative defect
should be fixed in the variant tree rather than excluded.

---

## Tier 6 — a second engine, and dead paths

`ihm/assembly/interactive_scene.py` is a **separate, hand-rolled reduced physics
engine** served at `/api/scene/sessions`, with its own declared scope:

    body_mechanics    'Canonical linked rigid translations and affine soft regions'
    body_rotations    'Reference orientations constrained'     <- the body cannot rotate
    body_environment  'no body-surface mattress or floor contact solve'
    body_object_contact  False                                 <- objects pass through
    clothing_contact     False
    physiology_feedback  False

The shipped UI hardcodes `bodyEndpoint('embodied')`
(`app/src/scene-interaction.js:63`), so this engine is **not** what a visitor
drives — but it is live on the API under a name a caller would reasonably
mistake for the body, and it would answer.

Similarly duplicated: **hair has two solvers**, `ihm/assembly/hair_dynamics.py`
(server, opt-in, frozen candidate) and `app/src/hair_dynamics.js` (client worker,
what the app actually runs, default off).

*To close:* delete the reduced engine or gate it behind a name that cannot be
mistaken for the body. Pick one hair solver.

### 6.1 Done, 18 September 2026 — the name now carries the truth

The reduced engine was **renamed, not deleted**, because deleting it would take
working, verified code and because the module is load-bearing for something else
(6.2). What is now true:

**The path says what it is.** The engine is served at
`/api/reduced-kinematics/sessions`. Nothing in that name can be read as the body.

**The old path is retired, not aliased.** `GET`/`POST` on `/api/scene/sessions`
and on every `/api/scene/sessions/{id}/…` sub-path return **410** with a message
that names the real endpoint, states the three disqualifying facts, and gives the
new reduced path for a caller who genuinely wanted it
(`RETIRED_SCENE_SESSIONS` in `ihm/app/__init__.py`). A silent redirect was
rejected deliberately: it would leave the caller believing the old name meant the
body. The 410 is answered before the `Content-Type` check, so a caller that gets
the body of the request wrong is still told what it actually asked for.

**Every response from the engine is disclosed at the top level.** `NOT_THE_BODY`
in `ihm/assembly/interactive_scene.py` is spread into every payload the engine
emits — create, step, insert, close, snapshot, and the immutable event journal —
**first and last**, so it reads before anything else and no payload key can
overwrite it:

    is_body_simulation    false
    engine                'reduced-kinematics'
    use_instead           '/api/embodied/sessions - the OpenSim/Simbody + BioGears body …'
    body_rotations        'Reference orientations constrained; this body cannot rotate'
    body_object_contact   false
    body_environment      'No body-surface mattress or floor contact solve'

The last three were already in `scope`; they are promoted because a caller who
has to open `scope` to find out has already been misled. `app/src/embodied-live.js`'s
`frameScope` now leads with that disclosure for any frame carrying
`is_body_simulation: false`, ahead of the scope prose.

**`/api/scene/catalog` deliberately keeps its name and is deliberately NOT
disclosed.** It is the environment/object/scene tile catalogue, not a physics
engine, and the **live embodied body draws its environment ids from it** — the
`environments` list is `ENVIRONMENTS` out of the reduced module. Stamping that
response `is_body_simulation: false` would be a false statement about the
catalogue. What it does instead: `ihm/app/scenes.py` rewrites the 11 derived
`selection[].endpoint` entries that still quote `POST /api/scene/sessions` to the
new path, marks each `is_body_simulation: false` with a `use_instead`, and adds a
note saying the rewrite happened and that neither path is the body.
`scripts/build_environment_catalogue.py` still emits the old string (it is
outside this change's territory); rebuilding it to emit the new one removes the
rewrite.

Verified: `python -m py_compile` on all three edited Python files (`ast.parse`
is not enough here); `scripts/verify_interactive_scene.py --light` passes its 13
checks unchanged; `cd app && npm test` 141/141; and a live server smoke test
confirmed 410 on all four retired spellings, 201/200 with the disclosure present
on create/step/current/close, and the catalogue rewrite.

### 6.2 What still needs the reduced module, and whether it can be deleted

**The HTTP session engine has no remaining caller. The module does.**

Nothing drives `/api/reduced-kinematics/sessions`: the shipped UI hardcodes
`bodyEndpoint('embodied')` (`app/src/scene-interaction.js:63`), and
`bodyEndpoint('reduced')` in `app/src/embodied-live.js` is reachable from no call
site. But `ihm/assembly/interactive_scene.py` is imported for three things that
are **not** the session engine:

| importer | what it takes | still needed |
|---|---|---|
| `ihm/app/scenes.py` → `/api/scene/catalog` | `SceneSessions.catalog()` → `ENVIRONMENTS` (gravity, axis, plane, supports), which the tile grid and the **live embodied** environment ids rest on | **yes** |
| `scripts/build_environment_catalogue.py:1062` | `ENVIRONMENTS`, and the `Sphere` defaults it cites as provenance for object mass/radius | **yes** |
| `scripts/verify_loaded_source.py` | `_loaded_source`, the loaded-bytecode/source identity check | **yes** |
| `scripts/verify_interactive_scene.py` | `Sphere`, `InteractiveScene`, `SceneSessions` — verification of the engine itself | only while the engine exists |
| `scripts/scene_object_contact.py` | names the module in prose only | no |

So: **the engine cannot be deleted outright today, and the file certainly cannot.**
The proposal, in order, for whoever takes it:

1. Move `ENVIRONMENTS` and `vector()` into a module that is not a physics engine
   (`ihm/assembly/environments.py`), and `_loaded_source`/`_code_hash` into a
   `ihm/assembly/source_identity.py`. Both are pure data/utility and neither
   belongs to the reduced solver.
2. Delete the HTTP surface: the two route blocks in `ihm/app/__init__.py`, the
   `reduced` branch of `bodyEndpoint`, and `server.scenes` as a session manager
   (`scenes.catalog()` becomes a direct read of the relocated `ENVIRONMENTS`).
   Keep the 410 on `/api/scene/sessions`.
3. Delete `InteractiveScene`, `Sphere`, `SceneSessions` and
   `scripts/verify_interactive_scene.py` together. The retained audit at
   `data/derived/audits/interactive-scene-wdowc0_3/verification.json` is the
   record of what it did; nothing in the live body depends on it.

Steps 1–2 are the ones worth doing whether or not step 3 ever happens: they end
the situation where the environment catalogue the real body uses is owned by a
file whose main export is a different physics.

---

## Tier 7 — appearance authored in the app

All of this is labelled in source. It is listed because "100% authentic" has to
reach it eventually, not because anything here is hidden.

| what | status |
|---|---|
| **colour, default palette** | `didactic`, **17 of 17 roles synthesized** |
| **colour, `realistic` palette** | 74 roles: **5 measured, 5 transferred, 64 synthesized** — and only **3** of the 5 measured meet the stated bar (§7.1) |
| skin deformation | linear blend skinning; *"no FEM or cloth"* (`surface-binding.js:3`) |
| garments | *"a geometric pattern/ease prior, not a cloth contact solve"* (`clothing.js:6`) |
| tissue relief/grain | procedural fbm; grey only, multiplies the structure's own colour |
| environment, scenery | catalogue geometry + procedural materials |
| upright/supine buttons | **rigid display rotations**, not simulation posture (`docs/APP.md`) |
| orientation gimbal | decimated from the body's own envelope — authentic |

### 7.1 The colorimetry campaign ran, and the literature is mostly not there (18 Sep 2026)

`docs/TISSUE_COLORIMETRY.md`, commit `87fb93d`. A value is `measured` only if it is a
measurement of that tissue in vivo or fresh, with a stated illuminant and observer,
converted by a stated transform. Held to that, the campaign moved **one** role:
`mucosa_pharyngeal` synthesized → transferred (Hosoki 2007 buccal mucosa, n = 62, D55).
Three more rest on better sources at the same tier (palmoplantar skin now Wang/Luo
2017 with D65/10° stated; palatal and lingual mucosa off the wrong donor). Adipose,
skeletal muscle and hair each had a candidate that was **read and rejected** —
carcass fat altered by a steam cabinet, unbloomed pork spanning 13 L\* units with no
observer, hair with no observer.

**The conversion is now code with known answers** (`ihm/colorimetry.py`,
`scripts/test_colorimetry.py`): Pascale 2006's ColorChecker tables reproduce to
**6/65535 counts at 16 bits and 1/255 (0.50 ΔE\*ab) at 8**; the Bradford matrices match
Lindbloom's published ones to 4.8e-8. Its first run failed at 204 counts because it
used Pascale's printed γ = 0.42 where his table was computed with 1/2.4; the
tolerance was not moved, and the 0.42 version stays in the test as a check that must
fail.

**Two of the five `measured` roles do not meet the bar.** Gingiva (Ho 2015) and nail
(Horibata 2025) state neither illuminant nor observer — both full texts re-read.
They are left at `measured` under the previous build's assumed D65/2°, and every
entry now records whether its illuminant was stated or assumed. **To the stated bar,
measured is 3.** Demoting them is a definitional call, not a new measurement, so it
is left to the owner.

**No usable published colorimetry exists** for the viscera, the glands, the male
reproductive tract, the heart, vessels and blood, the nervous system, bone,
cartilage, tendon, ligament, fascia, adipose, or any of the eight eye roles. That is
a finding about the literature, not a gap in the search. `docs/TISSUE_COLORIMETRY.md`
lists nine library requests that could move the count, led by a 1993 gallbladder and
bile-duct spectral study and a 2025 thyroid/parathyroid study whose spectra are
available from the authors.

*Still to close:* the skin blend closes with 4.2. Garments need a real cloth solve or
an honest label in the UI rather than only in the source.

---

## The register, ordered by what to do next

Rewritten 18 Sep 2026, end of day. "Done" means verified by a re-run here, not by a
worker's report.

**Closed today**

| | |
|---|---|
| §1.1 | plant kwargs forwarded — stops, real segment contact, 66 tissue elements, reachable from a live session and the UI |
| §1.3 | the `CoordinateLimitForce` damping was **57.3×** its declared value; fixed with the plant bit-identical |
| §4.1 | the verified `AnatomyPoser` is called; 22 segment motions, 8,111 bytes/frame, exact reconstruction |
| §6.1 | the reduced engine cannot be mistaken for the body; old path 410s |
| §0.3 | one source of truth for stature and mass — and `77.6122029` **has** a derivation: a BioGears patient's fed weight |
| §3.1 | a size-principle recruitment pool and the source's knee/hip reflexes, default unchanged |
| §3.4 | metabolic supply can limit force, opt-in, declared as accounting not physiology |
| §7.1 | the colorimetry campaign ran; the literature measures almost none of the body |
| — | the 25-body variant **runs end to end**: force frame, skin weights, display pose, all 23 checks |
| — | the skin sensor cap is the quadrature (**21,382**), not 128 |
| — | the bed skin foundation **loads again** — its manifest had been stale since 8 Sep |
| — | nociception: cutaneous and visceral transduction, every constant cited |
| — | soft tissue: a layer that deforms and pushes back on its bone (§2.1) |
| §2.1 | that layer is **in the plant's loop**: opt-in, off is bit-identical, and sub-cycling measured and refused |

**Open, in order**

1. **Finish the ankle-collapse mechanism.** The factor is the torso repartition
   (§0.2). Two mechanisms have been **excluded by purpose-built arms** — the subtalar,
   then the support plane. What survives: the proxy sphere radius (0.2577 → 0.3090 m,
   unchanged by pinning the floor) against the repartitioned inertia itself. Needs a
   per-body proxy-radius engine option, then `tweld` with the torso ball pinned at
   0.2577 m. **Queued for the next sequenced engine batch.**
2. **Make the soft-tissue coupling implicit, and affordable.** It is IN the loop now
   (§2.1), contact-gated, and free when unloaded. Two things it is not. Sub-cycling was
   measured against a per-step reference and **refused**: the reaction grows at
   1,473 N/s, so a 10% bar allows a 1.0 ms hold — a tenth of the plant step — and even
   at one solve per step an explicit reaction is 10x outside its own bar. That needs an
   implicit coupling, which needs a stiffness through the engine's force interface,
   which the interface does not take. And it costs 53.9 ms per 10 ms step on ONE
   14,874-DOF segment, which needs a compiled solver or a nonlinear reduced basis that
   passes K2/K3.
3. **The shoulder girdle.** MoBL-ARMS 4.1 is **non-commercial BSD-3 with two required
   citations** (established from primary sources), and is already published in a public
   repo — an owner decision, not an engineering one. The licence-clean alternative is
   Seth 2019 (CC BY 4.0), whose provenance and completeness are being established.
4. **Muscles for the joints that still have none.** The thoracic joint has no donor
   muscles (candidate attachments sit 32 mm above a joint centre itself uncertain to
   62 mm); the wrists wait on a forearm donor. The subtalar and neck are done, and the
   toe joint's geometry is now corrected in a variant (§0.2a).
5. **Decide the identified plant's toe joints** (§0.2a). Measured cost is negligible;
   the corrected foot is 27.8 mm longer, so adopting it means re-identifying the plant.
6. **A fatigue law.** Blocked on the engine: no published law runs on what BioGears
   exposes (no phosphate, pH or PCr). Until then supply-limiting is bookkeeping.
7. **Pain as a reward.** Transduction exists; the reward decision is deliberately not
   made, and a damage stimulus **during motion** does not exist — contact binds only to
   the supine reference pose.
8. **Posture in physiology**, and the upstream evaporative defect (Tier 5).
9. **Deformable skin with self-contact**, sliding for fascia and bursae, muscle as a
   volumetric activation-coupled solid (§2.1).
10. **Colorimetry:** 64 roles remain synthesized; progress needs the nine library
    requests, not more searching (§7.1).
11. **Two bodies of different stature and mass** stay two bodies (§0.3). Unifying them
    is a modelling decision; the seam is now declared in one place.

Items 1–5 are the body. Items 6–11 are the long tail. The milestone — the body pushing
itself up off the floor — is gated on 3 and 4.
