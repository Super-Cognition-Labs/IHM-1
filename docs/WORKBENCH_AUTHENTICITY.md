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

**What is actually missing, which is narrower and sharper:** no shoulder girdle
(no scapula, no clavicle), no deltoid, rotator cuff, pectoralis or latissimus; no
forearm rotation (the Arm26 donor fuses it — `pro_sup` has a 2.16 mm moment arm);
no wrist, no hand, no neck. The shoulder is three biceps and three triceps heads
per side. From the model's own masses, the 21.1 N·m shoulder-flexion bound holds
the arm's 41.2 N at 0.51 m — it can hold an arm out; **nothing here bounds
pressing the trunk up**, which is what the self-righting milestone needs.

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
  ankle? FAILED.** It still leaves its range by 1.140 rad right and 0.871 left,
  down from 1.356 and 1.352, against a 0.220 bar. So the missing arm was at most
  part of the cause, and the cause is now **open**. The run that would separate it
  — re-welding the subtalar alone — has not been made.
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

These are the toe extensors and flexors, and their real geometric arm about the
toe joint reaches 25.3 mm; the fit reads exactly zero, bit-identical on a repeat
call. Both toe joints are therefore unactuated hinges held only by passive terms.
Left untouched in the base plant deliberately — changing it changes every result
already recorded on it. Every foot contact, gait and crawl measurement on this
plant was made with toes that no muscle could move; that belongs beside any result
that depends on push-off.

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

### 1.3 A unit defect found by wiring this, in the shared engine

`scripts/native_mechanical_stream.cpp` builds each `CoordinateLimitForce` by
converting the limits, the stiffness and the transition width from radians to
degrees — and passing `damping` straight through. OpenSim's own header declares
that property as **`Nm/(degree/s)`** for a rotational coordinate
(`CoordinateLimitForce.h`, the `damping` property). So the 1.5 every stopped run in
this repo declares is applied as **1.5 × 180/π = 85.94 N·m·s/rad, 57.3× the
number**.

**The value has deliberately NOT been changed.** `scripts/crawl.py`'s measured
stiffness sweep, every stopped crawl, and the tissue-mechanics arms were all run
with this damping; changing the constant to "fix" the units would silently make a
stopped live body a different object from every stopped run already recorded, and
would invalidate the table this register quotes. What changed is that the field is
now named for what the engine actually consumes
(`damping_nm_per_deg_per_s`), and every resolution carries the effective radian
value and a note saying where the conversion is missing. The wire format is
byte-identical — `plant_options.joint_stops()` still reproduces
`crawl.joint_stops()` exactly.

Fixing the conversion is an engine change that invalidates every prior stopped run.
It should be done deliberately, with those runs redone, not as a side effect of
this work.

**And it cannot be done concurrently with anything else.** The engine build's
`manifest.json` hash-verifies `scripts/native_mechanical_stream.cpp` on every
single load — `raise ValueError('Stale native mechanical build: ' + path)`. So
editing that file does not just change the next build; it **immediately makes
every existing build stale and refuses to start any engine at all**, including
ones already relied on by work in flight. The fix therefore has to be sequenced
alone: edit, rebuild, verify the new build, flip
`data/runtime/mechanical-stream/latest.json`, then redo the stopped runs. Doing it
beside other work takes the plant away from that work without warning.

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

* **at most 128 skin sensor points** may be selected per session
  (`mechanical_stream.py:110`). A whole body's cutaneous afference is 128 numbers.
* no self-contact, no skin-on-skin, no sliding, no volume preservation.

*To close:* Tier 1's wiring gives upright the skin layer. Beyond that: a
deformable skin with self-contact, and an afferent allocation that is not capped
at 128 points.

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

*To close, in order:* upper-body effectors (0.1) → a recruitment model that is
not a linear rate gain → reflexes beyond the ankle → metabolic supply as a force
limit rather than an exception.

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
per-frame payload. The compact form is **8,077 bytes**, a 143× reduction, and the
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
| **colour, `realistic` palette** | 74 roles: **5 measured, 4 transferred, 65 synthesized** (each cites its source or says it has none) |
| skin deformation | linear blend skinning; *"no FEM or cloth"* (`surface-binding.js:3`) |
| garments | *"a geometric pattern/ease prior, not a cloth contact solve"* (`clothing.js:6`) |
| tissue relief/grain | procedural fbm; grey only, multiplies the structure's own colour |
| environment, scenery | catalogue geometry + procedural materials |
| upright/supine buttons | **rigid display rotations**, not simulation posture (`docs/APP.md`) |
| orientation gimbal | decimated from the body's own envelope — authentic |

*To close:* colorimetry for the remaining 65 roles is a literature campaign, not
a code change. The skin blend closes with 4.2. Garments need a real cloth solve
or an honest label in the UI rather than only in the source.

---

## The register, ordered by what to do next

Updated 18 Sep 2026. Items 1–3 are done and verified; the numbered order below is
what remains, with what each one now actually requires.

**Done**

1. ~~Forward the plant kwargs~~ — **done** (§1.1). Joint stops, real segment contact
   surfaces including both skin arms, and the admissible 66 tissue elements are
   reachable from a live session and from the workbench UI, with the bundle's own
   declared material and a disclosure on every frame.
   `scripts/verify_plant_fidelity.py`.
2. ~~Call `AnatomyPoser`~~ — **done** (§4.1). The frame carries 22 segment motions
   that reconstruct 3,995 entity poses exactly, at 8,077 bytes against 1,159,032.
   `scripts/verify_display_pose.py`.
3. ~~Delete or rename the reduced engine~~ — **done** (§6.1). It is
   `/api/reduced-kinematics/sessions`, the old path is a 410 naming the body, and
   every payload leads with `is_body_simulation: false`. Deleting the module
   outright is blocked only by three non-solver things it owns; §6.2 has the
   sequence.

**Next, in order**

4. **Give the existing joints muscles that can drive them.** This displaced
   "upper-body muscles": the plant already has 92–98 muscles with 12–18 above the
   pelvis (§0.1, withdrawn and corrected). What it does not have is a shoulder
   girdle, and — measured on the new spine variant — **no muscle has a moment arm
   about any newly freed coordinate, exactly 0**, because the 80 source paths are
   polynomials fitted in the coordinates that existed when they were fitted.
   Un-welding a joint does not give a muscle a way to move it.
5. **The shoulder girdle.** MoBL-ARMS 4.1 is already on disk; the blockers are
   girdle bodies, a one-time torso mass partition, a coordinate map and an
   **unresolved licence** that should be settled first
   (`docs/UPPER_BODY_ACTUATION.md`).
6. **Fix the `CoordinateLimitForce` damping conversion** (§1.3) and redo every
   stopped run. 57.3× is not a rounding error, and it is load-bearing for every
   result that used a stop.
7. **Extend the segment binding past 22 segments** so the display can follow the
   spine variant's new bodies; today the anatomy still rides `torso`.
8. **One stature and one mass** (§0.3).
9. **Recruitment instead of a rate gain; reflexes past the ankle** (§3.1).
10. **Deformable skin with self-contact; afference past 128 points** (Tier 2, §4.2).
11. **Metabolic supply as a force limit rather than an exception** (§3.4).
12. **Posture in physiology; fix the evaporative defect** (Tier 5).
13. **Colorimetry for the remaining 65 palette roles** (Tier 7).

Items 4–7 are the body. Items 8–13 are the long tail.
