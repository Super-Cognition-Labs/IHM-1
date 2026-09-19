# What actuates the upper body, measured

`docs/WORKBENCH_AUTHENTICITY.md` item 0.1 says *"There are no muscles above the
pelvis"* and names the plant's 13 `CoordinateActuator` torque motors as what drives
the torso and arms. **That is true of the source model and false of the plant that
runs.** This file records what the live plant actually carries, measured on
18 September 2026, and what is genuinely missing — which is a narrower and more
specific thing than 0.1 claims.

Everything below came from the native engine, not from reading XML. Scripts:

```sh
OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
  .venv/bin/python scripts/measure_upper_body_actuation.py \
    --registration <manifest> --output <fresh dir>

OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
  .venv/bin/python scripts/probe_upper_body_muscle_drive.py \
    --group <elbow_r|triceps_r|shoulder_r|lumbar_ext> --output <fresh dir>
```

One native session at a time.

---

## 1. The premise of 0.1 is about a file nobody loads

Three different bodies are reachable, and the bare source model is the only one
without upper-body muscle. Measured muscle counts from `snapshot()['muscles']`:

| plant | muscles | above pelvis | who loads it |
|---|---:|---:|---|
| `subject_walk_scaled.osim`, no augmentation | **80** | 0 | nothing on the live path; only bare `NativeMechanicalStream(...)` calls that pass no registration |
| `data/derived/mechanics/whole_body_arm26_v2` | **92** | **12** | **the `EmbodiedRuntime` default** (`embodied.py:365`, `augmented_registration or 'data/derived/mechanics/whole_body_arm26_v2/registration.json'`) |
| `data/models/engineering_stance_v1` | **98** | **18** | **forced** for every controller in `STANCE_KINDS` (`embodied.py`, `controller_selection.py`) |

The 12 are `arm26_{TRIlong,TRIlat,TRImed,BIClong,BICshort,BRA}_{r,l}`, Thelen2003,
registered from the Arm26 donor. The further 6 are
`gait2392_{ercspn,intobl,extobl}_{r,l}`, a Gait2392 trunk insert.

And the torque motors are not driving anything. `NativeMechanicalStream.advance`
documents that the `CoordinateActuator` ports *"start at zero"* and stay at the
previous command; `coordinate_actuation` is passed by exactly four offline scripts
(`crawl.py`, `walk_gait.py`, `scene_object_contact.py`, `measure_tissue_mechanics.py`)
and by **nothing** in `ihm/app/` or `ihm/assembly/`. So on the live
path all 13 torque ports sit at zero for the whole session and the upper body is
moved by muscle or by nothing.

**So 0.1's sentence should read: the plant has 12 upper-body muscles by default and
18 under the stance controllers, covering the elbow, three shoulder axes and three
lumbar axes — and it has none at all for the wrist, the hand, the neck (there is no
neck), the shoulder girdle (there is no scapula or clavicle) or forearm rotation.**

---

## 2. Moment-arm census: which coordinates a muscle actually crosses

`NativeMechanicalStream.moment_arms` queried against every rotational coordinate,
for all muscles, at the plant's own initial state.
`+bound`/`−bound` are `Σ Fmax·|r|` over the muscles pulling each way.

**92-muscle default plant** (`whole_body_arm26_v2`, zero arm pose):

| coordinate | muscles | +bound N·m | −bound N·m | authority | torque port |
|---|---:|---:|---:|---|---|
| arm_flex_{r,l} | 3 | 18.7 | 31.8 | bidirectional | shoulder_flex |
| arm_add_{r,l} | 3 | 37.0 | 10.7 | bidirectional | shoulder_add |
| arm_rot_{r,l} | 3 | 5.6 | 2.9 | bidirectional | shoulder_rot |
| elbow_flex_{r,l} | 6 | 22.9 | 43.7 | bidirectional | elbow_flex |
| pro_sup_{r,l} | 2 | 2.3 | 0.0 | **unidirectional, 2.16 mm arm** | pro_sup |
| lumbar_extension | **0** | 0 | 0 | **none** | lumbar_ext |
| lumbar_bending | **0** | 0 | 0 | **none** | lumbar_bend |
| lumbar_rotation | **0** | 0 | 0 | **none** | lumbar_rot |

**98-muscle stance bundle** (`engineering_stance_v1`, at its identified pose):

| coordinate | muscles | +bound N·m | −bound N·m | authority |
|---|---:|---:|---:|---|
| arm_flex_{r,l} | 3 | 21.1 | 31.0 | bidirectional |
| arm_add_{r,l} | 3 | 38.1 | 9.5 | bidirectional |
| arm_rot_{r,l} | 3 | 5.1 | 3.7 | bidirectional |
| elbow_flex_{r,l} | 6 | 44.7 | 52.4 | bidirectional |
| lumbar_extension | 6 | 213.4 | 208.1 | bidirectional |
| lumbar_bending | 6 | 235.7 | 235.7 | bidirectional |
| lumbar_rotation | 6 | 80.1 | 80.1 | bidirectional |
| pro_sup_{r,l} | **0** | 0 | 0 | **none** |

Reports: `data/derived/upper-body-actuation-arm26v2/report.json`,
`data/derived/upper-body-actuation-stance98/report.json`.

### The bound is an upper bound and only its zeros are decisive

`Fmax·|r|` at one pose ignores force–length, force–velocity, tendon state and
activation dynamics, and a per-axis sum is not a tension-feasible 3-D torque cone.
`docs/research/LUMBAR_SHOULDER_MUSCLE_COVERAGE.md` warns about exactly this
multiplication. **Read the zeros, not the magnitudes.** A column of zeros means no
path crosses the joint at all, which no modelling assumption can rescue.

### Two controls, both passed

* **Idempotence.** `moment_arms` was called twice at the identical state and the
  two dictionaries compared bit-equal, in both runs. A native query that is not a
  function of its arguments produces confident wrong conclusions (IHM-1 `CLAUDE.md`,
  *"Call it twice at the same input"*).
* **Known answer.** The ankle must read as a muscle-driven ankle. It does:
  `soleus_r` returns a plantarflexion moment arm of **−49.7 mm** (92-plant) and
  **−47.0 mm** (98-plant), inside the 40–60 mm the anatomy requires, with 11 muscles
  crossing each ankle. The lumbar arms also reproduce
  `LUMBAR_SHOULDER_MUSCLE_COVERAGE.md`'s independently computed table exactly —
  `ercspn` 42.69 mm, `intobl` −52.82 mm, `extobl` −62.79 mm — so the trunk transfer
  landed where its own audit said it would.

### `pro_sup` is effectively unactuated, and the two runs show why

The 92-plant reports 2 muscles at 2.16 mm; the 98-plant reports none. Both are
right: the moment arm is pose-dependent and the stance bundle starts at
`pro_sup = 0.164 rad`, where it falls under the 1 mm floor. A real biceps
supination moment arm is 15–20 mm. 2.16 mm is a **consequence of the Arm26 donor
fusing ulna, radius and hand into one body**; `UPPERBODY_EFFECTOR_REGISTRATION.md`
records that the biceps insertion was then assigned to the native radius by hand.
Treat forearm rotation as having no muscle actuation.

---

## 3. The upper body moves under muscle drive alone

A moment arm says a path crosses a joint; it does not say the joint moves. Each run
below holds the stance bundle's own 98 equilibrium excitations, drives one agonist
group to excitation 1.0, integrates 40 steps of 10 ms, and reports the excursion as
a **paired difference against an identical null run** with the equilibrium unchanged.
`coordinate_actuation` is never passed, so all 13 torque ports are at zero.

| drive | watched | driven Δ | contralateral Δ | ratio |
|---|---|---:|---:|---:|
| BIClong+BICshort+BRA (r) | elbow_flex_r | **+105.98°** | +0.10° | 1015× |
| TRIlong+TRIlat+TRImed (r) | elbow_flex_r | **−32.82°** | +0.27° | 120× |
| BICshort+TRIlong (r) | arm_add_r | **+10.38°** | −0.91° | 11× |
| ercspn (bilateral) | lumbar_extension | **+120.10°** | pelvis_tilt −6.70° | 18× |

The null arms moved 9×10⁻⁹ to 2×10⁻⁸ rad over the same interval, so the equilibrium
is genuinely static and the excursion is the drive. Agonist and antagonist move the
elbow in opposite directions with the right signs. Reports under
`data/derived/upper-body-drive-*/report.json`.

**This is the answer to "can a muscle move the upper body": yes, measured, with
every torque motor at zero.**

### Two cautions that ride with those numbers

* **The lumbar drive blew through the declared range and nothing stopped it.**
  `lumbar_extension` reached **2.096 rad = 120.1°** against a declared
  `range` of ±1.5708 rad with `clamped=true` — **30.1° outside**, because no
  `coordinate_limits` were passed and the model holds zero `CoordinateLimitForce`
  (`docs/NATIVE_JOINT_LIMITS.md`, IHM-1 `CLAUDE.md`). The excursion demonstrates
  torque authority and is **not** a physiological range of motion. This is exactly
  the Tier-1 `coordinate_limits` wiring the register already asks for.
* **The shoulder coordinates declare a meaningless range.** `arm_flex`, `arm_add`
  and `arm_rot` each carry `range = ±10 rad` — **±573°**. There is nothing to clamp
  them to even when clamping is switched on. Any self-righting search will find this.

### Scale: these muscles can hold an arm, not press a body up

The model's own right arm is **4.196 kg** (humerus 2.302, ulna 0.688, radius 0.688,
hand 0.518), so **41.2 N**. The measured shoulder-flexion bound of 21.1 N·m is
41.2 N at 0.51 m — comfortably more than the arm's own COM demands, and an upper
bound. Pressing the trunk up off the floor is a demand of a different order and
**nothing here bounds it**. Do not read §3 as "the body can push itself up".

---

## 4. What is missing, precisely

| capability | status |
|---|---|
| elbow flexion / extension | **muscle, bidirectional, verified moving** |
| shoulder flexion / adduction / rotation | **muscle, bidirectional** — but from 3 muscles per axis, all of them biceps/triceps heads |
| lumbar extension / bending / rotation | **muscle** in the 98-plant only; **absent** in the 92-muscle default |
| forearm pro/supination | **none** (2.16 mm at best; donor fuses the forearm) |
| deltoid, rotator cuff, pectoralis, latissimus | **none in any plant that existed when this was written**; present in `data/models/shoulder_girdle_v1` (§11) |
| shoulder girdle | **no scapula and no clavicle body** in the plants above; both, with sternoclavicular, scapulothoracic and acromioclavicular articulation, in `data/models/shoulder_girdle_v1` (§11) |
| wrist, fingers | **none**; `radius_hand_{l,r}` is a `WeldJoint` |
| neck, head | **no body at all** |

The honest one-line version: **the arms have an elbow and a crude sagittal
shoulder, the trunk has an optional three-axis one, and everything that makes a
shoulder a shoulder is missing.**

*Amended 18 Sep 2026.* That is still true of every plant in the table above. It
is no longer true of the repository: `data/models/shoulder_girdle_v1` carries a
scapula, a clavicle, the three girdle articulations and thirty donor muscles per
side, built from the licence-clean donor rather than from MoBL-ARMS. §11, and
`docs/SHOULDER_GIRDLE.md` in full.

---

## 5. What is on disk

A census of every `.osim` under `data/raw/` (215 files) plus `data/research/`:

| model | where | muscles | girdle | usable? |
|---|---|---|---|---|
| Arm26 | `data/raw/anatomy/opensim-models/source/Models/Arm26/` | 6 Thelen2003 | no | **already registered and running** (the 12) |
| Gait2392 trunk | `.../Models/Gait2392_Simbody/` | 92 Thelen2003 | no | **already registered** (the 6) |
| **MoBL-ARMS 4.1** (Saul/Murray 2015) | `data/research/shoulder_complement/MOBL_ARMS_41.osim` | **50 Millard2012** | **thorax, clavicle, scapula** | **acquired, extracted, not installed** |
| Thoracoscapular shoulder (Seth 2019) | `data/raw/mechanics/opensim-core/OpenSim/Tests/shared/ThoracoscapularShoulderModel.osim` | 33 Millard2012 | **clavicle + scapula, 4 scapular coordinates** | **examined, licence established, and INSTALLED** — `data/models/shoulder_girdle_v1`, §11 |
| `PushUpToesOnGroundWithMuscles.osim` | `data/raw/mechanics/opensim-core/OpenSim/Simulation/tests/resources/` | **100 Schutte1993_Deprecated + 54 Thelen2003**, 81 bodies, bilateral DELT1-3/SUPSP/INFSP/SUBSC/TMIN/TMAJ/PECM1-3/LAT1-3/CORB | **yes, full** | **deprecated muscle law and `<credits>Model authors names..</credits>` — a placeholder. Corroborating routes only; not a citable source** |
| Rajagopal2016 / RajagopalLaiUhlrich2023 | `.../Models/Rajagopal/` | 80–81 Millard2012, **0 above the pelvis** | no | 18 `CoordinateActuator`s incl. both wrists — **the same gap one joint further out, not a fix** |
| Hamner 2010 full body | `.../Models/Hamner/` | 93 Thelen2003 | no | lower limb + trunk only |
| WristModel (Gonzalez 1997) | `.../Models/WristModel/` | 25 Schutte1993_Deprecated | — | includes `ECU_pre-`/`post-surgery` alternatives; deprecated law |
| Neck3dof fixture | `.../OpenSim/Tools/tests/resources/` | 5 Schutte1993 | — | engine fixture, placeholder credits (`CERVICAL_MODEL_SOURCE_AUDIT.md`) |
| MASI, Mortensen2018 cervical | `data/research/cervical/` | 78 Thelen2003 / 72 Millard2012 | skull, clavicle, scapula | acquired; **geometry not acquired** (login wall) |

**Correction to a retained record.** `data/sources/sensorimotor/upperbody_sources.json`
still carries `unacquired_richer_source: {url: simtk.org/projects/upexdyn, status:
"web tool returned 403; no model bytes acquired; do not count as model data"}`.
That was true when written and is **stale**: `MOBL_ARMS_41.osim` was acquired on
2026-09-06 through the CEINMS-RT mirror and its 15 shoulder compartments are already
extracted to `data/research/shoulder_complement/shoulder_forces.xml`. The record
should be updated to point at the acquisition, with the mirror-vs-official caveat
that `data/sources/shoulder_complement.json` already carries
(`mirror_equals_official_release_verified: false`). *Done 18 Sep 2026 (§10.4).*

---

## 6. Verdict

The register offered three outcomes. The true answer is **(a) for the arms and the
trunk, and (b) for the shoulder** — and (c) is false: nothing needs acquiring.

* **(a) — already on disk and already wired.** 12 arm muscles are the live default
  and 18 run under every stance controller. They load, they integrate, they carry
  real Thelen states, and they move the joints they cross. §1–§3 are the evidence.
  Item 0.1's premise — that the upper body is torque-motor driven — does not
  describe the running plant.
* **(a), unwired until now.** The trunk muscles existed and the fail-closed loader
  could not read them. Fixed; see §7.
* **(b) — a suitable source exists and needs registration work.** MoBL-ARMS 4.1 is
  acquired, under a non-commercial licence whose application here is the owner's call
  (§10), and blocked on a shoulder
  girdle the plant does not have. Scoped in §8.
* **(c) — does not apply.** No acquisition is required for the shoulder. The one
  genuine acquisition gap left is *geometry*, for the cervical donors.

---

## 7. What was wired

`ihm/assembly/sensorimotor_catalog.py::whole_body_effector_catalog` accepted only
`ihm.upperbody-registration.v1`. Every 98-muscle body — `whole_body_lumbar_current`,
`engineering_stance_v1`, `engineering_supported_rest_v1` — declares
`ihm.lumbar-muscle-variant.v1` and was **rejected outright**:

```
 92 whole_body_arm26_v2            OK 92 rows   coverage OK 92
 98 whole_body_lumbar_current      REJECTED: Invalid registered effector manifest
 98 engineering_stance_v1          REJECTED: Invalid registered effector manifest
 98 engineering_supported_rest_v1  REJECTED: Invalid registered effector manifest
```

`ihm/assembly/embodied.py::_prepare_mechanical_registration` has accepted **both**
schemas since the variant existed. So the two readers disagreed about what a
registered body is, and the consequence was concrete: `build_peripheral_coverage`
could not audit the ports of the body the live stance path actually integrates — the
six trunk effectors had no coverage row anywhere.

The loader now accepts both schemas and, for a variant, additionally verifies the
two fields the variant schema adds: `insert_path`/`insert_sha256` is hashed, and
`base_model_path` is required to be one of the manifest's already-verified `sources`
rather than trusted from its own field. After the change all four load, with 92, 98,
98 and 98 rows and matching coverage reports; `lumbar_trunk` effectors now appear in
the audit.

It stays fail-closed. `scripts/verify_registered_effector_variants.py` pins it:
nine manifest mutations are each refused, including `base_model_path` swapped for
`whole_body_arm26_v2/subject_with_arms.osim` — **a real, parseable, genuinely
registered model that simply is not one of this manifest's sources**, so a loader
that checked only existence would pass it. A control that can pass for the wrong
reason is worse than none. Restoring the old one-schema tuple makes the new tests
fail, so the tests can fail.

```sh
OPENBLAS_NUM_THREADS=1 prlimit --as=2147483648 -- nice -n 10 \
  .venv/bin/python -m scripts.verify_registered_effector_variants   # 6 tests
  .venv/bin/python -m scripts.verify_upperbody_effectors            # 1 test, unchanged
  .venv/bin/python -m scripts.verify_peripheral_coverage            # 5 tests, unchanged
```

**What was NOT wired, deliberately.** The 98-muscle variant was not made the
default. Its own catalog rows say `native_control_ready: false` and
`default_excitation_assignment: null`, its manifest says `default_enabled: false`,
and §3 shows its trunk drive leaving the declared lumbar range by 30°. It is
reachable today as an explicit opt-in —
`EmbodiedRuntime.from_workspace(..., augmented_registration='data/derived/mechanics/whole_body_lumbar_current/registration.json')` —
and that is the right place for it until `coordinate_limits` are forwarded.

---

## 8. What installing the shoulder would take

Not an acquisition. A registration, and the mass partition is the hard part.

> **Written for MoBL-ARMS, and MoBL-ARMS is not the donor that was used.**
> §10 established that its SimTK notice restricts it to non-commercial use, so
> the girdle was built from Seth 2019 instead. Steps 1-5 below describe the job
> accurately and §11 records how each came out; step 6, the licence, is why the
> donor changed. Nothing in `data/models/shoulder_girdle_v1` comes from
> MoBL-ARMS.

1. **Girdle bodies.** MoBL-ARMS routes 44 fixed, 13 moving and 8 conditional path
   points across `thorax`, `clavicle`, `scapula` and `humerus`, over 23 wrap objects.
   The plant has no `clavicle` and no `scapula`. Either add them — with the donor's
   sternoclavicular/acromioclavicular couplers and its `shoulder0/1/2` decomposition —
   or weld them to `torso` and declare a fixed-girdle approximation in the UI, in
   which case pectoralis minor and every girdle-only path produce no joint torque.
2. **Mass partition, once.** Donor clavicle 0.156 kg and scapula 0.70396 kg per side.
   These must be **subtracted** from the torso's residual — mass, first moment and
   origin inertia, in a common frame, with a positive-definiteness check — not added.
   This has to be sequenced with the cervical/thoracic partition or the torso ledger
   is split twice. The donor's 0.0001 kg phantom bodies are numerical and must never
   become tissue.
3. **Kinematic reconciliation before any force.** The donor's default elevation is
   ~30° and its thorax sits ~−90° about y. The plant's three shoulder coordinates are
   not the donor's Euler decomposition. The orientation and angular-velocity Jacobian
   map must be defined and tested across the supported range including singularities.
4. **Reconcile with what is already there.** `arm26_TRIlong/BIClong/BICshort` already
   take registered `torso` origins. Adding DELT/PECM/LAT does not license duplicating
   or silently replacing them, and Millard (donor) and Thelen (Arm26) laws stay distinct.
5. **Gates.** Instantiate the unmodified donor first and reproduce its own path
   lengths, conditional transitions, wrapping and moment arms; then the registered
   variant against those; then a tension-feasible torque cone. `native_verified` stays
   false until each passes. Left side is a reflection **prior** and needs explicit
   validation of vectors, joint axes and wrap quadrants — not a name suffix.
6. **Licence.** `data/sources/shoulder_complement.json` records conflicting terms —
   SimTK noncommercial + BSD-3 wording, CEINMS mirror Apache-2.0, an older catalog
   listing MIT — deliberately unresolved, with `official_model_zip_held: false` and
   `mirror_equals_official_release_verified: false`. **Resolve the licence before
   anything from this donor ships**, not at integration time.
   *Primary-source findings, 18 Sep 2026: §10. The SimTK notice (BSD-3 restricted to
   non-commercial use, two citations required) is the operative statement; the donor is
   not cleared.*

`docs/research/SHOULDER_COMPLEMENT_SOURCE_DESIGN.md` holds the parameter tables and
the full gate list; this section is its summary with the plant-side blockers named.

**Ordering.** Steps 1–2 are the same work as 0.2's shoulder girdle, and the mass
partition is the same work as 0.3's single-mass reconciliation. Doing the shoulder
first and the girdle later means partitioning torso mass twice.

---

## 9. Changes this implies elsewhere

* `docs/WORKBENCH_AUTHENTICITY.md` 0.1 should be rewritten around §1 and §4: the
  claim is not "no muscles above the pelvis" but "no girdle, no forearm rotation, no
  wrist, no neck, and a shoulder built from three biceps/triceps heads per axis".
  Its "*to close*" note already points at `whole_body_arm26_v2` and
  `whole_body_effector_catalog` and was pointing at the right things.
* `data/sources/sensorimotor/upperbody_sources.json`'s `unacquired_richer_source`
  record is stale (§5).
* `data/derived/mechanics/whole_body_arm26_v2/registration.json` carries
  `native_verified: false` with the note that the mechanics owner must run a native
  load. That load has now run (§1–§3) and the promoted bundles built on it already
  carry `native_acceptance_complete: true`. The flag on the v2 manifest is stale;
  it is immutable by design, so the correction belongs here and in whatever
  supersedes it, not in an edit.
* The register's ordering puts "upper-body muscles" at item 4. On this evidence the
  arm and trunk halves are done, and what remains of 0.1 merges into 0.2's girdle
  and 0.3's mass reconciliation.

---

## 10. What licence governs the shoulder donors — primary sources, 18 September 2026

§8.6 said to resolve the MoBL-ARMS licence before anything from that donor ships. This
section records **what each primary source says, who issued it and what it covers**. It
is a statement of fact, not a legal opinion, and **it does not clear the donor to ship.**
Whether this programme's use falls inside the terms is the owner's decision; this is
the evidence for making it. Nothing was installed, no geometry was acquired, and no
simulation was run.

Every page below was read in full, not summarised by a tool. Pages already retained in
`data/research/shoulder_complement/` were re-fetched live on 18 September and compared.
The licence string on the SimTK download page is byte-identical in the 6 September copy
and the live page (sha256 of the string `847c2df4…0734`). The CEINMS `README.md` and `LICENSE`
on `main` are byte-identical to the pinned copies held in the repo.

### 10.1 MoBL-ARMS 4.1: five sources, and what each actually says

**(1) SimTK project `upexdyn` — the operative statement.**
The same licence string is attached to all three release packages at
<https://simtk.org/frs/?group_id=657>: the 2014-07-06 SIMM/OpenSim release, the
2016-03-11 "OpenSim 3.2+" release and the 2021-02-22 "OpenSim 4.1+" release
(`MobL_ARMS_OpenSim41_unimanual_tutorial.zip`). The last is the release whose name
`MOBL_ARMS_41.osim` matches, and the one the CEINMS README names as its source. The
project page carries the same string (`simtk_project.html`). Verbatim:

> Copyright (c) 2014-present, North Carolina State University, Northwestern University,
> Rehabilitation Institute of Chicago, Valparaiso University, Wake Forest University,
> Stanford University and VA HCS. All rights reserved.
>
> The MoBL-ARMS upper limb model has been open sourced solely for non-commercial
> purposes (including research, academic, evaluation and personal uses) under the BSD
> 3-Clause License below. By downloading or using this software, (1) you accept the
> terms and conditions of the aforementioned open source license, (2) acknowledge that
> your use of this software is non-commercial and commercial use requires a commercial
> license, and (3) accept that use of the model software must be acknowledged in all
> publications, presentations, or documents describing work in which the MoBL-ARMS upper
> limb model is used by citing the following work: [Saul et al. 2015, CMBBE 18:1445–58;
> McFarland et al. 2019, J Biomech Eng 141(5):051006]

That is followed by the standard BSD-3-Clause text: keep the notice in source
redistributions, reproduce it in the documentation of binary ones, and do not use the
holders' names for endorsement.

*Issuer and authority:* the project is run by the model's authors (team: Katherine Saul,
Wendy Murray), and the notice is written in the name of the seven institutions it lists
as copyright holders. It is the only statement found that comes from the rights holders.
*History:* the earliest Wayback capture that includes the licence text is
2018-09-12. It carries the same non-commercial notice, headed "Copyright (c) 2014-2015",
and it applies to both releases that existed then. Between the 2021-01-26 and 2021-04-12
captures three things changed: a copyright header was added at the top, the years became
"2014-present", and the McFarland 2019 citation was added. **The non-commercial clause
is the same in every capture.** No capture of the licence text from before 2018-09-12 was
found. The only earlier capture of the project is a 301 redirect from 2017.

**(2) The model file itself.** `MOBL_ARMS_41.osim` contains only
`<credits>Katherine R. Saul, Wendy M. Murray, Craig M. Goehler, Melissa Daly, Meghan E.
Vidt, Dustin L. Crouch</credits>` and `<publications>Comp Meth Biomech Biomed Eng
2014</publications>`. The file contains no licence, copyright or terms text anywhere.
Crouch is not among the authors SimTK lists, and nothing held says who made the 4.1
changes.

**(3) The originating paper.** Saul et al. 2015 (PMC4282829, held) states no terms. It
says only: *"Control inputs, simulation results, and the model itself will be publically
available via simtk.org ( https://simtk.org/home/upexdyn/ )."* It points to SimTK and
sets no terms of its own.

**(4) The CEINMS-RT mirror's Apache-2.0: the mirror does not claim it covers the model.**
The root `LICENSE` is Apache-2.0. It arrived in the repository's `Initial commit`
(`23f72c6ead`, 2022-09-01), before any MoBL file existed there. `MOBL_ARMS_41.osim` was
added later, in `5924cd622f` (2023-03-15, "Added the osim model."). GitHub's `license:
apache-2.0` field is GitHub's automatic detection of that root file, and the programme's
"Apache-2.0" reading came from that field. The mirror's own README says otherwise.
At the commit that added the model it said *"The MoBL license is non comercial (see
below)"* and gave the SimTK notice. Today it says *"CEINMS-rt is licensed under the
[Apache License](LICENSE)"* and then, under `### MoBL OpenSim model:`, reproduces the
SimTK notice verbatim. **The mirror never claimed to relicense MoBL-ARMS**, and CEINMS-RT
is not among the named copyright holders. It keeps the notice, which is what BSD-3 asks
of a redistributor.

**(5) The OpenSim catalog's "MIT" is a third-party wiki cell about the 2014 release.**
Confluence page 53090607 ("Musculoskeletal Models") has this row: *"Upper Extremity Dynamic
Model | … Katherine Saul, Xiao Hu, Craig Goehler Meghan Vidt, Melissa Daly, Anca Velisar,
Wendy Murray | Research-grade kinematics and dynamic simulation of shoulder and arm
movement. | MIT | July2014"*, and the licence cell links to opensource.org/licenses/MIT. The
page history (REST API) shows the row first appearing in **version 31, 2014-07-21**, edited
by an account displayed as "james". It is absent from version 30 (2014-03-05) and unchanged
through the current version 43 (2022-10-04). The row sits under the heading *"Models
contributed by members of the OpenSim community. These models are developed and maintained
by the authors listed, NOT the Stanford OpenSim team."* So the MIT cell is a catalog
editor's entry, not shown to come from the rights holders. It describes the July 2014
release (v31 adds "Compatible with 3.1") and was written 6½ years before the 4.1 release
existed. **Whether SimTK itself showed MIT in July 2014 cannot be established**, because no
capture exists. That question does not bear on the 4.1 release.

**No opensim-org distribution exists.** `opensim-models/Models/` holds no MoBL-ARMS, and a
GitHub code search for `MOBL_ARMS_41` (28 hits) returns no opensim-org repository. No
Stanford- or OpenSim-issued licence for this artefact exists to weigh.

### 10.2 Which statement governs, and how sure that is

**The SimTK notice governs, with high confidence.** It is the only statement issued in the
name of the copyright holders. It is the only one attached to the 4.1 release, and it has
been stable in substance for every year it can be observed. The Apache reading misreads
repository metadata that the mirror's own README contradicts. The MIT reading is a 2014
catalog cell about a different release, entered by someone other than the holders.

**What its text permits:** non-commercial use, which it says includes research,
academic, evaluation and personal use. It also permits redistribution and modification in
source or binary form for those purposes, provided the notice, conditions and disclaimer
travel with it.

**What its text does not permit:** commercial use without a separate commercial licence
(*"commercial use requires a commercial license"*); any use without citing Saul 2015 **and**
McFarland 2019 in every publication, presentation or document describing the work;
redistribution without the notice; use of the holders' names for endorsement. Nor does it
support labelling the donor, or anything derived from it, as Apache-2.0 or MIT.

**What stays ambiguous** (none of these is decided here):

* "Non-commercial" is defined only by its parenthetical. Whether this programme, its public
  site or any downstream counts as non-commercial is the owner's decision.
* A non-commercial restriction layered on BSD-3 is not an OSI open-source licence, whatever
  the notice calls it. The BSD text alone restricts nothing; the restriction sits in the
  preamble.
* The mirror's bytes are not verified equal to the SimTK 4.1 package, because SimTK
  requires a login. By its wording the notice applies to *"downloading or using this
  software"*, whatever route the bytes took.
* The terms before 2018-09-12 are unknown.
* Whether numerical parameters extracted from the model are covered is a legal question.
  This applies to `shoulder_forces.xml` and to the muscle and body tables in
  `data/sources/shoulder_complement.json`, and it is not answered here.

**A fact the owner needs, because it has already happened.** `MOBL_ARMS_41.osim`,
`shoulder_forces.xml` and `data/sources/shoulder_complement.json` (which carries donor
parameters) are tracked in git. They are on `origin/feat/integrated-human` of the
**public** repository `github.com/JacobFV/IHM-1`, first in commit `e99ad6e`. IHM-1 has no
repository licence. The same directory keeps the CEINMS `README.md`, which contains the
MoBL notice verbatim, and the CEINMS Apache `LICENSE`, which does not apply to MoBL. So the
donor is already being redistributed. Whether that redistribution is within the terms is
the owner's decision. It is recorded here so the decision is not made in ignorance of it.

### 10.3 The alternates in §5

**Thoracoscapular shoulder (Seth, Dong, Matias, Delp 2019).**
<https://simtk.org/projects/thoracoscapular> (group 1708) releases
`ThoracoscapularShoulderPaperMaterials.zip` (2019-06-28) under:

> Copyright (c) 2019, Stanford University and the authors. authors: Ajay Seth, Meilin Dong,
> Ricardo Matias, Scott Delp. This work is available under the Creative Commons Attribution
> 4.0 International Public License … You are free to: Share … Adapt — remix, transform, and
> build upon the material for any purpose, even commercially.

The paper's Data Availability Statement (Front. Neurorobot. 13:90, read from the NMBL PDF)
agrees: *"The model and simulation environment (OpenSim) are freely available, deployable,
and modifiable for any research or commercial use without restrictions from SimTK.org."*
The copy on disk is not that package. It is
`opensim-core/OpenSim/Tests/shared/ThoracoscapularShoulderModel.osim`, added to the
Apache-2.0 opensim-core repository by the OpenSim team in `62205cd879` (2021-03-12,
#2971). Its `<credits>` read *"Ajay Seth, Meilin Dong, Ricardo Matias, Scott Delp.
Parameters from van der Helm and Klein-Breteler"*, and it contains no licence text. Both
statements (CC BY 4.0 from the holders, Apache-2.0 on the repository) permit commercial use
with attribution. CC BY 4.0 also requires a link to the licence and a statement of changes.
*Open:* the on-disk bytes are not verified equal to the SimTK release. Its 7 `.vtp` meshes
are not on disk. The DSEM parameter sources' own terms were not examined. **Do not confuse
it with** Seth et al.'s 2016 *Scapulothoracic Joint* project (`simtk.org/projects/scapulothoracic`,
group 986). That is a different artefact with a different licence, **CC BY-NC 3.0**
("Copyright (C) 2015 Stanford University").

**The Holzbaur-derived fixture: `PushUpToesOnGroundWithMuscles.osim`.** Its
`<credits>` is the placeholder `Model authors names..`. Its `<publications>` ends:
*"NOTE: This model has been developed on the base of the 3DGaitModel2354.osim model
(developed by Delp S.L. et all) and the Stanford VA Upper Limb Model.osim model (developed
by Holzbaur KR et all). Not additional informations are present in this model that cannot
be found in the original ones. This model has been built by the author aiming only to
provide a starting model for a full body model."* That text and its citation list match
Andrea Menegolo's *Upper and Lower Body Model* (`simtk.org/projects/ulb_project`,
`ULB_Project_v02`). The match is textual; no recorded commit links the two. The ULB page
itself contradicts the "no additional information" claim, since it says the model *"comes
with assigned mass properties on the basis of the cited publications"*.

* The **ULB download page shows no licence at all**, and the catalog's licence cell for it
  reads `none` (Sep-11). The derivative's own author granted nothing that can be found.
* The fixture's only licence is opensim-core's repository-level Apache-2.0. That repository
  is Stanford's, and Stanford holds the upstreams, not Menegolo's contribution.
* Its **upstreams are permissive.** The Stanford-VA model (Holzbaur 2005) at
  `simtk.org/projects/up-ext-model` (2008-07-25) is under BSD-3 with a citation requirement
  and **no non-commercial clause**: *"The Stanford-VA upper limb model has been open sourced
  under the BSD 3-Clause License below … Copyright (c) 2005, Stanford University and VA
  Palo Alto Health Care System."* The catalog lists the same model as CC BY 3.0 (Jul-08),
  and so does Arm26's own `<credits>`. The Holzbaur 2005 paper states no terms (*"The
  computer model is available to researchers at http://nmbl.stanford.edu."*).
* It uses the deprecated Schutte muscle law, and §5 already rules it out as a citable
  source.

Relevant to the owner's choice, stated only as fact: Stanford-VA (Holzbaur 2005) is
MoBL-ARMS's catalogued forerunner. It is the only upper-limb musculature among these
sources released without a non-commercial restriction. The catalog also says of it: *"Due
to no inertial properties for the bodies, this model is inappropriate for dynamics
analysis."*

### 10.4 Records changed

* `data/sources/shoulder_complement.json` gains a `license_evidence` block with the
  verbatim notice, the issuer, what each statement covers, and the unresolved points. The
  original four-line `license` block is left unchanged. **`scripts/materialize_shoulder_complement.py`
  rewrites that file from a literal and will drop the new block.** Fix the script before
  re-running it.
* `data/sources/sensorimotor/upperbody_sources.json`'s `unacquired_richer_source` no longer
  says "no model bytes acquired". It now points to the 2026-09-06 mirror acquisition and
  keeps the original status as history (the correction §5 asked for).
* The fetched pages behind 10.1–10.3 were read into scratch and are **not** retained in the
  repo. Their URLs and dates are above, and the licence strings' sha256 values are in the
  JSON.
* Nothing here changes a gate. `registration.gates[0]` ("resolve donor provenance and reuse
  scope") stays open: the provenance facts are established, but the reuse scope is the
  owner's call.

---

## 11. The girdle, built — 18 September 2026

`data/models/shoulder_girdle_v1`. Full record and every number:
**`docs/SHOULDER_GIRDLE.md`**. This section is the part that belongs to *this*
file — which donor, why, and what §8's five steps came out as.

**The donor is Seth 2019, not MoBL-ARMS, and §10 is why.** §10.2 established that
the SimTK notice governing MoBL-ARMS 4.1 permits non-commercial use only, and
that the reuse scope is the owner's call. That call was not made, so the donor was
not used. The Thoracoscapular Shoulder Model (Seth, Dong, Matias, Delp 2019,
Front. Neurorobot. 13:90) is CC BY 4.0 from the authors on SimTK and Apache-2.0
as a file of the opensim-core repository these bytes were read from, and **both
statements permit commercial use with attribution** — they converge where
MoBL-ARMS's diverge, which is the whole difference. `registration.gates[0]` on
`data/sources/shoulder_complement.json` stays open, untouched: nothing here needs
it resolved.

**Provenance established** (`docs/SHOULDER_GIRDLE.md` §1): sha256, the single
opensim-core commit `62205cd879d0` (PR #2971, 2021-03-12) that added the file and
the fact that nothing has touched it since, both licence statements quoted from
primary sources read live, and the caution that the Apache reading is an
inference from the contributor agreement and not a statement attached to the
artefact — the same shape of reasoning that produced the wrong Apache reading of
MoBL-ARMS.

**Provenance NOT established, and it is the same wall:** byte-identity with the
SimTK CC BY release. `download_confirm.php` returns 216 bytes redirecting to
`/account/login.php`. The copy is corroborated as the same model — its 33 muscle
lines are exactly the sixteen groups the paper names, and both the file and the
package are pre-publication artefacts — but not verified as the same bytes.

**The seven `.vtp` meshes are absent and cost nothing mechanically.** Measured:
the engine loads the donor with seven `Couldn't find file` warnings and reports
33 actuators and a full assembly. They are display surfaces; the wrap surfaces
are analytic and in the file. This closes the "geometry not acquired" line in §5
for this donor: nothing mechanical is waiting on it.

### §8's five steps, as they came out

1. **Girdle bodies** — added, with the donor's own topology: `sternoclavicular`
   (CustomJoint, 2 coordinates), `scapulothoracic` (**ScapulothoracicJoint**, a
   Simbody ellipsoid mobilizer, 4), and the acromioclavicular `PointConstraint`
   that closes the loop, per side. `acromial_{r,l}` is re-parented onto the
   scapula, so §8.1's second option — welding the girdle to `torso` and losing
   every girdle-only path — was not needed. 22 → 25 → **29 bodies**, 48 → **60
   coordinates**.
2. **Mass partition, once** — done, and it composed exactly as
   `docs/ARTICULATED_SPINE.md` predicted: the girdle came out of the 17.938 kg
   residual core through the same `partition_body`, reconstruction residual
   **0.0 kg** and **7.1e-16 kg·m²**. The donor's masses (clavicle 0.1898 kg,
   scapula 0.5016 kg) are transferred unscaled and declared so. **The donor's own
   clavicle, scapula and radius inertia TENSORS are not physically realisable** —
   negative second-moment eigenvalues — and were replaced by convex-envelope
   priors of this body's own anatomy, never clipped.
3. **Kinematic reconciliation before any force** — done, and the map is two
   measured similarity factors plus a translation: `s_lat` = 1.31535 from the
   glenohumeral half-width for the girdle, `s_long` = 1.10236 from the
   glenohumeral-to-elbow length for humerus attachments, and a translation that
   puts the donor's glenohumeral centre exactly on this plant's. **The arm does
   not move**: measured against `articulated_spine_v1`, eight bodies × four
   stations, worst difference below 1e-9 m.
4. **Reconcile with what is already there** — `arm26_TRIlong`, `arm26_BIClong`
   and `arm26_BICshort` are kept and the donor's `TRIlong`, `BIC_long` and
   `BIC_brevis` are **not** transferred, so nothing is duplicated or silently
   replaced and the Millard and Thelen laws stay distinct. The arm26 origins are
   still on `torso` where their registration put them; moving them to the scapula
   is a re-registration that was not done.
5. **Gates** — `scripts/verify_shoulder_girdle.py`, 16 tests, all pass, including
   the two-sided moment-arm check (a scapulohumeral muscle must read exactly zero
   about every girdle coordinate, a thoracoscapular one must read nonzero about
   its own side and zero about the other) with this repository's standing
   controls: the query is a function, and `soleus_r` about `ankle_angle_r` still
   reads **−0.0497080 m**. `native_verified` stays **false**: there is no
   linearization, no stance acceptance and no controller.

   **§8.5's warning about the left side was right, and it cost the most.** The
   sternoclavicular mirror worked first try with this plant's own `−MIRROR·a`
   axis rule. The scapulothoracic one did not: its coordinates come from the
   ellipsoid mobilizer, not from declared axes, and the naive mirror put the left
   scapula **162.5 mm** off. All 48 candidates were built and loaded and scored
   on exact mirror symmetry; two are exact and the next is 6.4 mm out. **On the
   LEFT, `scapula_elevation_l` and `scapula_upward_rot_l` run opposite to their
   right-side namesakes** — written down here as well as there, because it is the
   gait2392-knee sign trap and it is now in this plant.

### The answer to §3's open question

§3 ends *"nothing here bounds"* pressing the trunk up. It does now, and the
girdle is what made it boundable rather than what made it possible: in
`articulated_spine_v1` the humerus is jointed straight to `torso`, so at the top
of a press-up the load reaches the trunk through joint reactions with **no muscle
in the path at all**. Measured at the rest pose, against a body weight of
761.1 N:

| | two arms | limited by |
|---|---:|---|
| arms straight, only the scapula to hold | **3265 N = 4.29 body weights** | girdle, `scapula_upward_rot` |
| worst-case lever (force ⟂ each segment, full length) | **345 N = 0.45 body weights** | **elbow**, 43.7 N·m |

**The girdle is not the limiting element. The elbow is**, and it is unchanged by
this build — three arm26 triceps heads, because the donor's `TRIlong` was
excluded as a duplicate. The shoulder itself went from 18.7/31.8, 37.0/10.7 and
5.6/2.9 N·m about flexion, adduction and rotation to **116.5/169.9, 151.6/151.5
and 81.3/90.2**, four to thirty-one times more and bidirectionally balanced for
the first time.

**This is capacity and not behaviour.** `Σ Fmax·|r|` at one pose is a ceiling,
the per-axis sum is not a tension-feasible torque cone — the caution §2 already
attaches to this exact quantity — and nothing here integrates. The body has not
been shown to press itself up.

---

## 12. The wrist is blocked on a licence, and the thoracic joint is not blocked any more — 18 September 2026

Full record and every number: **`docs/WRIST_AND_THORAX.md`**. This section is
the part that belongs to *this* file: what §5's census turned out to be missing,
and what §10's method found when it was pointed at the wrist.

### The wrist donor, and the verdict is STOP

§5's table lists `WristModel` (Gonzalez 1997) as *"25 Schutte1993_Deprecated;
includes `ECU_pre-`/`post-surgery` alternatives; deprecated law"* and says
nothing about its terms. **Its terms are the reason it cannot be used, and they
are stricter than MoBL-ARMS's.**

A whole-repository census — 1,641 `.osim` files, **414 unique by content** —
finds **four** files with a muscle whose path touches a carpal or metacarpal
body. Three are the same Gonzalez model (the `opensim-models` copy at schema
40000 and two `opensim-core` copies at 10905 whose credits are the placeholder
`Model authors names..`). The fourth is `PushUpToesOnGroundWithMuscles.osim`,
which §5 already ruled out and which could not do the job anyway: its
radiocarpal is a one-DOF `PinJoint` with **no deviation coordinate at all**, and
both wrist coordinates ship locked at −1.57079632 rad. Twenty-five further files
declare `wrist_flex_*`/`wrist_dev_*` and, checked by parsing every muscle block
rather than by reading coordinate names, **not one muscle path touches a hand or
carpal body** in any of them.

The SimTK project `wrist-model` (group 325) is run by the three people the model
credits. Its licence field reads **"Custom Use Agreement"** and the agreement
says, verbatim (sha256 of the string `51741011…3aa92`, read live 18 September):

> 4. **You may not copy or distribute this model.** If others are interested in
> using the model, please direct them to this website.
>
> 5. **This model may be used only for non-commercial, academic work. It may not
> be used in any commercial activity.** You may not sell the model or results or
> images generated with the model.

§11 records that the girdle donor was chosen because its two candidate
statements **converge** on commercial use with attribution. Here there is one
statement, issued by the rights holders, and it diverges. **Nothing was copied,
nothing was derived, no wrist muscle was installed**, and `data/raw/` is
gitignored — `git ls-files data/raw` returns zero tracked files — so this
repository does not redistribute it either.

**The Apache copy is not a second opinion.** Taking the byte-equivalent
`opensim-core` copy because that repository has a real Apache-2.0 `LICENSE.txt`
is **precisely the error §10.1(4) recorded for MoBL-ARMS**, and `opensim-core`'s
own `NOTICE` forbids the inference: *"If you use plugins, models, or other
components contributed by your fellow researchers, you must acknowledge their
work as described in the license that accompanies each of these files."* No
licence accompanies these files. The `opensim-models` tree has **no `LICENSE`,
`NOTICE` or `COPYING` anywhere in it**, and yet 33 of its `.osim` files carry an
explicit in-file grant — Arm26's CC BY 3.0 among them. `wrist.osim` carries
none. In a tree where granted models say so, silence is not an inheritance.

### §10.3's alternate, re-read, and it is the right donor

§10.3 recorded the Stanford-VA upper limb model (Holzbaur 2005) as the only
upper-limb musculature among these sources released without a non-commercial
restriction. Confirmed live at `simtk.org/frs/?group_id=324` (string sha256
`037e1fe3…caa92`): BSD-3, copyright Stanford University and VA Palo Alto Health
Care System, citation of Holzbaur 2005 required, **no non-commercial clause**.

**It is behind the same login wall.** `download_confirm.php` returns **197 bytes
of HTML** redirecting to `/account/login.php?triggered=1` — the 216-byte wall
`docs/SHOULDER_GIRDLE.md` §1.4 hit, and the one §10.1 hit for MoBL-ARMS. So §6's
*"(c) — does not apply. No acquisition is required"* is **no longer true of the
wrist**: this is a genuine acquisition gap, and it is one an account can close.

### And the registration would be cheap, which §8-style scoping missed

`docs/ARTICULATED_SPINE.md` scoped the wrist against MoBL-ARMS and concluded *"a
station copy is not possible"* because the frames differ. Against a
Rajagopal/Holzbaur-lineage donor they do not. This plant's arm bodies are
Rajagopal's own bones through `ScaleTool`, which records its factors:
`RajagopalLaiUhlrich2023.osim`'s `radius_hand_r` offset times this model's own
`radius_r` `<Mesh><scale_factors>` reproduces this plant's offset **to exactly
zero on all three components**. An attachment in a generic arm frame maps in by
componentwise multiplication — the rule `build_corrected_foot.py` check C3
already audits over 19 joint offsets. No fit, no landmark estimate, no residual.
The corroborating detail: this plant's wrist ranges and the Gonzalez model's
`flexion`/`deviation` ranges agree to six decimals. One lineage, two files.

**The licence is the blocker. Not the registration, and not the frames.**

### The thoracic joint: `data/models/thoracic_drive_v1`

§4's table says the trunk has *"lumbar extension / bending / rotation"* and
nothing about the thoracic joint, because when it was written nothing crossed
it. `docs/ARTICULATED_SPINE.md` blocked it on the six gait2392 trunk muscles,
whose torso insertion sits 32 mm above a joint centre known to ~62 mm.

That verdict was re-asked against the **sixty** Seth 2019 muscles §11 installed,
32 of which take a `torso` attachment, under a rule committed before it ran
(`1643a94`): an attachment moves to `thorax` only if its **nearest anatomical
surface in this body's own meshes** is one of the 48 structures the composition
plan debited (ribs, cartilages, sternum, intercostals, diaphragm — **no
vertebra**), and only if it sits more than `sqrt(2) × 62 mm` above the joint
centre.

**Ten muscles, five per side, passed**: pectoralis major thoracic I and M,
pectoralis minor, serratus anterior I and M. The refusals are what shows the
instrument works — trapezius, rhomboid and levator scapulae land **0.4–5.7 mm
from a vertebra or an intervertebral disc**, and latissimus dorsi comes back
MIXED with one point 6.4 mm *below* the joint centre. Two corroborations came
for free: `gait2392_ercspn` reads **32.2 mm**, reproducing the number
`ARTICULATED_SPINE.md` quoted from a different computation, and
`arm26_BIClong`/`TRIlong` come back nearest the **scapula**, which is §11's own
open re-registration.

The build moves those ten attachments and the two rib-cage wrap ellipsoids from
`torso` to `thorax` — a 128-line diff in 683 kB, no mass moved, no body or
coordinate added, no parameter touched. **10 gates pass, 1 is recorded FAILED.**
The known answer is that at the rest pose it is the same mechanical object as
its base: all 158 path lengths agree to better than 1e-12 m, every non-thoracic
moment arm to better than 1e-12 m, `soleus_r` still −0.0497080 m, the lumbar arms
still 42.69 / −52.82 / −62.79 mm. The control that can fail is re-run in the same
session: the base plant still reads **0.00e+00** about every thoracic coordinate.

**The result is not a controllable trunk, and the shape is the finding.**

| coordinate | crossing | + bound | − bound |
|---|---:|---:|---:|
| `thoracic_extension` | 10 | **903.7 N·m** | **0.0 N·m** |
| `thoracic_bending` | 10 | 370.2 | 368.8 |
| `thoracic_rotation` | 10 | 347.2 | 343.1 |

All ten pull the same way about extension and there is **no antagonist**,
because the antagonist is erector spinae and its origin is on T12, inside the
torso. Bending and rotation are balanced only by left–right mirroring, not by
opposition per side. So the joint went from topologically unmuscled to
one-directionally muscled at a 900 N·m ceiling, against a 30 N·m/rad stop at
±15°, which any search will drive to the bound. `default_enabled` stays false.
The 112–229 mm moment arms are a property of the single lumped thoracolumbar
joint, not a measurement of a human chest. Same `Σ Fmax·|r|` caution as §2.

### The press-up, and this is a correction to how §11's answer reads

Neither piece of work changes `docs/SHOULDER_GIRDLE.md` §7's numbers — 3,265 N
arms-straight, 345 N at the worst lever, **elbow-limited at 43.7 N·m** — and the
thoracic reassignment provably changes no arm about `elbow_flex`, `arm_*` or any
`scapula_*` coordinate, measured to 1e-12 m.

But the wrist changes what that sentence is worth. §11's whole argument was that
**with a girdle the question is bounded, because every newton reaches the trunk
through muscle.** The wrist is now the one link in that chain with no muscle in
it at all, and — measured from the model file — no elastic resistance either:
`wrist_flex_{r,l}` and `wrist_dev_{r,l}` carry only `-1.0*qdot` viscous damping
and a `CoordinateLimitForce` that by construction acts *outside* the declared
range. Inside the range the wrist is a free hinge. **So 43.7 N·m is the smallest
number among the links that have muscle; it is not the smallest number in the
chain.** Still capacity, still not behaviour.
