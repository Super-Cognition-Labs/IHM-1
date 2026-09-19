# the body has a shoulder girdle, and the elbow is now the weak link

`data/models/shoulder_girdle_v1` is a VARIANT of `data/models/articulated_spine_v1`
that adds the thing `docs/WORKBENCH_AUTHENTICITY.md` 0.2 and
`docs/UPPER_BODY_ACTUATION.md` §4 both name as the largest remaining actuation
gap: a scapula, a clavicle, the joints that carry them, and the muscles that move
them. Built and measured 18 September 2026.

`articulated_spine_v1` is untouched, and so is `engineering_stance_v1` under it.
This variant carries **no `linearization.npz`** and no stance acceptance, for the
reason that file already gives: the identified stance plant was solved on 22
bodies and 33 coordinates and does not transfer to 29 and 60.

    measure   .venv/bin/python -m scripts.measure_shoulder_donor_frames
    build     .venv/bin/python -m scripts.build_shoulder_girdle
    verify    OPENBLAS_NUM_THREADS=1 nice -n 10 \
                  .venv/bin/python -m unittest scripts.verify_shoulder_girdle -v
    capacity  OPENBLAS_NUM_THREADS=1 nice -n 10 \
                  .venv/bin/python -m scripts.measure_shoulder_girdle_capacity
    select    NativeMechanicalStream(..., augmented_registration=
                  'data/models/shoulder_girdle_v1/registration.json')

## What is there now

| | `engineering_stance_v1` | `articulated_spine_v1` | `shoulder_girdle_v1` |
|---|---|---|---|
| bodies | 22 | 25 | **29** |
| coordinates | 33 | 48 | **60** |
| muscles | 98 | 98 (148 muscled) | **158** |
| constraints closing a loop | 0 | 0 | **2** |
| total mass | 85.269848541731 kg unscaled | identical to 1e-9 | **identical to 1e-9** |

    torso ──sternoclavicular(2)── clavicle ─┐
    torso ──scapulothoracic(4)─── scapula ──┴─ acromioclavicular PointConstraint
    scapula ──acromial(3)── humerus ── elbow ── ...

Twelve new rotational coordinates, every one with a declared `<range>` and a
`CoordinateLimitForce` **in the model file** at exactly that range. Sixty new
muscles, thirty per side.

---

## 1. Provenance: what the donor is, and what could not be established

### 1.1 What is on disk

`data/raw/mechanics/opensim-core/OpenSim/Tests/shared/ThoracoscapularShoulderModel.osim`,
318,742 bytes, sha256
`892bd333a1e1e6a2b77902b54218f12727eefd2306ca82bbda81621c926271fe`,
`OpenSimDocument Version="40000"`. It holds:

| | |
|---|---|
| bodies | 7: thorax, clavicle, scapula, humerus, ulna, radius, hand |
| joints | `ground_thorax` (CustomJoint, 6), `sternoclavicular` (CustomJoint, 2), `scapulothoracic` (**ScapulothoracicJoint**, 4), `GlenoHumeral` (CustomJoint, 3), `elbow`, `radioulnar`, `rc` (weld) |
| coordinates | 17 |
| constraints | one `PointConstraint`, the acromioclavicular joint |
| muscles | 33 `Millard2012EquilibriumMuscle` |
| path points | 73 fixed, 5 moving; **no conditional points** |
| wrap objects | 10 (2 thorax, 4 scapula, 3 humerus, 1 ulna); 9 are referenced by a muscle |
| markers | 21 |
| geometry | seven `.vtp` references, **none of the files present** |

Its `<credits>` read *"Ajay Seth, Meilin Dong, Ricardo Matias, Scott Delp.
Parameters from van der Helm and Klein-Breteler"* and its `<publications>`
*"Frontiers in Neurorobotics: in In Press DSEM: van der Helm 1994 Klein-Breteler
et al. 1996."* **The file contains no licence, copyright or terms text anywhere.**

### 1.2 Where those bytes came from, established

The file sits in a git checkout of `opensim-org/opensim-core` at
`86b30588374650fbaf012a345a836a64f6855522` — the same revision the native engine
in `data/runtime/opensim` is built from (`docs/NATIVE_OPENSIM.md`) — with a clean
working tree for that path. The checkout is shallow, so its history was read from
the GitHub API instead: the commits endpoint for that path returns **exactly one
commit**, `62205cd879d0`, 2021-03-12T20:48:58Z, author `carmichaelong`, from PR
**#2971** *"Fix IDTool and IDSolver to account for extra q's in state"*, whose
body says *"new ID test with ThorascapularShoulderModel"*. Nothing has touched it
since, so the bytes at `HEAD` are the bytes that commit added. This independently
confirms what `docs/UPPER_BODY_ACTUATION.md` §10.3 recorded.

### 1.3 Licence: two statements, and they agree

**The authors' own, on SimTK.** Project `thoracoscapular`, group 1708, read live
on 18 September 2026 from `https://simtk.org/frs/?group_id=1708`:

> **FrontiersPubMaterials**
> Copyright (c) 2019, Stanford University and the authors.
>
> authors: Ajay Seth, Meilin Dong, Ricardo Matias, Scott Delp
>
> This work is available under the Creative Commons Attribution 4.0
> International Public License, summarized below. For full legal text, please see
> http://creativecommons.org/licenses/by/4.0/legalcode
>
> … Share — copy and redistribute the material in any medium or format. Adapt —
> remix, transform, and build upon the material for any purpose, even
> commercially. … Attribution — You must give appropriate credit, provide a link
> to the license, and indicate if changes were made. … No additional
> restrictions.

The package that licence is attached to is `ThoracoscapularShoulderPaperMaterials.zip`,
196 MB, 2019-06-28, and the page's own metadata labels it
`"version": "Pre-publication Release"` with the citation ending *"Frontiers in
Neurorobotics, in review."*

The published paper agrees. Seth, Dong, Matias and Delp, *Muscle Contributions to
Upper-Extremity Movement and Work From a Musculoskeletal Model of the Human
Shoulder*, Front. Neurorobot. **13**:90 (2019), Data Availability Statement:
*"The model and simulation environment (OpenSim) are freely available,
deployable, and modifiable for any research or commercial use without
restrictions from SimTK.org at https://simtk.org/projects/thoracoscapular"*. The
article itself is CC BY.

**The repository's, on the bytes actually held.** `opensim-core/LICENSE.txt` is
Apache-2.0, and `CONTRIBUTING.md` under *Contributor License Agreement* says
*"by contributing you are agreeing to the following terms: … You grant the
OpenSim project, developers, and users a nonexclusive, irrevocable license to use
your submission and any necessary intellectual property, under terms of the
Apache 2.0 license. … You are capable of granting these rights for the
contribution."* The `ScapulothoracicJoint` class this model needs carries its own
explicit header: *"Copyright (c) 2005-2020 Stanford University, TU Delft, and the
Authors. Author(s): Ajay Seth. Licensed under the Apache License, Version 2.0"*.

**The one caution on the Apache reading, and it is the MoBL-ARMS lesson.**
`opensim-core/NOTICE` ends: *"If you use plugins, models, or other components
contributed by your fellow researchers, you must acknowledge their work as
described in the license that accompanies each of these files."* **No licence
accompanies this file.** So the Apache grant here is an inference from the CLA
and not a statement attached to the artefact, exactly the shape of the reasoning
that produced the wrong "Apache-2.0" reading of MoBL-ARMS. The difference, and it
is the whole difference, is that **both candidate statements converge**: the
authors' CC BY 4.0 and the repository's Apache-2.0 each permit use, modification
and redistribution, commercially, with attribution. For MoBL-ARMS the candidate
statements *diverge* on precisely the commercial question, which is why that
donor is not usable and this one is.

**Obligations this build takes on.** Cite Seth et al. 2019 in anything describing
work that uses it; under CC BY also provide a link to the licence and state that
changes were made. `data/models/shoulder_girdle_v1/registration.json` carries the
citation, the licence statement and the list of changes, and every catalog row
for the sixty transferred muscles names the donor file and its hash.

**MoBL-ARMS 4.1 was not used, directly or indirectly.** No number in the model,
the builder, the verifier or this document comes from
`data/research/shoulder_complement/MOBL_ARMS_41.osim` or `shoulder_forces.xml`.

### 1.4 What could NOT be established, and was checked rather than assumed

* **Byte-identity with the SimTK release. NO.** `https://simtk.org/frs/download_confirm.php/file/5812/ThoracoscapularShoulderPaperMaterials.zip?group_id=1708`
  returns **216 bytes of HTML** redirecting to `/account/login.php?triggered=1`.
  The package is behind the same login wall as MoBL-ARMS and the bytes cannot be
  fetched. So the on-disk copy is **corroborated as the same model but not
  verified as the same bytes**, and this document does not claim otherwise.

  What corroborates it: its 33 muscle lines group into exactly the sixteen muscle
  groups the paper names (trapezius, serratus anterior, rhomboideus, levator
  scapulae, coracobrachialis, deltoideus, latissimus dorsi, pectoralis major,
  teres major, infraspinatus, pectoralis minor, teres minor, subscapularis,
  supraspinatus, triceps long, biceps); its credits name the four authors and the
  DSEM parameter source; and the on-disk file and the SimTK package are both
  pre-publication artefacts — *"in In Press"* in the file, *"Pre-publication
  Release"* and *"in review"* on the page.

* **Geometry: the seven `.vtp` files are absent, and it costs nothing
  mechanically.** Measured, not assumed: the engine loads the donor model with
  seven `[warn] Couldn't find file` lines and then reports 33 actuators and a
  full assembly. The meshes are display surfaces; the model's wrap surfaces are
  analytic and are in the file. Acquiring them means the same login-walled
  196 MB package. **The girdle built here attaches no display geometry**, and the
  workbench consequence is in §7.

* **A second copy exists in the same repository and DIFFERS.**
  `OpenSim/Tests/Wrapping/testGeometryPathWrapping_ShoulderWrappingWithPathSprings.osim`
  carries the same six girdle coordinates; its `shoulder_elv` range is
  `-1.5 … 3.14159` where this one declares `-0.78 … 2.53`. Not used here, and
  noted so that "the thoracoscapular model in opensim-core" is never treated as
  one artefact.

* **The DSEM parameter sources' own terms** (van der Helm 1994, Klein-Breteler
  et al. 1996) were not examined. Unchanged from `docs/UPPER_BODY_ACTUATION.md`
  §10.3.

### 1.5 A defect in the donor, found by a check that exists for this

**The donor's clavicle, scapula and radius inertia tensors are not physically
realisable.** For each, the second-moment matrix `tr(I)/2·I − I` has a negative
eigenvalue — `−3.07e-05`, `−9.79e-05` and `−1.40e-04` kg·m² — so the principal
moments violate the triangle inequality, and the radius has a negative principal
moment outright (`−4.54e-06` kg·m²). No rigid mass distribution has those
tensors. `ihm.assembly.cervical_inertia._physical` refuses them, which is how
this surfaced: the build threw instead of shipping.

They are **replaced, never clipped**, which is this repository's standing
decision for exactly this case — every cervical body in
`data/research/cervical_inertia/v2/manifest.json` carries the string *"Replace
donor tensor with geometry-conditioned prior; no eigenvalue clipping and no
Body.cpp fallback"*. The replacement is the convex-envelope prior of **this
body's own** anatomical clavicle and scapula (`body-bp3d-FJ3362`, `FJ3384`,
`FJ3237`, `FJ3279`) at the donor's stated mass, carried into the new body's axes
by the same canonical→torso rotation the cervical partition used. A *rotation*,
not a translation: a central second moment does not depend on where the mesh
sits, so that registration's 62 mm positional residual does not enter. Each side
uses its own side's mesh, so there is no mirror prior in the inertia at all.

| | donor principal moments (kg·m²) | replacement |
|---|---|---|
| clavicle | 1.44e-05, 2.47e-04, 3.23e-04 | 2.14e-05, 3.01e-04, 3.12e-04 |
| scapula | 5.09e-04, 6.88e-04, 1.39e-03 | 3.03e-04, 8.50e-04, 9.87e-04 |

Same order of magnitude, and the replacement is admissible. `verify_shoulder_girdle`
checks both halves: the new tensors are admissible **and** the donor's still are
not, because a test that would pass if the donor had been fine all along is not a
test.

---

## 2. The map between the donor's frames and ours

Two similarity factors and one translation, each from a correspondence both
models declare in their own files, and nothing fitted.

| | value | from |
|---|---|---|
| `s_lat` | **1.3153489519566615** | this plant's glenohumeral half-width (`acromial_r`'s parent offset on `torso`, z = 0.198768021792 m) divided by the donor's glenohumeral centre in its thorax frame (z = 0.15111428909896535 m) |
| `s_long` | **1.102364795064231** | this plant's glenohumeral-to-elbow distance (0.3204855851242984 m) divided by the donor's (0.2907255262135116 m), both from the models' own joint frames |
| `delta` | **(0.058290596211699276, 0.40654277409938316, 0.0)** m | x and y put the donor's glenohumeral centre exactly on this plant's; **z is exactly zero** |

**Why two factors and not one.** The girdle spans the midline to the glenoid, and
`s_lat` is the ratio of exactly that distance — so the clavicle and the scapula
are scaled by the measurement that sets their size. Because `delta_z` is zero,
the donor's midline maps to this plant's midline, and a trapezius or rhomboid
origin on a spinous process stays on the midline. A single factor anchored on the
glenohumeral centre would have pushed every midline attachment **32 mm** off it,
which for moment arms of 20–60 mm is an order-unity error. The humerus is a
different case: it is this plant's own body and is not being replaced, so donor
attachments that land on it are scaled by `s_long`, the ratio of the two models'
humerus lengths, and rotated by the scapula's reference rotation because the
donor's humerus frame at zero glenohumeral coordinates is scapula-aligned while
this plant's is torso-aligned.

**The donor's reference pose is not readable from its XML.** The scapulothoracic
joint is a Simbody ellipsoid mobilizer and the girdle loop is closed by a
`PointConstraint`, so the scapula's pose is the solution of an assembly, not a
stored transform. Reimplementing the mobilizer in Python would be a second
implementation of the thing being measured, so it is measured with the engine
instead (`scripts/measure_shoulder_donor_frames.py`): sixteen two-point
`PathActuator` rulers per body, from four non-coplanar thorax stations to four
stations on the target body, trilaterated. A `PathActuator` with two points is a
ruler, and path length is the one quantity the sampler reports.

Its known answers, all passed:

| | |
|---|---|
| recovered rotations orthonormal | max abs error **1.8e-15**, determinant +1 |
| the sixteen distances reproduced by the stored rigid transform | **1.7e-16 m** |
| the clavicle origin equals the sternoclavicular offset the XML declares | exactly |
| the **independently measured** humerus origin equals the glenohumeral centre computed from the scapula frame and the XML | **3.5e-17 m** |

The last one is the one that matters: two routes to the same point, no shared
arithmetic.

**The arm does not move.** `acromial_{r,l}` is re-parented from `torso` onto
`scapula_{r,l}` through an offset frame at the donor's glenohumeral centre in the
scapula, carrying the **inverse** of the scapula's reference rotation, so at zero
arm coordinates the humerus frame is identical to the base plant's. That is
arithmetic; the verifier measures it anyway, four stations per arm body against
`articulated_spine_v1` in the same protocol, and the worst difference over eight
bodies is **below 1e-9 m**.

---

## 3. The one-time torso partition

`docs/ARTICULATED_SPINE.md` is explicit that its own partition is the first half
of one job: its thoracic share is rib cage only — 2.298 kg of ribs, cartilages,
sternum, intercostals and diaphragm — with no mass in it for a scapula, so the
girdle has to come out of the **17.938 kg residual core still inside `torso`**,
and it composes because both debits go through the same instrument.

They do. `ihm.assembly.cervical_inertia.partition_body`, XML (unscaled) units:

| body | before | after |
|---|---:|---:|
| `torso` | 19.7075366916046 kg | **18.188364076006206 kg** |
| `clavicle_r`, `clavicle_l` | — | 0.20849694475759684 kg each |
| `scapula_r`, `scapula_l` | — | 0.5510893630416 kg each |

Reconstruction residual **0.0 kg**, first moment 6.8e-17 kg·m, inertia Frobenius
**7.1e-16 kg·m²**. Total model mass conserved to 1e-9 with no body but `torso`
changed, and the verifier reassembles the five parts from the **written model
file** rather than replaying the builder's arithmetic.

**The masses are the donor's and are NOT scaled to this body**, because the donor
model states no body mass to scale from. That is the stance this repository
already takes for donor muscle forces (`catalog.json`, *"Unscaled retained donor
parameters"*), and it is declared rather than hidden. The mass centres are the
donor's, scaled by `s_lat`. The inertia tensors are not the donor's at all (§1.5).

---

## 4. The left side had to be measured, and assembly was no test

`docs/UPPER_BODY_ACTUATION.md` §8.5 warns that a reflection prior *"needs explicit
validation of vectors, joint axes and wrap quadrants — not a name suffix"*. The
donor has one arm. Here is what the validation found.

The **sternoclavicular** joint is a `CustomJoint` with declared axes, so the rule
this plant already uses works: each left rotation axis is `−MIRROR·a`, which is
exactly the relation this model's own `acromial_l` axes have to `acromial_r`'s.
Measured result: the left clavicle's rest pose is the z-mirror of the right's to
**0.000000 m**, first try.

The **scapulothoracic** joint is not. Its coordinates are not user-supplied axes;
they come from a Simbody ellipsoid mobilizer whose conventions are fixed inside
`ScapulothoracicJoint.cpp`, so mirroring the parent frame does not tell you what
the mirrored coordinates mean. The naive choice — conjugate the rotation by
`diag(1,1,−1)`, keep the donor's defaults — put the left scapula **162.5 mm**
away from the mirror of the right and left the left acromioclavicular constraint
unsatisfiable by **83 mm**, while the right side assembled to 1e-10 with nothing
moved.

So all 48 candidates were built and loaded: three ways to restore a right-handed
mirrored frame (right-multiply by `diag(−1,1,1)`, `diag(1,−1,1)` or
`diag(1,1,−1)`) times sixteen coordinate sign patterns, each scored on the only
thing that is not a matter of taste — is the left girdle's rest pose the exact
z-mirror of the right's?

| candidate | origin residual | rotation residual |
|---|---:|---:|
| `diag(1,−1,1)` with `(+1, −1, −1, +1)` — **used** | 1.4e-16 m | 1.3e-15 |
| `diag(−1,1,1)` with `(−1, +1, −1, −1)` | identical — the same configuration named twice | |
| next best of the 48 | 6.4e-03 m | 3.9e-02 |

Two exact answers and a 6.4 mm gap to the third. This is a discrete fact about
the mobilizer, not a tuned threshold.

**Assembly success was useless as a test, and that is worth keeping.** Sixteen of
the 48 assembled, wrong ones included, because the acromioclavicular constraint
is three equations on a six-coordinate girdle and the solver simply moves the
coordinates off their declared defaults until it closes. A control that passes
for a reason unrelated to what it tests.

> **Sign convention, written down before it can catch anyone.** On the LEFT,
> `scapula_elevation_l` and `scapula_upward_rot_l` run **opposite** to their
> right-side namesakes. `scapula_abduction_l`, `scapula_winging_l`,
> `clav_prot_l` and `clav_elev_l` run the same way. This is the same shape as the
> gait2392 knee in `docs/NATIVE_JOINT_LIMITS.md`: two artefacts that declare the
> same quantity with opposite signs, and a gate that only ever looks at one of
> them cannot see it.

---

## 5. Ranges and stops

Twelve new coordinates, each with a declared `<range>` and a
`CoordinateLimitForce` **in the model file** at exactly that range, 30 N·m/rad,
transition 0.35 rad — `scripts/crawl.py`'s constants and the only stiffness this
repository has swept (`docs/NATIVE_JOINT_LIMITS.md`). The same damping-unit
disagreement `docs/ARTICULATED_SPINE.md` records applies here too and is written
the same way, so the two variants agree with each other.

**The ranges have no measured provenance and are labelled so, in the builder, in
the registration and here.** Each is the donor's own default value ± 0.349066 rad
(20°). The donor declares ±π for five of its six girdle coordinates and ±π/2 for
the sixth — the same placeholder as this plant's ±10 rad shoulder coordinates,
not a range — and the motion file shipped beside it in opensim-core is synthetic
(every column is 0.01, 0.02, 0.03 … held constant over six rows). No other donor
on disk declares a scapular or clavicular range of motion. Narrow is the
conservative side for a stop whose only job is to keep a search out of absurd
configurations, which is the position `docs/ARTICULATED_SPINE.md` takes for its
thoracic range.

*What would replace it:* recorded scapular kinematics, or the donor's own
experimental data — which is inside the login-walled 196 MB package (§1.4).

**The rest pose is the donor's own girdle configuration, not zero**, because the
acromioclavicular constraint is satisfied exactly there and at no value this
build could invent. All twelve rest strictly inside their declared range.

---

## 6. What the checks print

`scripts/verify_shoulder_girdle.py`, **16 tests, 6.6 s, all pass.** Two are the
ones to read.

| check | result |
|---|---|
| total mass conserved against `articulated_spine_v1` | 0 to 1e-9; no body but `torso` changed |
| the five parts recombine to the original one torso | mass **0.0 kg**, COM < 1e-12 m, inertia Frobenius **7.1e-16 kg·m²** |
| every new body's inertia is physically admissible | pass — **and the donor's still are not**, re-checked in the same run |
| every new coordinate declares a range and rests strictly inside it | pass, 12/12 |
| every new coordinate carries a stop AT that range, and no old one gained one | pass, 12/12 |
| the donor declares no usable range | pass — all six are ±π or ±π/2 |
| the model loads in the engine | 29 bodies, 158 muscles, 60 coordinates, 77.6122029 kg |
| **the arm did not move when it was re-parented** | **< 1e-9 m** over 8 bodies × 4 stations, against the base plant |
| **the left girdle is the exact z-mirror of the right** | **< 1e-9 m**, position and orientation |
| a single girdle coordinate cannot be set on its own | refused, and that is the closed loop |
| each new joint moves only its distal entities | driven side 3.10e-3 … 5.54e-2 m; everything else ≤ 3.44e-6 m |
| the builder is a function (run twice → identical bytes) | pass |
| moment arms, with both controls | `soleus_r` about `ankle_angle_r` **−0.0497080 m**, repeat call bit-equal |

**The closed loop refuses a single coordinate, and that is correct.** Each girdle
is `torso → clavicle` (2 coordinates) plus `torso → scapula` (4) closed by a
3-equation `PointConstraint`, so only three are free. The engine answers
`assembly changed requested initial pose` rather than silently returning a
different pose. The distal-only test therefore goes through OpenSim's own
assembler instead, measuring in the **torso frame** by trilateration.

**The distal-only bar had to be measured, and the first one was below the noise.**
Displacements split into two populations 900× apart: driven side 3.10e-3 to
5.54e-2 m, everything else 0 to 3.44e-6 m — and all but one of those is below
3.2e-11 m. The exception is the **contralateral hand**, 2.47e-6 and 3.44e-6 m,
while the contralateral humerus moves 3e-16: the assembler re-solves the whole
body (two knee couplers plus two acromioclavicular constraints) and nudges the far
wrist. A bar at 1e-6, which is what the test was first written with, sits below
that and fails on it — it did, and it read as a topology error. Same lesson
`docs/ARTICULATED_SPINE.md` already records for the same reason.

### Moment arms, measured

Sixty muscles, and the two-sided check is what makes the numbers mean anything: a
**scapulohumeral** muscle must read exactly zero about every girdle coordinate,
because it spans scapula to humerus and crosses no girdle joint, and a
**thoracoscapular** muscle must read nonzero about its own side and zero about
the other. A build that put the muscles on the wrong bodies breaks one; a build
that put them all on one body breaks both.

| muscle | about | m |
|---|---|---:|
| `seth_SerratusAnterior_M_r` | `scapula_upward_rot_r` | **+0.0643** |
| `seth_DeltoideusScapula_M_r` | `scapula_upward_rot_r` | **0** (spans scapula→humerus) |
| `seth_DeltoideusScapula_M_r` | `arm_flex_r` | +0.0103 |
| `soleus_r` (control) | `ankle_angle_r` | **−0.0497080** |

---

## 7. Can the body press its own trunk up?

`docs/UPPER_BODY_ACTUATION.md` §3 ends: *"Pressing the trunk up off the floor is
a demand of a different order and **nothing here bounds it**."* The girdle is what
makes the question answerable, and not because it adds strength.

**Without a girdle the question has no answer.** In `articulated_spine_v1` the
humerus is jointed straight to `torso`, so at the top of a press-up — trunk
horizontal, arms straight, hands under the shoulders — the ground reaction runs
hand → radius → ulna → humerus → glenohumeral joint → torso entirely through
**joint reactions**. No muscle is in the load path and the skeleton carries any
force whatever. The quantity is not unbounded because the plant is strong; it is
undefined because nothing carries it.

**With the girdle it is bounded**, because the scapula has no bony attachment to
the axial skeleton except the clavicle: every newton the hand pushes with reaches
the trunk through the thoracoscapular muscles.

Measured (`scripts/measure_shoulder_girdle_capacity.py`), at each model's own
rest pose, from `Σ Fmax·|moment arm|` read out of the engine:

| joint | capacity, weaker direction | lever | force, one arm |
|---|---:|---:|---:|
| girdle, `scapula_upward_rot_r` | 175.4 N·m | 0.1075 m | **1633 N** |
| shoulder, `arm_flex_r` | 116.5 N·m | 0.5701 m | 204 N |
| elbow, `elbow_flex_r` extension | **43.7 N·m** | 0.2534 m | **173 N** |

Body weight is **761.1 N** (77.6122029 kg × the model's own gravity). Two
readings bracket the answer:

* **Arms straight, hands under the shoulders** — the position the trunk has to be
  held in. The load line runs through the elbow and the glenohumeral centre, so
  their levers go to zero and no arm muscle is in the path; the only muscular
  demand left is holding the scapula. Two arms: **3265 N = 4.29 body weights.**
  **The girdle is not the limiting element, by a factor of four.**
* **Worst-case lever** — the hand force perpendicular to each segment at the
  segment's full length, the hardest a bent arm can ever be loaded. Two arms:
  **345 N = 0.45 body weights**, and the binding element is the **elbow**.

### The answer

**Yes for holding the top of a press-up, no for pushing up from a bent elbow, and
the girdle is not what stops it.** The thoracoscapular muscles can carry 4.29
body weights at the scapula. The elbow extensors — three arm26 triceps heads,
43.7 N·m, unchanged by this build — are the smallest number in the chain. The
donor's own `TRIlong` was deliberately not transferred, because it would have
duplicated `arm26_TRIlong` (`docs/UPPER_BODY_ACTUATION.md` §8.4), so the elbow
is exactly as strong as it was before the girdle.

**And this is capacity, not demonstrated behaviour.** Nothing here integrates,
nothing here is a controller, and the body has not been shown to press itself up.
`Σ Fmax·|r|` is full activation at optimal fibre length with no force–velocity,
no tendon state and no activation dynamics: a ceiling. A per-axis sum is not a
tension-feasible torque cone — the caution `docs/UPPER_BODY_ACTUATION.md` §2
already attaches to this exact quantity. Moment arms are pose dependent and these
are the rest pose's; no press-up pose was evaluated. The levers are straight-line
distances, not perpendicular distances to a line of action.

### What the girdle did change, measured

`Σ Fmax·r` in each direction, this plant against `articulated_spine_v1`, same
protocol, same run:

| coordinate | before (+ / −) | after (+ / −) |
|---|---:|---:|
| `arm_flex_{r,l}` | 18.7 / 31.8 | **116.5 / 169.9** |
| `arm_add_{r,l}` | 37.0 / 10.7 | **151.6 / 151.5** |
| `arm_rot_{r,l}` | 5.6 / 2.9 | **81.3 / 90.2** |
| `elbow_flex_{r,l}` | 22.9 / 43.7 | 22.9 / 43.7 — unchanged, by design |

The shoulder is four to thirty-one times stronger and, for the first time,
**bidirectionally balanced**: before, adduction had 37.0 N·m against 10.7 in the
other direction, from three biceps and triceps heads doing duty as a shoulder.

---

## 8. What this variant does not do

* **No `linearization.npz`, no stance acceptance, `default_enabled: false`.**
  Nothing in the app selects it; `ihm/assembly/embodied.py` hard-codes
  `engineering_stance_v1` when a controller is chosen. A caller can ask for it
  today through `augmented_registration=`.
* **No display change.** The girdle attaches no geometry, and
  `data/derived/anatomy-segment-binding/` assigns every anatomical entity to one
  of the original 22 segments, so both scapulae and both clavicles still ride
  `torso` in the workbench. Closing that is the same job
  `docs/ARTICULATED_SPINE.md` leaves open for the head and the rib cage, plus the
  seven donor `.vtp` meshes, which are behind the login wall in §1.4.
* **The three latissimus dorsi via-points are frozen.** In the donor they are
  `MovingPathPoint`s driven by its `shoulder_elv`, a coordinate this plant does
  not have — its shoulder is a flexion/adduction/rotation decomposition, not an
  elevation plane. They are frozen at the donor's own value at `shoulder_elv = 0`,
  which is the arm hanging and therefore exact at this model's rest pose. The
  donor moves them by at most **16.2 mm** over 180° of elevation (21.3 mm after
  `s_lat`); the per-axis travel of all five is in `registration.json`.
* **Three donor muscles were not transferred**: `TRIlong`, `BIC_long` and
  `BIC_brevis`, which attach to the donor's ulna or radius and would duplicate
  `arm26_TRIlong`, `arm26_BIClong` and `arm26_BICshort` already in the plant.
  Millard (donor) and Thelen (arm26) laws stay distinct.
* **`arm26_TRIlong`, `arm26_BIClong` and `arm26_BICshort` still originate on
  `torso`**, where their registration put them, although anatomically they
  originate on the scapula. Moving them is a re-registration with its own
  provenance and was not done here. They now cross the girdle joints, so they
  acquire girdle moment arms that a scapular origin would change.
* **The ten girdle-anchored MASI neck muscles are now unblocked but not
  transferred.** `docs/ARTICULATED_SPINE.md` excluded cleidomastoid,
  cleido-occipital, trapezius clavicular and acromial, and levator scapulae for
  want of a clavicle and a scapula. This model has both. Transferring them means
  reconciling the MASI registration's girdle placement with this donor-anchored
  one, and the two disagree (§9) — a separate job.
* **No actuator port on any new coordinate**, and no controller. The girdle is
  muscle-driven or nothing.
* **`equilibrium_excitations.json` is the base plant's and names only its 98
  muscles.** It is carried across unchanged, as `articulated_spine_v1` carries it
  across for its own 148-muscle plant, so a caller that drives "every muscle at
  its equilibrium excitation" will leave all sixty new ones at zero. That is the
  right default for an unaccepted variant and it is not an equilibrium claim.
* **No wrist or forearm muscle.** Unchanged: those are blocked on MoBL-ARMS and
  its licence, not on the girdle (`docs/ARTICULATED_SPINE.md`).

---

## 9. One disagreement worth recording

Anchoring the donor girdle on this plant's own glenohumeral centre puts the
sternoclavicular joint at **y = 0.41566 m** in the torso frame, which is
**41.4 mm ABOVE** this plant's C7/T1 `neck` joint centre at 0.37427 m. Anatomically the
sternoclavicular joint sits below C7/T1, at about the T2–T3 level.

The two placements come from two different registrations of two different donors:
the neck joint from `data/research/cervical_registration/v2/recipe.json` (rms
14.9 mm over eight bone-envelope centroid proxies) and the girdle from the
glenohumeral correspondence used here. The disagreement is the same order as the
62 mm RMS residual the torso registration already declares
(`data/research/cervical_inertia/v2/manifest.json`), so it is not evidence that
either is wrong — but it is evidence that **the upper thorax of this plant is
known to a few centimetres and no better**, and anything that needs the two
girdle and neck registrations to be consistent with each other has to reconcile
them first. That includes the ten neck muscles in §8.
