# the wrist is blocked on a licence, and the thoracic joint is not blocked any more

Two joints were recorded as having no muscle that can drive them
(`docs/ARTICULATED_SPINE.md`, "What is still blocked, and on what"). This file
re-asks both, on 18 September 2026, and they come out opposite ways.

* **The wrist. STILL BLOCKED, and the blocker is a licence, not a registration.**
  Exactly one wrist-muscled donor is on disk in any form. Its authors' own
  licence says *"You may not copy or distribute this model"* and *"This model
  may be used only for non-commercial, academic work."* Nothing was copied,
  nothing was derived from it, and no wrist muscle was installed. §1–§5.
* **The thoracic joint. UNBLOCKED, for ten muscles the girdle build brought in.**
  `data/models/thoracic_drive_v1`: five muscles per side reassigned from `torso`
  to `thorax` under a rule fixed before it ran. 10 gates pass, 1 is recorded
  FAILED. §6–§8.
* **Neither changes the press-up.** The elbow's 43.7 N·m is still the limiting
  element, and the wrist is the reason that number is not even the whole story.
  §9.

      measure   OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
                    .venv/bin/python -m scripts.measure_thoracic_reassignability
      build     .venv/bin/python -m scripts.build_thoracic_drive [--check]
      verify    OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
                    .venv/bin/python -m unittest scripts.verify_thoracic_drive -v
      capacity  OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
                    .venv/bin/python -m scripts.measure_thoracic_drive_capacity
      select    NativeMechanicalStream(..., augmented_registration=
                    'data/models/thoracic_drive_v1/registration.json')

`data/models/engineering_stance_v1`, `articulated_spine_v1` and
`shoulder_girdle_v1` are read and never written.

---

# PART ONE — THE WRIST

## 1. What is on disk, and it is a short list

A census of every `.osim` in the repository: **1,641 files, 414 unique by
content** (the rest are copies under `data/derived/` and `data/runtime/`).
`data/raw` holds 215. Of the 414, **four** have a muscle whose path touches a
carpal or metacarpal body.

| file | bytes | doc | muscle classes | wrist muscles |
|---|---:|---:|---|---:|
| `data/raw/anatomy/opensim-models/source/Models/WristModel/wrist.osim` | 419,036 | 40000 | 25 `Schutte1993Muscle_Deprecated` | 24 |
| `data/raw/mechanics/opensim-core/Bindings/Java/OpenSimJNI/Test/wrist.osim` | 381,258 | 10905 | 26 `Schutte1993Muscle` | 25 |
| `data/raw/mechanics/opensim-core/OpenSim/Simulation/tests/resources/wrist_mass.osim` | 382,251 | 10905 | 26 `Schutte1993Muscle` | 25 |
| `data/raw/mechanics/opensim-core/OpenSim/Simulation/tests/resources/PushUpToesOnGroundWithMuscles.osim` | 989,833 | 30000 | 100 `Schutte1993Muscle_Deprecated` + 54 `Thelen2003Muscle` | 48 |

sha256: `66488a78…c1996`, `3a702ad1…b60b`, `78f07837…e619`, `5b8344e1…95d6`.

**The first three are one model.** `data/raw/anatomy/.../WristModel/wrist.osim`
credits `Gonzalez, R.V., Buchanan, T.S., Delp, S.L.` and publishes *"How muscle
architecture and moment arms affect wrist flexion-extension moments. Journal of
Biomechanics, vol. 30, pp. 705-712, 1997."* 28 bodies, 14 coordinates, the
carpals welded to a two-body flexion/deviation carrier chain, and `flexion`
(±1.22173048 rad) and `deviation` (−0.436332 … +0.610865 rad). The two
`opensim-core` copies are the same model at schema 10905 with one extra muscle
(a plain `ECU` beside the `ECU_pre-surgery`/`ECU_post-surgery` pair) and
**placeholder credits** — `<credits>Model authors names..</credits>`.

> **A census trap worth keeping.** The two `opensim-core` copies use the legacy
> `<Schutte1993Muscle>` tag **without** the `_Deprecated` suffix. A search for
> the three modern class names reports them as having *zero* forces. They have
> twenty-six.

**The fourth cannot do the job even if it were clear.** Its radiocarpal is a
one-DOF `PinJoint` — flexion only, **no deviation coordinate at all** — and both
wrist coordinates ship `<locked>true</locked>` at −1.57079632 rad, the push-up
pose. Its `<credits>` is the placeholder `Model authors names..`.
`docs/UPPER_BODY_ACTUATION.md` §5 and §10.3 already ruled it out as a citable
source and nothing here changes that.

**Everything else has the joint and not the muscle.** Twenty-five further files
declare `wrist_flex_*`/`wrist_dev_*` — the whole Rajagopal family, the
Hamner/Delp gait family, this repository's own four variants — and in every one
of them, checked by parsing each muscle block rather than by reading coordinate
names, **not one muscle path touches a hand or carpal body.** `hand_{r,l}` is an
inert end segment.

## 2. The licence, from the rights holders, verbatim

The SimTK project **Wrist Model** (`simtk.org/projects/wrist-model`, group 325)
is run by the three people the model credits: Roger Gonzalez, Thomas S Buchanan,
Scott Delp. Its licence field reads **"Custom Use Agreement"**, the package is
`WristModel.zip` (73 KB, released 2008-07-25, last updated 2008-08-19), and the
agreement attached to it, read live on 18 September 2026 and reproduced in full
(sha256 of the raw string `51741011d62f393ae91ed347463238634ea0c038bced97c8529f4508bd23aa92`):

> **Wrist Model**
>
> Thank you for your interest in using our model of the wrist. You are welcome
> to use and modify this model for our noncommercial research. We hope it will
> be helpful to you and that you will continue to improve and apply the model to
> new research applications.
>
> There are five conditions on the use of this model.
>
> 1. This model was developed and tested over several years of collaborative
> work. An overview of the model is published in the references cited on this
> website. You must acknowledge the source of the model in any presentation,
> publication, or grant application that includes an image or description of the
> model or results generated with the model.
>
> 2. You must test the accuracy of the model in the context of your application.
> We have developed the model for a specific purpose and understand its strengths
> and limitations in this specific application. You are responsible for testing
> its accuracy in the context of your planned application. We make no guarantees
> as to its suitability for your particular application.
>
> 3. We ask you to provide improvements you make to the model to us, so that we,
> and others, may benefit from your work.
>
> 4. **You may not copy or distribute this model.** If others are interested in
> using the model, please direct them to this website.
>
> 5. **This model may be used only for non-commercial, academic work. It may not
> be used in any commercial activity.** You may not sell the model or results or
> images generated with the model.

*Issuer and authority:* the project team is the model's three authors, and it is
the only statement found that is issued by anyone who holds rights in it. The
bold is added here; every word is otherwise as published.

**Verdict: NOT USABLE, and by a wider margin than MoBL-ARMS.** MoBL-ARMS's
notice is BSD-3 restricted to non-commercial use, which at least *permits
redistribution with the notice*; this one forbids copying and distribution
outright, and forbids commercial use in terms that leave nothing to read around.
`docs/UPPER_BODY_ACTUATION.md` §11 records that the girdle donor was chosen
because its two candidate statements **converge** on commercial use with
attribution. Here there is one statement and it diverges from what this
programme would need.

**Nothing was taken.** No muscle, no path point, no parameter and no number from
this model appears in any file this work produced. The model file is read-only
on disk and `data/raw/` is gitignored — `git ls-files data/raw` returns **zero
tracked files** — so this repository does not redistribute it either.

## 3. The Apache copy is not a second opinion, and that is the recorded trap

The obvious move is to take the byte-equivalent copy from `opensim-core`, which
has a real `LICENSE.txt` (verbatim Apache-2.0) and would look clean. **That is
the exact error `docs/UPPER_BODY_ACTUATION.md` §10.1(4) recorded for
MoBL-ARMS** — reading a repository's own licence as covering a model somebody
else contributed to it — and the repository itself says not to. `opensim-core`'s
`NOTICE` ends:

> If you use plugins, models, or other components contributed by your fellow
> researchers, you must acknowledge their work as described in the license that
> accompanies each of these files.

No licence accompanies these files. Worse, the two `opensim-core` copies have
**placeholder credits**, so they do not even carry the attribution the authors'
condition 1 requires; a user of those bytes alone would not know whose model it
is.

**The asymmetry inside `opensim-models` is the other half of the argument, and
it points the same way.** That checkout has **no `LICENSE`, `NOTICE` or
`COPYING` file anywhere in the tree**; its only `README.md` is two lines and two
AppVeyor badges. Yet 33 of its `.osim` files carry an explicit in-file grant.
Arm26, two directories away, says so in its own `<credits>`:

> License: Creative Commons (CCBY 3.0). You are free to distribute, remix, tweak,
> and build upon this work, even commercially, as long as you credit us for the
> original creation. http://creativecommons.org/licenses/by/3.0/

`WristModel/wrist.osim` carries no such line: a grep for
`licen|copyright|creative|commons|BSD|MIT|apache` across all four donor files
returns **zero** hits. In a tree where models that are granted say so, silence is
not an inheritance. It is silence.

## 4. The donor that WOULD be clean, and why it is not here

`docs/UPPER_BODY_ACTUATION.md` §10.3 recorded that the Stanford-VA upper limb
model (Holzbaur, Murray, Delp 2005) is the one upper-limb musculature among
these sources released **without** a non-commercial restriction. Re-read live on
18 September 2026 from `simtk.org/frs/?group_id=324` ("Upper Extremity Kinematic
Model", package *Upper Extremity Model*, sha256 of the raw licence string
`037e1fe374fc9a1f6006d06e8c0a1bf93cff274ed019ac9f551152ff5a3caa7b`):

> The Stanford-VA upper limb model has been open sourced under the BSD 3-Clause
> License below. By downloading or using this software, (1) you accept the terms
> and conditions of the aforementioned open source license, (2) accept that use
> of the model software must be acknowledged in all publications, presentations,
> or documents describing work in which the Stanford-VA upper limb model is used
> by citing the following work: Holzbaur KRS, Murray WM, Delp SL. A model of the
> upper extremity for simulating musculoskeletal surgery and analyzing
> neuromuscular control. Annals of Biomedical Engineering 2005; 33: 829-840.
>
> Copyright (c) 2005, Stanford University and VA Palo Alto Health Care System.
> All rights reserved.

Then the standard BSD-3 text. **No non-commercial clause anywhere**, and the
attribution requirement is a citation this repository already carries (Arm26 is
its educational reduction and names the same paper).

**It is behind the same login wall as every other SimTK package.**
`https://simtk.org/frs/download_confirm.php/file/1131/UpperExtremityModel.zip?group_id=324`
returns **197 bytes of HTML** redirecting to `/account/login.php?triggered=1` —
the identical wall `docs/SHOULDER_GIRDLE.md` §1.4 hit at 216 bytes for the
thoracoscapular package and §10.1 hit for MoBL-ARMS. **So the acquisition gap is
real and it is not an acquisition this agent can perform.** Whether the model
contains the wrist muscles at the fidelity needed cannot be established from
here either; the 2005 paper describes 50 muscle compartments spanning shoulder
to wrist, and the on-disk catalog note that it has "no inertial properties for
the bodies" does not matter, because this plant supplies its own radius, ulna
and hand inertia and would take only paths and parameters.

## 5. What the registration would cost once a donor is cleared — and it is cheap

Worth recording, because `docs/ARTICULATED_SPINE.md` scoped the wrist job
against MoBL-ARMS and concluded *"a station copy is not possible"*. Against a
**Rajagopal/Holzbaur-lineage** donor it very nearly is, and here is the
measurement that says so.

This plant's arm bodies are Rajagopal's own bones put through OpenSim's
`ScaleTool`, and `ScaleTool` records the factors it used. Taking the
`radius_hand_r` joint offset — the wrist itself — from
`RajagopalLaiUhlrich2023.osim` and multiplying componentwise by this model's own
`radius_r` `<Mesh><scale_factors>`:

| | generic | × scale factor | product | this plant | Δ |
|---|---:|---:|---:|---:|---:|
| x | −0.008797 | 1.0693436226 | −0.009407015848012 | −0.009407015848012 | **0** |
| y | −0.235841 | 1.0010835523 | −0.236096546057984 | −0.236096546057984 | **0** |
| z | +0.013610 | 1.0916949738 | +0.014857968593418 | +0.014857968593418 | **0** |

**Zero to the last bit on all three components.** So an attachment expressed in
the generic humerus, ulna, radius or hand frame maps into this plant by
componentwise multiplication by that body's own recorded factors — the same rule
`scripts/build_corrected_foot.py` check C3 already audits over 19 joint offsets.
No similarity fit, no landmark estimate, no residual. The corroborating detail:
this plant's wrist ranges, inherited from `RajagopalLaiUhlrich2023.osim`'s own
`radius_hand`, are ±1.22173048 and −0.43633231 … +0.61086524 rad, and the
Gonzalez model declares the *same two ranges* to six decimals — one lineage, two
files.

What would still be real work: the Gonzalez wrist is two coupled carrier bodies
plus welded carpals where this plant is a single `UniversalJoint`, so the two
DOFs have to be identified rather than copied; the donor law is Schutte1993, so
either the deprecated law travels or the parameters are re-expressed under a
current one and that is declared; and the left side is a reflection prior needing
the explicit validation `docs/SHOULDER_GIRDLE.md` §4 shows is not optional.

**None of that is the blocker. The licence is.**

---

# PART TWO — THE THORACIC JOINT

## 6. The verdict was re-asked, under a rule fixed before it ran

`docs/ARTICULATED_SPINE.md` blocked the thoracic joint because

> "The only candidates in the plant are the gait2392 trunk muscles. Their
> `torso` insertion (y = 0.11 m) sits 32 mm above the thoracic joint centre,
> which is itself inherited from a registration with 62 mm RMS proxy residual,
> so reassigning them to `thorax` is not decidable from this data."

That is still true of those six muscles. It was written before
`shoulder_girdle_v1` added **sixty** Seth 2019 donor muscles, **32 of which take
an attachment on `torso`**, so the candidate set changed and the verdict had to
be re-measured rather than inherited.

`scripts/measure_thoracic_reassignability.py` holds the rule, and it was
**committed before the script was first run** (`1643a94`). Both halves:

* **R1, the endpoint.** `thorax` is not "the upper trunk". It is exactly the 48
  structures `native_composition_v1/plan.json` debited — ribs 1–12, costal
  cartilages 1–7, manubrium, body of sternum, xiphoid, three intercostal layers,
  diaphragm. **It holds no thoracic vertebra**: T1–T12 stayed in the 17.938 kg
  residual core inside `torso`. So an attachment belongs on `thorax` only if its
  **nearest anatomical surface**, measured in the canonical frame against this
  body's own meshes, is one of those 48. Nearest wins, no tolerance. The
  muscle's name is never consulted — IHM-1 `CLAUDE.md`, *"check each ENDPOINT
  against the structure its own label names"*.
* **R2, the margin.** More than **sqrt(2) × 62 mm = 87.7 mm** above the joint
  centre. Both sides of that comparison carry registration error of the same
  order — the joint centre from the cervical registration, the attachment from
  the girdle map, which `docs/SHOULDER_GIRDLE.md` §9 measures at the same few
  centimetres — so the margin must clear about `sqrt(2)·sd`, not a bare `sd`.
  This is the half the original verdict failed, at 32 mm.

A muscle passes only if **every** torso attachment passes; anything split is
reported MIXED and not moved.

### What it found

**Ten muscles, five per side. Verdict: UNBLOCKED.**

| muscle | nearest surface | distance | next non-thorax bone | margin above the joint centre |
|---|---|---:|---:|---:|
| `seth_PectoralisMajorThorax_I_{r,l}` | body of sternum | 31.6 / 31.9 mm | clavicle, 126 mm | **176.6 mm** |
| `seth_PectoralisMajorThorax_M_{r,l}` | manubrium | 39.9 / 40.4 mm | clavicle, 49 mm | **278.8 mm** |
| `seth_PectoralisMinor_{r,l}` | 4th rib / intercostal | 62.9 / 63.9 mm | 132 / 136 mm | **192.0 mm** |
| `seth_SerratusAnterior_I_{r,l}` | external intercostal | 66.4 / 67.9 mm | 166 / 168 mm | **104.5 mm** |
| `seth_SerratusAnterior_M_{r,l}` | external intercostal | 72.2 / 73.5 mm | 95 / 97 mm | **228.4 mm** |

**The refusals are the part that shows the instrument works.** Trapezius,
rhomboid and levator scapulae land **0.4–5.7 mm from a vertebra or an
intervertebral disc** — exactly where a spinous-process or transverse-process
origin should be — and the vertebrae are in `torso`, so R1 refuses all of them
however high they sit (up to 488 mm). Latissimus dorsi comes back **MIXED**,
with its inferior point **6.4 mm BELOW** the joint centre, which is what a
thoracolumbar-fascia-and-iliac-crest origin should look like. Serratus anterior's
superior slip is refused because its torso point is nearest the **scapula**.

**Two corroborations nobody asked for.** `gait2392_ercspn` reads **32.2 mm**
above the joint centre — reproducing the "32 mm" `ARTICULATED_SPINE.md` quoted
from an entirely separate computation. And `arm26_BIClong`/`TRIlong` come back
nearest the **scapula** (15.7–37.5 mm), which is the re-registration
`docs/SHOULDER_GIRDLE.md` §8 lists as open. Two known answers the script was not
written to produce.

### Two things recorded against this measurement

* **The candidate set had a defect.** The substring `'scapula'` also matched
  *subscapularis* (a muscle), *subscapular vein* and *infrascapular region*.
  They polluted the non-thorax set. Fixed to bones only (`role == 'rigid_bone'`)
  and re-run: **no verdict moved**, 10 reassignable before and after, same ten.
  The rule was not touched; only the set of structures it searches. Fifty label
  matches are now rejected as not-a-bone and the count is printed.
* **Two of the five won R1 by 8–15 mm.** `PectoralisMajorThorax_M` beats the
  clavicle by **9.1 mm** and `SerratusAnterior_M` beats its runner-up by
  **14.4 mm** — well inside the 62 mm residual, so those two assignments are not
  established anatomy, they are a nearest-wins test passing narrowly. **They were
  not rescored**: the rule was fixed in advance and they pass it. The margins are
  reported here and in the model's own `thoracic_reassignability.json` so nobody
  has to take the verdict on trust. *What would settle them:* a registration of
  the girdle donor with a residual smaller than the 9 mm in question, which does
  not exist here.

## 7. The build, and what it is allowed to be

`data/models/thoracic_drive_v1`, from `shoulder_girdle_v1`. **Two changes and
nothing else:**

1. the ten muscles' torso attachments move from `torso` to `thorax`;
2. the `Thorax_r` and `Thorax_l` rib-cage wrap ellipsoids move with them.

No mass moves, no body or coordinate is added, no muscle parameter is touched,
and **no number enters that was not already in the base model.** The `thoracic`
joint's parent frame is a pure translation on `torso` with zero orientation and
its child frame is the `thorax` origin — the builder **asserts both** before
using them — so `p_thorax = p_torso − t` is exact and the whole thing is
arithmetic. The wrap moves because it *is* the rib cage, and leaving the rib
cage's surface on one body and its attachments on another would put them ~44 mm
apart at the declared range limit; `hat_ribs_scap.vtp` was moved to `thorax` on
the same reasoning by `build_articulated_spine.py`.

`diff` against a round-trip of the base model is **128 lines**: two wrap
ellipsoids and ten path points. Nothing else in 683 kB.

**Ranges and stops.** Nothing new is needed and nothing new was invented.
`thoracic_extension` (±0.261799 rad), `thoracic_bending` and `thoracic_rotation`
(±0.349066 rad) already carry a declared `<range>` and an `articulated_stop_*`
`CoordinateLimitForce` **at exactly that range**, inherited unchanged. The
builder asserts all three — and the first version of that assertion **fired**,
because `CoordinateLimitForce` stores its limits in OpenSim's degree convention
while `<range>` is radians. The assertion converts; the same
degrees-vs-radians convention, and the damping term that is *not* converted, are
`docs/ARTICULATED_SPINE.md`'s existing record and are unchanged here.

### Gates: 10 pass, 1 recorded FAILED

`scripts/verify_thoracic_drive.py`, 11 tests, 3.0 s.

| check | result |
|---|---|
| no body's mass, mass centre or inertia changed | pass, all 29 |
| the rib-cage wrap did not move **in the world** (body and coordinates must cancel) | pass, < 1e-15 m |
| every `PathWrap` still resolves to a wrap object some body owns | pass |
| exactly the ten measured muscles changed body | pass |
| every thoracic coordinate declares a range, carries a stop at it, rests inside it | pass, 3/3 |
| the builder is a function (run twice → identical bytes) | pass, 6 files |
| the model loads | 29 bodies, 158 muscles, 60 coordinates, 77.6122029 kg |
| **every path length unchanged at the rest pose** | **< 1e-12 m over 158 muscles** |
| **every moment arm about every non-thoracic coordinate unchanged** | **< 1e-12 m** |
| standing controls: `moment_arms` bit-equal on repeat; `soleus_r` about `ankle_angle_r` | **−0.0497080 m**; `ercspn`/`intobl`/`extobl` about `lumbar_extension` 42.69 / −52.82 / −62.79 mm |
| **the control that can fail:** the BASE plant must still read ZERO about every thoracic coordinate | **0.00e+00**, re-measured in the same run |
| the reassignment adds no motion the base plant does not have | **paired difference exactly 0.000e+00 at all 17 bodies** |
| ~~distal-only, unpaired, bar 1e-9 m~~ | **FAILED — `tibia_r` moved 2.949e-07 m** |

**The failed gate stays failed and the bar stays where it was written.** It was
the bar that was wrong, and wrong in a way this repository has already recorded:
`docs/SHOULDER_GIRDLE.md` §6 says *"the distal-only bar had to be measured, and
the first one was below the noise"*, because this plant's assembler re-solves two
knee couplers and two acromioclavicular point constraints whenever any
coordinate is set, and that re-solve nudges far bodies by ~1e-7 m whatever was
perturbed. The gate is marked as an expected failure in the test file so it can
never quietly pass, and the question it was *for* is answered by a **different
quantity**, pre-registered rather than rescored: the paired difference against
the base plant under the identical perturbation. The base plant reproduces the
2.949e-07 m tibia nudge **bit for bit** (variant 2.949e-07, base 2.949e-07,
difference 0.000e+00), along with `calcn_r` 1.367e-06, `hand_r` 8.386e-07 and
`radius_r` 1.434e-07. Head, cervical and thorax move 93.4, 69.6 and 14.1 mm.

## 8. What it can do, and the shape of it is the finding

`Σ Fmax·|r|` per direction, at the rest pose, engine-read
(`data/models/thoracic_drive_v1/thoracic_capacity.json`):

| coordinate | muscles crossing | + bound | − bound |
|---|---:|---:|---:|
| `thoracic_extension` | 10 | **903.7 N·m** | **0.0 N·m** |
| `thoracic_bending` | 10 | 370.2 N·m | 368.8 N·m |
| `thoracic_rotation` | 10 | 347.2 N·m | 343.1 N·m |

**The joint is one-directional and that is the headline, not the 903.7.** All
ten muscles pull the same way about extension and there is **no antagonist at
all** — because the antagonist is erector spinae, whose origin the measurement
put on **T12**, inside the torso, 32.2 mm above the joint centre and therefore
refused. Bending and rotation look balanced, but only by left–right mirroring:
the right-side muscles bend right and the left bend left. There is no opposition
on either side.

So the honest statement of what changed is: **the thoracic joint went from
topologically unmuscled to one-directionally muscled at a 900 N·m ceiling,
against a stop of 30 N·m/rad at ±15°.** Anything free to search will drive it
straight to that bound, which is the `NATIVE_JOINT_LIMITS.md` hazard arriving at
a new coordinate. `default_enabled` stays **false** and there is no controller.

**And this is a ceiling, not behaviour.** `Σ Fmax·|r|` is full activation at
optimal fibre length with no force–length, no force–velocity, no tendon state
and no activation dynamics; a per-axis sum is not a tension-feasible torque cone;
the arms are one pose's; nothing here integrates.
`docs/UPPER_BODY_ACTUATION.md` §2 attaches that caution to this exact quantity.

**One number to read sceptically.** The moment arms are **112–229 mm**, which is
very large for muscle. That is not anatomy, it is the reduction: the plant's
thoracic joint is a *single lumped thoracolumbar articulation* at T12/L1, so a
pectoral attachment 28 cm higher sits on a 28 cm lever. A real thoracic spine
distributes that over twelve levels with short arms at each. The lever is a
property of the model's own topology and should not be quoted as a measurement
of a human chest.

**The wrap-object worry was measured rather than argued, and it was small.** The
ellipsoid that moved is used by **19 muscles per side**, not five. Of the 28 that
keep both attachments where they were, only **2** acquire a thoracic moment arm
at all, worth **0.3–0.6 N·m**. That population is reported as its own row in
`thoracic_capacity.json` and is never added to the ten: a moment arm that comes
from a wrap surface alone is a much weaker claim than one that comes from an
endpoint.

---

## 9. The press-up: nothing changed, and the wrist is why the elbow number is not the whole story

`docs/SHOULDER_GIRDLE.md` §7 measured, at the rest pose against a body weight of
761.1 N: arms straight gives **3,265 N (4.29 body weights)**, bounded by the
scapula; at the worst-case lever it gives **345 N (0.45 body weights)**, bounded
by the **elbow's 43.7 N·m**.

**Neither number moves.** No wrist muscle was installed, and the ten thoracic
reassignments change no moment arm about `elbow_flex`, `arm_flex`,
`arm_add`, `arm_rot` or any `scapula_*` coordinate — measured, to better than
1e-12 m, in the gate above. **The elbow is still the limiting element.**

But the wrist changes what that sentence is worth, and in the direction of
making it weaker:

* The load path in a press-up is hand → radius → ulna → humerus → scapula →
  trunk. `docs/SHOULDER_GIRDLE.md` §7's whole point was that **with a girdle the
  question is bounded, because every newton reaches the trunk through muscle**.
  The wrist is now the one link in that chain with **no muscle in it at all**.
* Inside its declared range the wrist has **no elastic resistance whatsoever**.
  Measured from the model file: `wrist_flex_{r,l}` and `wrist_dev_{r,l}` carry a
  `CoordinateLimitForce` at ±70° and −25°…+35°, an explicit engineering constant
  at 30 N·m/rad with a 0.35 rad transition, and an
  `ExpressionBasedCoordinateForce` of `-1.0*qdot`, which is **viscous damping and
  not a spring**. There is nothing else. Inside the range the wrist is a free
  hinge, and the only thing that would ever hold a hand flat against the floor is
  a stop that by construction only acts *outside* the range the anatomy declares.
* So **the elbow's 43.7 N·m is an upper bound on a chain whose weakest link is
  not yet measurable.** It is the smallest number among the links that have
  muscle. It is not the smallest number in the chain.

That is the whole consequence, and it is a correction to how §7's answer should
be read rather than a new number. And it remains **capacity, not demonstrated
behaviour**: nothing integrates, no controller exists, and the body has still not
been shown to press itself up.

## 10. What this does not do

* **No wrist muscle, and none is coming from anything on disk.** §2.
* **`data/sources/sensorimotor/upperbody_sources.json`'s `wrist` entry records
  no licence field.** That absence was accurate — the file carries none — and it
  is now filled in from the primary source, with the verdict, in
  `data/sources/wrist_donor.json`.
* **The thoracic variant has no `linearization.npz`, no stance acceptance, no
  controller, no actuator port on any thoracic coordinate, and
  `default_enabled: false`.** Nothing in the app selects it;
  `ihm/assembly/embodied.py` hard-codes `engineering_stance_v1` when a controller
  is chosen. A caller can ask for it today through `augmented_registration=`.
* **No display change.** `data/derived/anatomy-segment-binding/` still assigns
  every anatomical entity to one of the original 22 segments, so the rib cage
  rides `torso` in the workbench whatever the mechanics say. Same open job
  `docs/ARTICULATED_SPINE.md` and `docs/SHOULDER_GIRDLE.md` §8 already list.
* **The erector spinae was not moved**, so the thoracic joint has no extensor.
  Giving it one needs the thoracic vertebrae to become part of a body above the
  joint, which is the *"no thoracic vertebrae as bodies"* item
  `docs/ARTICULATED_SPINE.md` already carries — a per-level chain, not a
  reassignment.
* **`equilibrium_excitations.json` is the base plant's** and names only its 98
  muscles, carried across unchanged as `shoulder_girdle_v1` carries it. A caller
  that drives "every muscle at its equilibrium excitation" leaves all sixty
  girdle muscles, and therefore all ten reassigned ones, at zero. That is the
  right default for an unaccepted variant and it is not an equilibrium claim.
