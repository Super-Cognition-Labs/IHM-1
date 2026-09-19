# the body can nod, and it still cannot hold its head up

`data/models/articulated_spine_v1` is a VARIANT of `engineering_stance_v1` that
adds the articulation `docs/WORKBENCH_AUTHENTICITY.md` item 0.2 says is missing:
a head, a neck, a thoracic segment, both wrists and both subtalar joints.
Built and measured 18 September 2026.

`engineering_stance_v1` is untouched. It is the identified plant every existing
result was measured on, and a gate is never rescored under a changed instrument.
This variant carries **no `linearization.npz`** for the same reason: the stance
plant was solved on 22 bodies and 33 coordinates and does not transfer to 25 and
48.

    build     .venv/bin/python -m scripts.build_articulated_spine
    verify    OPENBLAS_NUM_THREADS=1 nice -n 10 \
                  .venv/bin/python -m unittest scripts.verify_articulated_spine -v
    select    NativeMechanicalStream(..., augmented_registration=
                  'data/models/articulated_spine_v1/registration.json')

## What is there now

| | `engineering_stance_v1` | `articulated_spine_v1` |
|---|---|---|
| bodies | 22 | **25** |
| coordinates | 33 | **48** |
| `WeldJoint` | 4 | **0** |
| muscles | 98 | 98 |
| upright contact elements | 28 | **31** |
| total mass | 85.269848541731 kg unscaled | **identical to 1e-9** |

    torso ──thoracic(3)── thorax ──neck(3)── cervical ──atlantooccipital(3)── head
    radius_hand_{l,r}   WeldJoint → UniversalJoint   wrist_flex, wrist_dev
    subtalar_{l,r}      WeldJoint → PinJoint         subtalar_angle

Fifteen new rotational coordinates, every one with a declared `<range>` and a
`CoordinateLimitForce` **in the model file** at exactly that range. The stop
travels with the model rather than with the caller, because a stop a caller can
forget is the trap `docs/NATIVE_JOINT_LIMITS.md` is about.

## Where every number came from

Nothing here is typed. The torso was **repartitioned, never added to**, using
work this repository had already measured and left as `inertial_inputs_only`:

* `data/research/cervical_inertia/v2/manifest.json` — convex-envelope inertial
  priors for C1–C7, skull and jaw, already expressed in the current torso frame
  and already debited from the mass-scaled torso.
  Method and limits: `docs/research/CERVICAL_INERTIAL_PRIORS.md`.
* `data/research/thoracic_mechanism/native_composition_v1/plan.json` — the same
  torso's 48 measured thoracic material shares (ribs 1–12, costal cartilages
  1–7, manubrium / body / xiphoid, intercostals, diaphragm) and the residual
  core that is left when those and the nine cervical bodies are taken out.

Those records are in **mass-scaled** units (the engine multiplies every body mass
by `target_mass_kg / 85.26984854173146` at load); they are divided back out so
the XML carries unscaled masses and the engine's own scaling reproduces them.

| body | scaled kg | unscaled kg | COM in its own frame (m) |
|---|---:|---:|---|
| `torso` (was 27.654677 / 30.383239) | 17.937704 | 19.707537 | −0.02717, 0.25837, −0.00143 |
| `thorax` | 2.298053 | 2.524792 | 0.02386, 0.12559, 0.00109 |
| `cervical` (C1–C7 lumped) | 3.418920 | 3.756249 | −0.01655, 0.06346, 0.00050 |
| `head` (skull + jaw lumped) | 4.000000 | 4.394662 | −0.00762, 0.06865, 0.00030 |

**Joint centres**, torso frame, each read from a measured record rather than
placed by hand:

| joint | m | from |
|---|---|---|
| `thoracic` | (−0.02249, 0.07777, 0.00119) | midpoint of the T12 and L1 vertebral centroids (`anatomy.json` FJ3156/FJ3157) through the cervical manifest's own canonical→torso transform |
| `neck` | (−0.03454, 0.37427, 0.00205) | the composition plan's `cerv7` body-frame origin, i.e. the donor's C7/T1 |
| `atlantooccipital` | (−0.03505, 0.49413, 0.00262) | the plan's `skull` body-frame origin, i.e. the donor's atlanto-occipital level |

All three new body frames are **axis-aligned with the torso**. Positions come
from the registration and axes come from the torso; the registration's own
10.65° tilt is deliberately not carried into the joint axes, so the new joints
use the `back` joint's convention unchanged (rot1 z = extension, rot2 x =
bending, rot3 y = axial rotation).

**Ranges**, every one a donor's own declaration except the thoracic:

| coordinate | rad | source |
|---|---|---|
| `neck_extension` | −0.575959 … +0.837758 | MASI `auxt1jnt` `pitch2` |
| `neck_bending` | ±0.575959 | MASI `roll2` |
| `neck_rotation` | ±0.471239 | MASI `yaw2` |
| `head_extension` | −0.279253 … +0.418879 | MASI `aux2jnt` `pitch1` |
| `head_bending` | ±0.104720 | MASI `roll1` |
| `head_rotation` | ±0.663225 | MASI `yaw1` |
| `wrist_flex_{l,r}` | ±1.221730 | `RajagopalLaiUhlrich2023.osim`'s own `radius_hand` |
| `wrist_dev_{l,r}` | −0.436332 … +0.610865 | same |
| `subtalar_angle_{l,r}` | ±0.610865 | `RajagopalLaiUhlrich2023.osim`'s own `subtalar` |
| `thoracic_*` | ±0.261799 / ±0.349066 | **NOT MEASURED — see below** |

MASI is `data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim`, the donor
the cervical registration work already acquired and audited. Its two FREE
cervical joints are `auxt1jnt` (T1–C7) and `aux2jnt` (C2–C1); the other five
levels are driven by couplers. This model lumps C1–C7 into one body, so the
upper-cervical ranges are applied at the atlanto-occipital joint rather than at
C1/C2. That is a reduction, and it is why the two joints between them permit
−49°…+72° flexion/extension, ±39° lateral bending and ±65° axial rotation
instead of the donor's seven-level chain.

**The one number with no measured provenance** is the thoracic range, and it is
labelled as such in the builder. No donor here declares one: MASI welds its rib
cage (`ribcagejnt`), the Rajagopal family has no thoracic coordinate, and
`data/research/thoracic_mechanism/v2/manifest.json` states its own 0.05 rad rib
bound as "a declared kinematic exploration bound, not physiological ROM". The
values are deliberately narrower than any published whole-thoracic range, because
narrow is the conservative side for a stop whose only job is to keep a search out
of absurd configurations. *What would replace it:* a per-level thoracic ROM
source, or a thoracic coordinate measured from recorded motion — and the pose
corpus cannot supply the second, because it is in gait2392 space and gait2392 has
no thoracic degree of freedom either.

Stop constants are `scripts/crawl.py`'s: 30 N·m/rad, damping 1.5, transition
0.35 rad — the only stiffness this repository has swept
(`docs/NATIVE_JOINT_LIMITS.md`). Viscous damping is transferred from the source
model's own `ExpressionBasedCoordinateForceSet` by limb: every upper-body
coordinate there carries `-1.0*qdot` and every lower-limb one `-0.1*qdot`.

## What the known-answer checks print

`scripts/verify_articulated_spine.py`, 15 tests, 43-49 s (two 2-second
native runs dominate it). **14 pass, 1 is a
pre-registered gate that FAILS and is recorded as failed.**

| check | result |
|---|---|
| total mass conserved against `engineering_stance_v1` | **0 to 1e-9**; no body but `torso` changed, and what `torso` lost equals what the three new bodies gained |
| the four parts recombine to the original ONE torso body | mass **3.6e-15 kg**, COM **< 1e-12 m**, inertia Frobenius **6.7e-16 kg·m²** |
| every new body's inertia admits the engine's inscribed contact sphere | pass, all 25 |
| every new coordinate declares a range, and zero is inside it | pass, 15/15 |
| every new coordinate carries a stop AT that range | pass, 15/15; and **no existing coordinate gained one** |
| the model loads in the engine | 25 bodies, 98 muscles, 48 coordinates, total mass 77.6122029 kg to 1e-9 |
| each new joint moves only its distal entities | pass; smallest distal displacement **9.2 mm**, worst proximal **1.46 µm**, a factor of 6,300 |
| the rest pose is strictly inside every new range | pass, all 15 read exactly 0 |
| no muscle has a moment arm about any new coordinate | **0 exactly**, for all 5 probed muscles × 15 coordinates |
| the builder is a function (run twice → identical bytes) | pass |
| **GATE: no new coordinate leaves its range by more than the base model leaves its own** | **FAIL** |

Two of those deserve the detail.

**The distal-only test failed first, for the reason a control is supposed to
fail.** It probed each body's ORIGIN, and the origin of `thorax`, `cervical` and
`head` IS its own proximal joint centre — so a 15° rotation moved it by 8e-12 m
and the test read "the joint does nothing". A pivot does not move under its own
rotation. It now probes a station at (0.05, 0.05, 0.05) m, off every rotation
axis. The proximal bar also had to move: it is **not** zero, because this model
carries `CoordinateCouplerConstraint`s on both knees, so OpenSim re-assembles the
whole body whenever any coordinate is set and the solver's answer depends on
where it started — 8e-12 m for the spine joints, up to 1.46e-6 m for the wrist
and subtalar. A bar at 1e-12 would have failed on assembly noise and read as a
topology error.

**The gate, in full.** Protocol taken unchanged from
`docs/NATIVE_JOINT_LIMITS.md`: supine, 0.02 tonic excitation on every muscle,
2 s, worst excursion past each declared range. The bar is the **base model's**
own worst under the identical protocol, measured in the same run.

    base     0.2202 rad  knee_angle_r   (knee_angle_l 0.2195, pro_sup 0.1809,
                                         ankle_angle 0.1229)
    new      thoracic_extension 0.3940   ← over the bar, 1.79x
             thoracic_bending   0.2099
             neck_extension     0.2084
             subtalar_angle_l   0.1716
             head_extension     0.1261
             subtalar_angle_r   0.1263
             neck_bending       0.0829
             the other eight stay inside their ranges

**VERDICT: FAIL on `thoracic_extension`.** The threshold is not moved and the
stiffness is not tuned; the number was seen after the bar was fixed, and changing
either afterwards is indistinguishable from choosing the one that passes.

The bar had to be measured on the base model, not read off the variant's own
non-new coordinates — that would have been the wrong population, and a
spectacularly wrong one: see the next section.

## What this variant costs, measured

**Un-welding the subtalar makes the ankle far worse.** Same protocol, same
excitation, base against variant:

| coordinate | base excursion | variant | change |
|---|---:|---:|---:|
| `ankle_angle_r` | 0.1229 | **1.3563** | +1.2334 |
| `ankle_angle_l` | 0.1207 | 1.3522 | +1.2315 |
| `lumbar_extension` | — | — | +0.4001 |
| `hip_adduction_r` | — | — | +0.1413 |
| `pelvis_tilt` | — | — | +0.1407 |

Sixteen existing coordinates get worse by more than 0.05 rad. The cause is the
next finding, and it is the most important sentence in this document.

**No muscle in this plant has a moment arm about any of the fifteen new
coordinates. Measured, exactly zero.** The 80 source muscle paths are fitted
polynomials in the coordinates that existed when they were fitted
(`subject_walk_scaled_FunctionBasedPathSet.xml`, substituted in by
`ModelFactory::replacePathsWithFunctionBasedPaths`), so **un-welding a joint does
not give the muscles a moment arm about it**. `soleus_r` still has −0.0497 m
about `ankle_angle_r`, and 0 — the integer, in the engine's own JSON — about
`subtalar_angle_r`. The 18 added Thelen muscles keep a real `GeometryPath`, but
none of them spans a new joint either: `gait2392_ercspn_r` has +0.0427 m about
`lumbar_extension` and 0 about `thoracic_extension`.

So the subtalar is a free hinge that the plantarflexors load and cannot control,
and the foot folds, and the ankle follows it. The nine spine coordinates are the
same: a 4.39 kg head on a joint with a soft stop, viscous damping and nothing
else. **This variant gives the body joints, not the ability to use them.** The
muscles that would use them are `docs/UPPER_BODY_ACTUATION.md`'s subject, not a
model file's: this plant has 18 muscles above the pelvis (six lumbar, twelve
arm26) and none of them spans a neck, a wrist or a subtalar.

**Wall clock**, supine, 2 s at 0.02 tonic: base **0.123 s** per 10 ms advance,
variant **0.144 s** (+17%). Upright at default excitation: **0.160 s**. Fifteen
new free coordinates cost less than expected, which is the one pleasant number
here.

**Contact.** The upright environment inscribes one sphere in each body's inertia
ellipsoid at its mass centre, so three new bodies are three new balls (16 → 19
proxies, **28 → 31 contact elements**; the other 12 are the source foot
contacts). The new `torso`'s own proxy also grew, from **0.2577 m to 0.3090 m**
radius, because taking the head and rib cage out left a rounder residual; the new
balls are 0.1356 m (`thorax`), 0.0599 m (`cervical`) and 0.0758 m (`head`).
This is a real change to the contact model and it is a CONFOUND on the gate
above: in supine, part of the thoracic moment is the new head and thorax balls
pressing on the plane, not the stop being soft. The two are not separable with
this contact model.

## What I could not build, and why

* **No scapula, no clavicle, no shoulder girdle.** The thoracic material
  partition is ribs, cartilages, sternum, intercostals and diaphragm — 2.298 kg
  of rib cage. The scapulae and clavicles are not in it, so this partition has no
  mass to give them. Both `acromial` joints, the two `arm26_TRIlongglen` wrap
  cylinders and the six arm26 muscle attachments therefore still sit on `torso`,
  not on `thorax`; moving them onto a `thorax` whose inertia does not include the
  girdle would have put the arms on the wrong body.
  `docs/UPPER_BODY_ACTUATION.md` warns that "doing the neck first and the girdle
  later means partitioning torso mass twice", and this variant is the neck-first
  half of exactly that. It is not a problem: the partition is done with
  `ihm.assembly.cervical_inertia.partition_body`, which debits mass, first moment
  and full inertia exactly (reconstruction residual 3.6e-15 kg and 6.7e-16 kg·m²
  here), so a second debit of MoBL-ARMS' donor clavicle (0.156 kg) and scapula
  (0.70396 kg) per side composes with this one instead of competing with it. The
  girdle must come out of `torso`, whose 17.938 kg residual core still contains
  it.
* **No fingers.** `hand_{l,r}` is still one rigid body; the variant gives it a
  wrist, not a hand.
* **C1–C7 are lumped into one `cervical` body, and the jaw into `head`.** The
  measured records for all nine exist and the full chain is registered
  (`data/research/cervical_registration/v2/recipe.json`: 24 coordinates, 18
  coupler constraints, 78 muscles). Building it means converting a SIMM-era
  donor's `CoordinateCouplerConstraint`s and 78 muscle paths into this model's
  schema and re-anchoring their frames, which is the recipe's own
  `required_before_native_acceptance` list and is a separate job. Two joints and
  two bodies is the reduction that fits in one model file.
* **No thoracic vertebrae as bodies.** The thoracic joint is a single lumped
  thoracolumbar articulation, not a T1–T12 chain.
* **No actuator port on any new coordinate.** Adding `CoordinateActuator`s would
  have given the new joints drivable torque ports, but a port held at zero is the
  same torque as no port, and the 13 declared ports are counted elsewhere. The
  real answer is muscles, above.
* **No `linearization.npz`, no stance acceptance, `default_enabled: false`.**
  Nothing in the app selects this variant; `ihm/assembly/embodied.py:271`
  hard-codes `engineering_stance_v1` when a controller is chosen. Selection IS
  parameterised — `augmented_registration=` reaches `NativeMechanicalStream`
  through `ArticulatedBodyPlant` and `EmbodiedRuntime.from_workspace` — so a
  caller can ask for this model today; the app's own default cannot, and changing
  that is someone else's file.
* **The display still rides `torso` in the app.** The `.osim`'s own
  `attached_geometry` was moved (`hat_skull.vtp` and `hat_jaw.vtp` to `head`,
  `hat_ribs_scap.vtp` to `thorax`, each behind an offset frame that puts the
  torso-frame vertices back where they were), but the workbench poses anatomy
  through `data/derived/anatomy-segment-binding/`, which assigns every entity to
  one of the 22 original segments. Until that binding knows about `head`,
  `cervical` and `thorax`, the skull in the display still rides `torso`.

## One thing found on the way, which is not about this model

`scripts/native_mechanical_stream.cpp:293` builds its opt-in stop as

    new CoordinateLimitForce(name, upper*to_degrees, stiffness/to_degrees,
                             lower*to_degrees, stiffness/to_degrees,
                             damping, transition*to_degrees, true)

The limits, the stiffness and the transition are all converted to OpenSim's
degree convention. **`damping` is not.** `CoordinateLimitForce` reads rotational
damping in N·m/(degree/s), so `scripts/crawl.py`'s declared
`damping_nm_s_per_rad: 1.5` is applied as 1.5 N·m/(deg/s) = **85.9 N·m/(rad/s)**,
57.3× the stated value. It is a damping term inside a stop that is already an
explicit engineering constant, so nothing published rests on it — but the number
in `execution.json` is not the number the plant used, and a sweep over damping
would have been reading the wrong axis. Not fixed here: the engine is a shared
binary and changing it changes the plant for every caller and every prior run.
The same constants are written into this variant's own
`CoordinateLimitForce`s **with the same convention**, so the two paths agree with
each other and both disagree with their declared units.
