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
              (kinematic, as built; registration_foot_paths.json and
               registration_muscled.json add muscles -- see "Making the joints
               drivable" below, which also withdraws part of this document's
               explanation of the ankle)

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

Stop constants are `scripts/crawl.py`'s: 30 N·m/rad, damping 1.5 **written into the
`.osim` in OpenSim's own Nm/(degree/s)** — i.e. 85.94 N·m·s/rad, the same physical
damping the base plant's engine-built stops apply (see `NATIVE_JOINT_LIMITS.md`,
18 Sep) — transition
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
and the foot folds, and the ankle follows it. *(Partly withdrawn 18 Sep: with
real subtalar arms the ankle still leaves its range by 1.14 and 0.87 rad —
gate G-S below. The missing arm is at most part of the cause.)* The nine spine coordinates are the
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

## Making the joints drivable: why the arms were zero, and what changed

Added 18 September 2026. The kinematic registration above is **unchanged**: its
model file, its F2 gate and every number in the sections above still describe
it, and F2 still runs on it and still fails. The new plants are separate,
opt-in registrations in the same directory:

    registration_foot_paths.json   model.osim + muscle_paths.xml (foot paths only)
    registration_muscled.json      model_muscled.osim + muscle_paths.xml (+50 neck muscles)

    diagnose/install  .venv/bin/python -m scripts.refit_muscle_paths_articulated_spine --install
    transfer          .venv/bin/python -m scripts.transfer_neck_muscle_paths
    report            data/models/articulated_spine_v1/muscle_paths_report.json

### "Exactly zero" was two different problems

The sentence above ("fitted polynomials in the coordinates that existed when they
were fitted") is true and does not separate two cases that need opposite fixes.
Measured through OpenSim on the model's **own `GeometryPath`s**, wrap objects
included, against the shipped fitted set at the same poses, and cross-checked by
reading path-point bodies off the XML (two routes, no shared code; they agree
on all 21 coordinates probed):

| coordinates | muscles crossing in the geometry | fitted set | class |
|---|---|---|---|
| `subtalar_angle_{l,r}` | **11 per side**, peak 32.7 mm (perbrev) | 0 | **representation** |
| `mtp_angle_{l,r}` | **4 per side**, peak 25.3 mm (ehl) | 0 | **representation, and it predates this variant** |
| nine spine/neck, four wrist | **none** | 0 | **topology** |
| `ankle_angle_r` (control) | 11 | 11 | expressed |

*Representation*: the muscle's real path crosses the joint and the polynomial has
no argument for it. The 11 per side are tibpost, perlong, perbrev, fdl, fhl, edl,
ehl, tibant, soleus, gasmed, gaslat. All run tibia → calcn.

*The mtp finding is about the base plant, not this variant.* `mtp_angle` is free and
unlocked in `engineering_stance_v1`. The shipped path set was fitted upstream on a
body with the toes **welded**: `exampleMocoInverse.cpp:51` and `exampleMocoTrack.cpp:51`
apply `ModOpReplaceJointsWithWelds({"mtp_r","mtp_l"})` before
`ModOpReplacePathsWithFunctionBasedPaths`. `native_mechanical_stream.cpp` does not
weld them. So the identified plant has carried two free toe hinges that no muscle
can control since the engine was built. This is not fixed there (do not touch that
plant), but it is fixed in both new registrations here.

*Topology*: no muscle in the plant has a path point above `torso` or beyond `radius`.
No fit can change that. Only a new muscle can.

### Subtalar: two refits FAILED, and what is installed instead

`PolynomialPathFitter` discovers each path's coordinates from the geometry
(`findIndependentCoordinates`), so a refit on this model picks subtalar up
unaided. Every fit made here recovered all 11 + 4 arms to within 1 mm. The gate
that failed was the one that stops a refit from recovering subtalar by damaging
what already worked. It was written before the first fit ran, and its bar is 1.5×
the shipped set's own worst peak-arm error per coordinate, on 27 held-out frames:

| attempt | instrument | `knee_angle_r` worst | bar | verdict |
|---|---|---|---:|---|
| 1 | one refit of all 80 | gasmed_r 1.66 mm | 1.36 mm | **FAIL** |
| 2 (pre-registered after 1) | shipped for 58 unchanged muscles; mean of 4 NEW refits for the 22 | gaslat_r **4.50 mm** | 1.36 mm | **FAIL** |

Both are recorded and neither is rescored. The fitter's own log explains why
attempt 2 was worse. Of its four fits, two discarded **706 and 1,170 of 1,485**
sample rows as NaN. The one that broke the mean fitted on ~21% of its samples in
51 s, against ~190 s for the others, and put gaslat_r's knee arm 20.3 mm off. A seeded probe of
675 fitter-like poses through the model's own paths gave **0 NaN**, so the cause
lies in the fitter's own design and was not isolated. Even the three clean fits put
the right gastrocnemius knee error at 1.66, 1.64 or 0.21 mm: bimodal. **The fitter
is not a stable instrument for wrapped biarticular paths with an added dimension,**
and a third refit would have been a third draw.

**Attempt 3 (pre-registered after 2, installed) is not a fit.** `muscle_paths.xml`
is the shipped set with those 22 muscles *removed*. The engine replaces only the
paths a set names, so the 22 run on the model's own `GeometryPath`s. That is the
representation this plant already uses for its 18 arm and trunk muscles. 18 of the
22 are via-point paths. The 4 gastrocnemii carry two wrap objects each. The other 58 keep
the shipped coefficients string for string (M2). Its accuracy gates pass **by
construction**, because the candidate is the truth they are scored against, and the
pre-registration says in advance that they are not evidence.

### Neck: the donor's own muscles, moved by the donor's own recipe

MASI (`data/research/cervical/MASI_HMaleMuscle_HMaleMassDistr.osim`, 78 Thelen2003,
plain path points, no wraps) is the donor the variant's neck ranges already came
from. `data/research/cervical_registration/v2/recipe.json` already holds the rigid
transform of every donor body into this torso frame, and its `cerv7` and `skull`
origins **are** this variant's `neck` and `atlantooccipital` centres (to 1e-12 m,
N1). The body map is kinematic: donor `spine`/`torso` → `thorax`, `cerv1–7` →
`cervical`, `skull`/`jaw` → `head`. Nothing is authored. Force, optimal fibre
length, tendon slack and pennation are the donor's text verbatim. The donor
declares no other Thelen property and has no `<defaults>`, so OpenSim's class
defaults apply, which is how the donor itself runs. Forces are unscaled
50th-percentile male, as the gait2392 trunk insert's are.

Which muscles, by rule:

| | donor muscles (both sides) | why |
|---|---:|---|
| **transferred** | **50** | cross at least one of `neck`, `atlantooccipital` and touch no absent body. Includes sternocleidomastoid (sternal), splenius capitis and cervicis, semispinalis capitis and cervicis, longus colli and capitis, scalenes, longissimus, iliocostalis cervicis, suboccipitals, T-level multifidi |
| excluded, girdle | 10 | cleidomastoid, cleido-occipital, trapezius (clavicular, acromial), levator scapulae: no clavicle or scapula body here. **Blocked on the same girdle as the shoulder.** |
| excluded, internal | 18 | every point on the one lumped `cervical` body (C–C multifidi, obliquus capitis inferior, longus colli C1–C5, T2–T1). They cross no joint here and would add inert states |

**Known answer (N3): every transferred path length equals the donor's own**,
evaluated by OpenSim on the donor model at its reference pose, to **1.5e-11 m** over
all 50. A wrong transform direction, a wrong body origin or a wrong body map would
each move lengths by centimetres. (The donor file does not load in OpenSim 4 as
shipped, because it has `/` in 36 component names. The verifier evaluates a copy with
only those names changed and asserts that no other byte differs.)

**Measured, not gated (N4): signs against the donor.** The donor's arm about
`pitch2` is a generalized arm summed over seven coupled levels. Ours about
`neck_extension` is one joint at C7/T1. These are different quantities, so no bar on
their agreement is derivable. Where both exceed 0.1 mm, **160 of 168** pairs agree in
sign. The 8 that disagree are 4 left/right pairs: the scalenus medius and posterior
about neck extension, longus colli C5–T about neck rotation (donor |arm| ≤ 1.7 mm,
i.e. near-cancelling across levels), and longissimus capitis about head rotation
(4.4 vs 4.8 mm: the donor's joint is C2/C1, ours the occiput). **28 donor arms are
lost to the lumping.** The largest is splenius capitis C6→skull's 40.4 mm about the
lower neck, which our lumped `cervical` absorbs.

**In the engine (D)**, `registration_muscled.json`: 25 bodies, **148 muscles**, 48
coordinates, 77.6122029 kg, and a repeat `moment_arms` call is bit-equal. Controls:
`soleus_r` about `ankle_angle_r` **−0.04977 m** (now its own GeometryPath; the fitted
value was −0.0497). `ercspn`/`intobl`/`extobl` about `lumbar_extension` are
42.69 / −52.82 / −62.79 mm, reproducing `LUMBAR_SHOULDER_MUSCLE_COVERAGE.md`'s
table. What changed, all engine-read:

| muscle | `neck_extension` | `head_extension` | anatomical reading |
|---|---:|---:|---|
| sternocleidomastoid (sternal) | −41.4 mm | +3.2 mm | flexes the lower neck, slightly extends the head at the occiput |
| splenius capitis (T→skull) | +41.0 | +39.3 | extensor at both |
| semispinalis capitis (T→skull) | +28.1 | +70.6 | extensor at both |
| semispinalis cervicis (T→C3) | +31.5 | 0 | lower neck only |
| rectus capitis post. major | 0 | +50.2 | occiput only |

Every one of the 50 has an arm about the neck or head. The 22 foot muscles have
one about subtalar and 4 per side about mtp. For all 148 muscles, every
`thoracic_*` and `wrist_*` arm is below 1e-9 m.

**Capacity, derived from the model's own records** (Σ Fmax × positive arm at the rest
pose: full activation at optimal length, so an **upper bound**), against the
gravitational moment of head plus cervical with the neck horizontal (prone), from
the partition's own masses and centres:

| joint | extension capacity | demand | flexion capacity |
|---|---:|---:|---:|
| `neck` | 53.2 N·m | 9.53 N·m | 17.4 N·m |
| `atlantooccipital` | 50.8 N·m | 2.69 N·m | **1.12 N·m** |

The one weak direction is **flexion at the occiput**. The only donor muscle
flexing there is longus capitis, and 1.12 N·m is below the 2.69 N·m it would take to
lift the head from supine at that joint alone. Rectus capitis anterior and
lateralis are not in the donor. Holding the head up has not been demonstrated
dynamically. This is capacity, not behaviour.

### GATE G-S: the foot paths do NOT repair the ankle — FAILED

Pre-registered before any plant carrying the foot paths was integrated. The
protocol and bar are F2's, unchanged: supine, 0.02 tonic, 2 s, and the base
model's own worst excursion measured in the same run (0.2202 rad). The plant is
`registration_foot_paths.json`, so the result is attributable to the foot paths
alone.

| coordinate | base | kinematic variant | foot paths | muscled |
|---|---:|---:|---:|---:|
| `ankle_angle_r` | 0.1229 | 1.3563 | **1.1399** | 1.1366 |
| `ankle_angle_l` | 0.1207 | 1.3522 | **0.8709** | 0.8710 |
| `subtalar_angle_r` | — | 0.1263 | 0.1865 | 0.1866 |
| `subtalar_angle_l` | — | 0.1716 | 0.1889 | 0.1888 |
| `thoracic_extension` | — | 0.3940 | 0.3970 | 0.3979 |
| `neck_extension` | — | 0.2084 | 0.2095 | 0.2130 |

**VERDICT: FAIL on both ankles, 5.2× and 4.0× the bar.** The muscles now have real
subtalar arms, and giving them those arms removes 16% and 36% of the ankle excursion,
not the regression. **So "the plantarflexors load a hinge they cannot control" is at
most part of why the ankle collapses.** The section above offered it as the cause;
that is now withdrawn as a sole explanation, and what else causes it is **not
measured**. The separating control not yet run is the kinematic variant with *only*
the subtalar re-welded: if the ankle recovers, the cause is the subtalar's freedom
itself (for example, contact geometry about a new axis). If it does not, the cause
is elsewhere in the variant. The neck muscles at 0.02 tonic do not reduce supine
neck excursion either. That is expected, since a tonic 2% drive is not posture
control, and it is reported rather than read as a finding.

Wall clock per 10 ms advance, same run: base **0.291 s**, foot paths **0.405 s**,
muscled **0.412 s**. The 22 GeometryPaths cost more than the 50 neck muscles. This
machine was at load 30–40 during the run, so the ratio is the usable number. No bar
was set, because none has provenance here.

### What is still blocked, and on what

* **`thoracic_*`: topology class, no donor muscle transferred.** The only candidates
  in the plant are the gait2392 trunk muscles. Their `torso` insertion (y = 0.11 m)
  sits 32 mm above the thoracic joint centre, which is itself inherited from a
  registration with 62 mm RMS proxy residual, so reassigning them to `thorax` is not
  decidable from this data. It needs a thoracolumbar donor with per-level
  attachments. None is registered here.
* **Wrists: topology class, blocked on MoBL-ARMS 4.1, but not on the girdle.** The
  wrist muscles (ECRL, ECRB, ECU, FCR, FCU, PL) originate on the humerus, ulna and
  radius and insert on the hand, so they need no scapula or clavicle. They do need a
  per-body registration of MoBL's humerus/ulna/radius/hand onto this model's:
  the frames are **not** the same frames. The elbow joint frame is rotated
  (−0.023, 0.228, 0.005) rad here and identity in MoBL, and MoBL's wrist is two
  coupled joints through a `proximal_row` body where this one is a single
  `UniversalJoint`. So a station copy is not possible. **It is the same donor as the
  shoulder, with the same unresolved licence** (`docs/UPPER_BODY_ACTUATION.md` §5),
  and a separable registration job. Forearm rotation stays effectively unactuated
  (`pro_sup` 2.16 mm, the Arm26 fusion) until the same job is done.
* **Girdle-anchored neck muscles** (10): the same girdle as the shoulder.
* **Lower-neck extension lost to lumping** (splenius capitis 40 mm about `pitch2`,
  and 27 more). This needs the full C1–C7 chain, which is the recipe's own
  `required_before_native_acceptance`.
* **The ankle regression**: cause open (G-S above).

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
