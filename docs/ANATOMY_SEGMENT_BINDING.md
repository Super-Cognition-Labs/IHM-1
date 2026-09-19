# the anatomy and the segments the controller drives

`docs/DISCONNECTS.md` item 1, in IBM-1:

| | |
|---|---|
| declared / rendered | 3,816 BodyParts3D entities in the simulated body |
| actually driven | 22 OpenSim bodies, 80 muscles |
| mapping between them | **none** |

There is now a mapping, and everything below is measured on it.

    scripts/bind_anatomy_to_segments.py   -> data/derived/anatomy-segment-binding/binding.json
                                          -> data/derived/anatomy-segment-binding/report.json
    scripts/render_anatomical_motion.py   -> a render, an app-format trajectory, a motion report

`data/derived/` is gitignored. The scripts are the artifact; regenerating the
binding takes about a minute.

## What the binding is

Each of the 4,000 atlas entities is assigned one of the 22 OpenSim bodies.

334 of them are bones and are assigned by the atlas' own label — "right femur"
is `femur_r`, every rib and vertebra and the whole skull is `torso`, because
`torso` is the model's HAT body and carries `hat_skull`, `hat_jaw`,
`hat_spine` and `hat_ribs_scap`. Those 334 bones then *are* the reference
shape of each segment.

The other 3,666 are assigned geometrically: each entity's surface is sampled,
every sample takes the segment whose bone group is nearest, and the entity
takes the majority. The vote share is kept as `coherence`. An entity that
straddles a joint has a low one, and a rigid binding will tear it — 220
entities are below 0.6 and the renderer does not pretend otherwise.

Five entities are **not** bound and are excluded by name and reason: the skin,
its three layers, and the lymphatic network graph. No single rigid segment
carries a whole-body surface. The skin already has a continuous
graph-regularised linear-blend attachment over the same 22 segments
(`data/derived/canonical/continuous_surface_binding.json.gz`) and the renderer
uses it, so the integument stretches across a joint instead of tearing at it.

**3,995 of 4,000 bound, 99.88%.**

## The registration, and what it costs

The OpenSim model is a scaled Rajagopal subject; the atlas is the BodyParts3D
adult male; and the model's default pose is not the atlas pose. So the
reference pose (33 coordinates) and one similarity (scale, rotation,
translation) are fitted together by least squares against the 22 bone groups.
Scale 0.963. What is left over:

| | |
|---|---|
| bone-group centroid residual | **24.7 mm RMS**, 51.8 mm max (humerus_r) |
| joint centre vs the anatomical joint | **38.4 mm median**, 101.4 mm max (shoulder) |

Neither is measurement error and neither is subject variation. It is the
residual of putting one specimen's skeleton onto another's.

The joint-centre offset is the number that matters at playback, because the
anatomy is rotated about those points. Its consequence is measured directly:
at the reference pose the 64 closest surface point pairs across a joint sit a
fixed distance apart; bind the two sides to different rigid segments and that
distance changes. The largest change over a clip is exactly the gap or
interpenetration the binding introduces.

## Gates

None of these is the assignment rule restated.

| gate | result |
|---|---|
| held-out bones — drop a bone from its own segment's reference shape and re-assign it | **317 / 322, 1.6% misassigned** |
| native muscle attachments — the bound segment must lie on the kinematic chain between the muscle's declared attachment bodies | **78 / 80, 2.5% misassigned** |
| hand-checked cases with an answer anatomy fixes | **21 / 21** |
| second opinion — same rule, run against the registered OpenSim bone meshes instead of the atlas' own bones | **93.5% agreement** |

Only 322 of the 334 bones can be held out; the other twelve are the only bone
in their segment, so removing one removes the segment. All five held-out
failures are bones that are more than 30% of their own segment's surface —
calcaneus, both hip bones, right tibia — where the same thing happens in
miniature.

Endpoint-only scoring gives the muscle gate 64/80, but that is the wrong test:
a two-joint muscle's belly lies on the bone *between* its attachments, and
`attachment_bodies` names the endpoints. Hamstrings bind to the femur,
gastrocnemius to the tibia; those are right, not wrong.

Both remaining misses are **psoas major**, and they are a real disagreement
rather than a binding error: OpenSim attaches psoas to the pelvis, and the
atlas originates it on the lumbar spine, which is `torso`. The anatomy is
right and the model's attachment list is the approximation.

## How far the anatomy actually moves

The comparison that matters is the body's own canonical run, which moves
**6.35 mm** at most across 30 s. It is perfusing and breathing, not moving.

| trajectory | duration | max entity travel | median | worst joint opening | median opening |
|---|---|---|---|---|---|
| `unified-world-ofb2dp7z/full` — trained cortical kernel, in the loop | 4.98 s | 262 mm | 100 mm | 44 mm | 5 mm |
| `gait-best` — the gait controller: one genuine step, then a fall | 1.49 s | 1,028 mm | 397 mm | 73 mm | 6 mm |
| `forced-gait-corpus` — measured human gait, played in place | 2.39 s | 1,193 mm | 1,001 mm | 69 mm | 15 mm |
| `unified-world-ofb2dp7z/sever` — cortex severed, falls | 2.88 s | 1,835 mm | 949 mm | 142 mm | 22 mm |

**Travel is not skill.** The severed run travels furthest of all because it
collapses, and the falling gait controller beats the standing cortical one.
The only thing this column establishes is that the anatomy is now driven at
all, which is what item 1 says it was not. Skill belongs to
`data/derived/gait-best/report.json` and to the ablations, not here.

## What this does not do

- One rigid segment per entity. There is no soft-tissue deformation, no
  volume preservation, no sliding, and the joint openings above are the price.
- The model has no neck, no shoulder girdle and no spine joint. Skull, jaw,
  teeth, every vertebra, every rib, both scapulae and both clavicles ride
  `torso` as one body, because that is what the model is. The head cannot nod.
- A replay carries no physiology and no internal state. `/api/body/segment-bound`
  lists these separately from the canonical trajectory and says so on every
  entry, so a kinematic replay can never read as the body's own physics.

## The per-frame call (18 September 2026)

`ihm/assembly/anatomy_pose.py` turns this binding into something a running simulation can call
every frame. `scripts/render_anatomical_motion.py` had used it only offline.

    poser = AnatomyPoser.from_workspace(root)          # ~0.5 s, checks frames and provenance
    pose = poser.pose_from_native(native_frame)         # or .pose({body: 4x4}), .pose_from_coordinates(q)
    pose.rotation[i], pose.translation[i]               # entity pose.entity_ids[i]: x = R x_rest + t
    skin = poser.skin_vertices(pose)                    # FJ2810, by the continuous linear blend

It costs 1.3 ms a frame, and 57 ms with all 102,467 skin vertices blended.
`scripts/verify_anatomy_pose.py` gates it, and every gate passes:
- FK against Simbody: 7.8e-16 over 12 native frames.
- The reference pose reproduces rest: 5e-16 m, centroids and vertices.
- Each of the 31 coordinates moves exactly its distal entities, with 0 violations in either pivot
  mode.
- Output is bit-identical when called twice.
- Each frame guard is made to fire: millimetres, the atlas frame, the 15.7% display-box rescale,
  a missing body, and `z-anatomy-display-normalized` vertices are all refused.
- Coverage: 3,995 rigid, the skin blended, its 3 layers following, and 1 not posed (the lymphatic
  network graph, listed on every frame).

Two corrections the module makes, and says so:
- **FK.** It evaluates `SimmSpline` as OpenSim does (Forsythe–Malcolm–Moler).
  `scripts/render_body_3d.py` uses a natural spline, and the patellae there sit up to 6.5 mm off
  Simbody.
- **Symmetry.** 11 of the 1,342 mirrored pairs in this binding sit on unmirrored segments. At load,
  each pair takes its higher-coherence side. All 11 overrides are in `poser.symmetry_overrides`,
  and the fix at source is to score both sides in `bind_anatomy_to_segments.py`.

`pivot='anatomical'` re-seats each joint at the closest bone surfaces instead of at OpenSim's
joint centre, keeping every orientation exact. It cuts the worst joint opening from 108 to 32 mm
at the knee and from 105 to 14 mm at the lumbar joint. In exchange it moves the anatomy off the
simulated segments by 11 mm median and 50 mm max over gait-best. The default is `opensim`.

The full account is in IBM-1 `docs/DISCONNECTS.md` §1.

## The 25-body variant: a second binding (18 September 2026)

`data/models/articulated_spine_v1` (`docs/ARTICULATED_SPINE.md`) adds `thorax`, `cervical` and
`head` and un-welds both wrists and both subtalars. The binding above knew only the 22 original
bodies, so in the display the skull, every vertebra and every rib still rode `torso`, and
`AnatomyPoser.check_bodies` refused the variant's 25 bodies outright. There is now a second
binding. **The 22-segment one is untouched**; every result above was measured on it.

    .venv/bin/python scripts/bind_anatomy_to_segments.py --variant articulated_spine_v1
        -> data/derived/anatomy-segment-binding-articulated-spine-v1/{binding,report}.json   (~50 s)
    .venv/bin/python scripts/record_articulated_spine_native_frames.py      (2 engine sessions, one at a time)
        -> .../native_frames.json
    poser = AnatomyPoser.from_workspace(root, 'articulated_spine_v1')       # default is unchanged
    .venv/bin/python scripts/verify_anatomy_pose.py --plant articulated_spine_v1

### What it is: a refinement, not a refit

* **The registration carries over exactly.** Where the variant's fifteen new coordinates are
  zero, its 22 shared bodies are the base body: their transforms at the base reference pose
  differ by **0.0** (the builder checks this and refuses to continue if it does not hold). So
  the similarity A (scale 0.963), the reference pose and every residual in the table above
  still apply. The fifteen new coordinates have a reference value of 0. That is not a fitted
  value. It is the pose the variant was built to equal the base body at.
* **Only what rode `torso` is re-voted.** An entity the base bound to any other segment keeps
  it. Hand, talus and calcaneus were already segments, so un-welding the wrists and subtalars
  needs no reassignment. A base-`torso` entity is re-voted with the same rule
  (`vote`: nearest bone group per sampled vertex, majority wins), restricted to the four groups
  that make up the torso family. A flat 25-way vote was computed as a comparison and not used.
  It can differ only on base-`torso` entities, and it would have sent **13** of them out of the
  family to a limb or the pelvis because splitting the torso's vertex share four ways lets a
  runner-up win. The 13 are the clavicular pectoralis, iliocostalis lumborum, thoracolumbar
  fascia layers and four mesocolic/omental structures.
* **Torso-family bones are named from the records that define the variant's mass.** These are
  the cervical prior (`data/research/cervical_inertia/v2/manifest.json`: skull = 15 cranial
  bones, jaw = mandible, cerv1–7) and the thoracic plan
  (`data/research/thoracic_mechanism/native_composition_v1/plan.json`: ribs, sternum parts).
  Two classes are named by anatomy and not by a record, and the binding says so (`named_basis`):
  * thoracic vertebrae go to `thorax`, because the thoracic joint sits at T12/L1;
  * teeth, ossicles and the small facial bones go to `head`.
  The hyoid is claimed by no record, so it is voted (it went to `head`, 0.74). Scapulae and
  clavicles stay on `torso`, because both acromial joints still hang from `torso` in the variant.

### How many entities moved where

Of the base binding's 2,440 `torso` entities (`report.json`, `base_torso_entities_now_on`):

| segment | entities | named bones among them |
|---|---:|---:|
| `thorax` | **1,021** | 43 |
| `head` | **885** | 89 |
| `torso` (stays) | **372** | 13 |
| `cervical` | **162** | 10 |

After the poser's load-time symmetrisation, the posed counts are thorax 1,022, head 885,
torso 372 and cervical 162. Symmetry overrides go from **11 to 26**. The 15 new ones are pairs
that the four-way vote split across the thorax/torso or head/cervical boundary, for example
the left and right common carotid, the adrenal glands and `rectus capitis lateralis`.

### Gates on the binding

| gate | result |
|---|---|
| held-out bones: drop each named torso-family bone from its group and re-vote among the four | **150 / 155** |
| mass-record consistency: every entity the variant's mass partition names rides the body its mass is in | named 50/50 (a restatement); **voted 21/21**, the real test: intercostals, diaphragm, costal cartilages 1–7 |
| hand-checked cases, written before the first build | **32 / 33, FAILED on one** |
| native muscles against the VARIANT's kinematic chain | 78 / 80, the same two psoas as the base |
| builder run twice with a pinned `IBM_GIT_SHA` | byte-identical `binding.json` |

Held-out misses: both scapulae (30–32% of their segment's surface) and both clavicles go to
`thorax`. The girdle lies on the rib cage but kinematically hangs from `torso`, and removing a
scapula leaves the torso group as mostly lumbar spine. C7 also goes to `thorax` (it is adjacent
to T1).

**Hand-checked miss, recorded FAILED:** `left sternocleidomastoid` was pre-registered as
`{cervical}` and went to `head`, with family coherence 0.33 and `torso` as runner-up. It is the
least coherent entity in the family. Its belly lies over the cervical spine, but its surface is
split three ways between the mastoid, the neck and the sternum. The acceptance set is not
widened.

### Where the new joints sit relative to the anatomy

The joint centre is carried through A and compared with the bones that form the joint.

| joint | offset | note |
|---|---:|---|
| thoracic | **8.1 mm** from L1/T12 | measured on whole segments it read **187.0 mm**, which is where the scapulae touch the ribs, not the joint |
| neck | 30.7 mm from T1/C7 | |
| atlantooccipital | 28.8 mm from C1/occipital | |
| subtalar r / l | **8.1 / 7.6 mm from the pin's AXIS** | 97.8 / 97.9 mm from its frame origin |
| wrist r / l | 7.5 / 8.9 mm | |

**Two instrument corrections came out of this table, and both apply to the base report too:**

* **Whole-segment proxies measure the wrong contact.** Torso and thorax touch most at the
  scapula–rib interface. The builder now names the bones that form each spine joint and the
  subtalar (talus with calcaneus, not navicular). It writes them to the binding as
  `joint_pivot_bones`, and `AnatomyPoser.pivots()` uses them. The base binding names none, so
  the base pivots are unchanged.
* **For a PinJoint, the distance to the frame origin is not the joint error.** Every point on
  the axis is the same pin. The base report's **101 mm** subtalar offset (and the 86–91 mm
  radioulnar ones, and the other pins) are point distances. The subtalar axis passes within
  8 mm of the talocalcaneal articulation.

**The foot.** The named split (talus → `talus`; calcaneus, navicular, cuboid, cuneiforms,
metatarsals, sesamoids → `calcn`) is the model's own geometry partition. `r_talus.vtp` rides
`talus_r` with 52 mm extent, and `r_foot.vtp` rides `calcn_r` with 212 mm extent from heel to
metatarsal heads. The OpenSim-mesh second opinion disagrees on 9 named foot bones: both calcanei
and both naviculars go to talus, and 3 sesamoids go to toes. That instrument is density-biased,
not evidence against the split. It measures distance to mesh VERTICES: 99 on the talus against
1,000 spread over the foot. The atlas' right calcaneus is 8.8 mm from the nearest talus vertex
and 15.7 mm from the nearest foot vertex, which is the same long-facet-tail hazard as CLAUDE.md's
centroid-seeded search. 17 soft entities ride a talus, all ligaments and tendon sheaths that
cross the ankle or subtalar joint (coherence 0.43–0.67). Each will open at the subtalar joint
now that it moves.

### The verified pose, both plants

`scripts/verify_anatomy_pose.py` runs the same battery on each plant. The variant runs it on its
own native frames: 12 frames, 2 sessions, every new coordinate off zero, one session on each side
of its range.

| gate | `engineering_stance_v1` | `articulated_spine_v1` |
|---|---|---|
| FK vs Simbody `transform_ground` | 7.77e-16, 12 frames | **9.99e-16**, 12 frames (both UniversalJoint wrists, both subtalar pins, all three spine joints) |
| reference pose = atlas rest | 5.05e-16 m | 5.05e-16 m |
| distal only, OpenSim / anatomical pivots | 31 coordinates, 0 violations | **46 coordinates, 0 violations**, both modes |
| idempotent (twice, and via a fresh poser) | bit-identical | bit-identical |
| frame guards fire | 7/7 refused | 7/7 refused |
| coverage | 3,995 rigid + skin + 3 layers + 1 listed | the same |
| mirrored pairs (new) | 0 label / 0 motion violations; control 11 / 15 | 0 / 0; control 26 / 18 |
| variant at new-coordinate zero reproduces the base poser | — | **3,995 entities × 150 gait frames, worst 6.7e-16 m, rotation 0.0** |
| gait-best posed, ms per pose | 0.92 | 1.07 |

Coordinates in the variant move exactly the entities downstream of them:
`thoracic_extension` 2,069 (thorax + cervical + head), `neck_extension` 1,047, `head_*` 885,
`subtalar_angle_r` 172, `subtalar_angle_l` 171, `wrist_flex_r` 158, and `wrist_dev_l` 156.

**The mirrored-pair gate FAILED in its first form, on the base plant, and is recorded as such.**
The first version compared every moved entity, by mirrored name, between each left and right
coordinate. It read **36** violations on the base plant (45 on the variant). Every one was an
entity with no unique twin in the atlas: unsided names that exist once and sit on one side, such
as `ulnopisiform ligament` and `set of plantar digital arteries proper`, and one-sided vessels
such as `right anterior tibial vein`. No binding can mirror an entity whose twin does not exist,
so v1 counted the wrong population. The gate is now restricted to the 2,684 names whose twin
exists exactly once. The 51 names outside that population are listed in the report and not gated.
It also has a control that must fire, the same counts with `symmetric=False`, and it does: 11/15
on the base plant and 26/18 on the variant. The bar did not move. The population changed, and v1's
count is still printed on every run.

The base battery's ten original gates print the same numbers as before this change. Separately,
`OsimKinematics` had built its topological order with every joint TWICE (41 entries for 22
joints; see the comment in the code). That is fixed, and the base poser was checked against the
committed module on all 150 gait-best frames in both pivot modes. The rotation, translation,
centroid and segment-motion arrays and the skin were bit-identical.

### The live plant runs the variant (18 September 2026)

`ArticulatedBodyPlant(root, out, environment='upright', target_mass_kg=77.6122029,
augmented_registration='data/models/articulated_spine_v1/registration.json', display_pose='opensim')`
now constructs, steps, and projects frames with 25 segments. Its display pose moves the skull
with `head`. The runtime had been built layer by layer for exactly 22 bodies, and two layers
refused the variant, in this order:

| layer | refusal | now uses | commit |
|---|---|---|---|
| force-frame registration (`CanonicalRegistration`) | 'Native segment has no retained canonical bone anchors: thorax' | bodies `canonical/mechanics.json` does not register are anchored on this binding's `segment_named_bones` (thorax 43, head 89, cervical 10), and those bones leave `torso` | `7ded91b` |
| display poser | would have refused the 25 bodies in `check_bodies` | `AnatomyPoser.from_workspace(root, plant)` for the identified plant | `7ded91b` |
| skin blend (`ContinuousSurfaceBinding.from_root`) | 'Surface binding registration supports changed' | a second blend over 25 segments, selected by `articulated.SURFACE_BINDINGS[plant]` | this change |

After the third fix, nothing else in that construction refused. The variant also runs with
`display_pose='anatomical'` and with `environment='supine'`. `enable_garments=True` fails on the
BASE plant as well ('Garment constructor hash mismatch'), so that failure is not the variant's.
The plant does not claim the variant stands. Why its ankle collapses is being diagnosed
separately.

The base plant is bit-identical to `7ded91b`: 7 frames and 6 output files, 0 differences, and a
1-ulp control that fires. The display-pose battery is unchanged. Everything is in
`scripts/verify_variant_runtime.py`, and its receipt is
`data/derived/variant-runtime-verification.json`.

### What still does not follow

* **The skin: done for the plant, not yet for the poser.** This bullet used to say that the
  skin's linear blend covered only the 22 base segments, so the variant's head and neck skin
  rode `torso` while the skull nodded underneath it. A 25-segment blend now exists, built by the
  same method (see `DISTRIBUTED_SURFACE_BINDING.md`, "The 25-body variant"). The live plant's
  force-frame skin and `surface_transforms` use it. `AnatomyPoser.skin_vertices` still reads
  `SKIN_BINDING`, which is the 22-segment blend, for every plant. So in the variant the POSED
  skin still follows `torso`. That is a one-line change in `anatomy_pose.py`: select the blend by
  plant, as `PLANTS` does.
* **The app: done.** This bullet used to say that `ArticulatedBodyPlant` called
  `AnatomyPoser.from_workspace` with no plant. Since `7ded91b` it passes the plant it identified
  from the engine's body set.
* **Mass and surface disagree.** 1,128 entities, including the thoracic vertebrae and the lung lobes, are
  bound to `thorax` or `cervical` but have their mass in the torso residual core, because the
  variant's partition is rib cage only.
* **The girdle.** Scapulae, clavicles and the arms ride `torso`, not `thorax`, as the model
  does. When the thorax rotates, the ribs slide under the scapulae. The worst cross-joint
  opening at `thoracic` (76.2 mm OpenSim, 71.7 mm anatomical pivots) is exactly that scapula–rib
  pair, not the T12/L1 joint.
* **Tearing at the new joints.** 127 re-voted entities have family coherence below 0.6: neck
  fascia, platysma, sternocleidomastoid, rhomboid minor, the infrahyoids. These are the
  structures that span the new joints, and a rigid binding tears them.
* **Worst openings at the new joints** (64 closest bone pairs, through each range):
  * neck: 38.4 mm with OpenSim pivots, 14.2 mm with anatomical pivots;
  * atlantooccipital: 18.5 → 4.5 mm;
  * subtalar: 17.1 → 15.1 mm;
  * wrist: 11.5 → 11.8 mm.
