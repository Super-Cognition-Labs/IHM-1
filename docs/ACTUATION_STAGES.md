# how the body comes to move: the staged plan, and what is scaffold

Read this before working on anything that moves the body. It exists because the
low-resolution body is easy to mistake for a component of the design when it is
a **scaffold with a planned disposal**, and because this programme's recorded
failure mode is drifting toward whatever is measurable on the proxy.

The end state, stated first so nothing below is mistaken for it:

> **The brain directly controls individual muscles. Those muscles' contractions
> control the skeleton. The whole thing is a highly complex rigid AND soft body
> system, actuated by the dynamic elasticity of the muscles, which is determined
> by the brain's outputs.**

**Contact with the world is never bone against world.** It is bone, mediated
through the interlayering of soft bodies — fatty tissue, muscle, skin and the
rest — and it is the skin that meets the floor. A bone mesh is a collider against
*other bones and its own soft tissue*, not against the ground.

**And the end state is not a later stage. It is to be fully implemented now.**
The body must be selectable between two modes:

| mode | meaning |
|---|---|
| **driven** | the body is posed from the crude scaffold. Fast, and what generates training corpora. |
| **fully present participant** | the body is dynamically simulated in its own right: muscles actuate, soft tissue deforms and mediates contact, the skeleton moves because muscles pulled it. |

Both exist, and the caller chooses. The stages below describe how work is
sequenced and what the scaffold is *for* — they do not license shipping only the
scaffold and calling the participant mode future work.

---

## The two bodies, and which one is real

| | crude body | real body |
|---|---|---|
| what | 22 rigid bodies, 80 muscles, 33 coordinates | 3,816 anatomical entities |
| driven how | dynamically integrated by the engine | **kinematically posed** from the crude body through the segment binding |
| touches the world | yes — all contact happens here | no |
| status | **scaffold** | the thing being built |

The binding is `data/derived/anatomy-segment-binding/` (99.88% of entities
assigned, gates in `report.json`). `ihm/assembly/body.py` states the real body's
transforms come from the run and never hold a last value.

**The crude body has always been a scaffold.** It is not a simplification we
have settled for; it is a stage. Work that makes the scaffold better *as a
scaffold* is worth doing. Work that deepens our dependence on it is not.

---

## The stages

### 1. Now — the scaffold actuates, and generates training data

The crude body actuates the real body's muscles and bones. It determines rough
contact with the floor and with objects. That motion, and the afference it
produces, is **collected as training data for teaching the brain to actuate the
real body**.

This is why forced motion is legitimate here and only here: a pose trajectory
imposed on the body produces what real muscles and skin experience, and that
corpus bootstraps circuits that would otherwise have nothing to condition on.
`scripts/collect_pose_corpus.py` (68 motions, 17,622 frames) and
`scripts/collect_forced_gait_corpus.py` are that step.

Contact is rough on purpose at this stage. It only has to be good enough to
generate plausible load and afference.

### 2. Possibly — the brain learns to actuate the scaffold too

Letting the brain drive the low-res body may earn its place as a curriculum
stage: a lower-dimensional control problem to solve before the real one. Optional,
and justified only by whether it helps the transfer.

### 3. The scaffold's jobs move to the real body

Not because the participant mode arrives later — it is to be built now — but
because each job the scaffold still holds is a place the real body is not yet
sovereign. Contact first, then actuation. When a job has moved, the scaffold
keeps it only as the `driven` mode's fast path.

### 4. The end state

Brain output → individual muscle activation → muscle contraction → skeleton
motion, in a coupled rigid and soft body system where the muscles' dynamic
elasticity is what does the work.

---

## How the brain learns to drive it

The point of all of the above, stated so the body work is not mistaken for an end
in itself: **watching the IBM-1 brain control a real naked human body as it learns
to pick itself up.** "Pick itself up" is the first concrete milestone — it needs
no locomotion, it is reachable from prone, and it is unambiguous on video.

The learning shifts in two phases:

**Supervised, from the scaffold.** The low-resolution OpenSim model supplies pose
trajectories; those become forced motion on the real body; the afference they
produce is the training signal. This is what `scripts/collect_pose_corpus.py`
(68 motions, 17,622 frames) exists for, and it is why the scaffold earns its
keep.

**Then pain-guided reinforcement learning and self-supervised prediction, on the
real body.** Control is learned from consequence rather than demonstration.

### What pain-guided RL requires, after the receptor was built (18 Sep 2026)

**The receptor exists now.** This section used to say there was no nociceptor
transduction component and that pain-guided RL was blocked on a receptor, not a
pathway. Both halves are built and gated — `docs/NOCICEPTION.md`:

- **Cutaneous** (`ihm/assembly/nociception.py`, `scripts/gate_nociception.py`,
  19 of 19): contact pressure in, Aδ and C firing out, at all 1,326 innervated
  patches. Threshold 133.3 kPa (4 N over 30 mm², Adriaensen et al. 1984),
  half-saturation 466.7 kPa (14 N over the same area, Schmidt et al. 2000), ceiling
  10 Hz (Van Hees & Gybels 1981, the highest C rate seen under mechanical
  stimulation, and a lower bound on the true maximum). Both anchors share one
  contact area, so no number crosses a change of units. Aδ arrives first on every
  patch, by 56–1,122 ms.
- **Visceral** (IBM-1 `ibm/interoception.py`): `transduction.visceral_nociceptor`
  exists on the viscera support, and the five splanchnic high-threshold channels
  are bound to it instead of to `transduction.baroreceptor` with a nociceptive tag.
  `interoception.check()` refuses any channel whose flag and receptor disagree.

**What it still lacks, so the gap is not re-closed by a summary line:**

- **a reward decision** — deliberately absent. Transduction is not a reward, and
  choosing one is a separate decision;
- **a damage stimulus during motion** — contact is bound to the supine reference
  pose only, and the scaffold's contact spheres carry a force but no area, so they
  cannot be converted to pressure;
- **a visceral firing law** — the channels are bound but IHM emits them in mL and
  human visceral pain thresholds are published as pressures; the conversion needs
  an organ compliance the body does not declare;
- sub-threshold firing, adaptation, sensitisation, a separately measured Aδ law,
  per-region thresholds, and heat, cold and chemical nociception;
- any central pain pathway past the relay, and measured nerve routes — every
  delay is a schematic lower bound.

**A conflict the owner of IBM-1 has to resolve.** IBM-1 already carried
`nociceptor_polymodal` (`ibm/processes/transduction.py`, 9 Sep) with an 8 N
*force* threshold and no contact area, a 100 Hz ceiling — ten times the human C
maximum above — and a docstring describing its rate as usable as a reward signal.
Two nociceptors now disagree by an order of magnitude at the ceiling and in their
units at the threshold. It has not been touched here: it is another lane's
component, and re-basing it changes that lane's results.

---

## What this means for decisions you are about to make

**The participant mode needs the real body to be dynamic, and that is the point.**
Earlier guidance here said not to give the anatomical body its own dynamics; that
was written when the end state was mistaken for a later stage, and it is
withdrawn. What remains true is that the `driven` mode must stay cheap — it is
the corpus generator — so the two modes are different code paths over the same
anatomy, not one compromise between them.

**Contact belongs to skin, not bone.** A design that puts the ground collider on
a bone mesh has skipped the fat, muscle and skin that actually meet the floor.
The 21,381-point skin contact quadrature and `surface_contact_manifest` already
exist and are the right place to start.

**Do not optimise the scaffold's fidelity as an end.** Making the 22-segment
body's gait beautiful is not progress toward a brain controlling real muscles.
Making its contact good enough to generate honest afference is.

**Do not report scaffold behaviour as body behaviour.** A trajectory the crude
body executed is a trajectory the crude body executed. The real body was posed
from it. Both are true; only one is what the programme claims to be building.

**Forced motion is never a result.** It is stage-1 data generation. This is
recorded in `docs/DIRECTION.md` in the IBM-1 repo and repeated here because it is
the single easiest thing to misreport.

---

## Where the pieces are

- segment binding — `data/derived/anatomy-segment-binding/`, `docs/ANATOMY_SEGMENT_BINDING.md`
- pose corpora — `scripts/collect_pose_corpus.py`, `scripts/collect_forced_gait_corpus.py`
- the plant and its contact geometry — `scripts/native_mechanical_stream.cpp`
- tissue as force elements, and the three different reasons the rest is blocked —
  `docs/TISSUE_MECHANICS.md`, `scripts/build_tissue_force_elements.py`
- innervation, so the brain reaches muscle and skin only through nerves —
  `scripts/measure_innervation_coverage.py` in IBM-1
- the standing direction and its corrections — `docs/DIRECTION.md` in IBM-1
- where declared and running models differ — `docs/DISCONNECTS.md` in IBM-1
