# nociception: the receptor pain-guided RL was blocked on

`docs/ACTUATION_STAGES.md` records that pain-guided reinforcement learning was blocked on
a receptor, not a pathway: the fibre classes and the routes existed, and nothing
transduced damage at the skin. This is that receptor, in both halves, and a precise
account of what it is and is not.

**What it delivers is transduction: stimulus in, afferent firing out, with provenance.**
It is not a pain scale and it is not a reward. Turning these rates into a signal to
optimise is a separate decision that has not been made. The human recordings the
receptor is built from say plainly why the two must not be confused: Van Hees and Gybels
(1981) recorded mechanically evoked C-nociceptor discharge "up to more than 10 spikes/s"
that was "not necessarily accompanied by pain sensation", and Adriaensen et al. (1984)
titled their paper *a paradox* because pain rose through a 120 s squeeze while the
polymodal C discharge adapted away.

| | where | what |
|---|---|---|
| cutaneous transduction | `ihm/assembly/nociception.py` | contact pressure → A-delta and C rate at each of the 1,326 patches |
| exact stimulus binding | `scripts/build_nociceptor_patch_binding.py` | skin triangle → patch, recovered exactly |
| gate | `scripts/gate_nociception.py` → `out/gate_nociception.json` | 18 checks, plus one that all 18 ran |
| visceral component | IBM-1 `ibm/fields/transduction.py`, `ibm/interoception.py` | `transduction.visceral_nociceptor`, with the splanchnic rows bound to it |
| routes | IBM-1 `ibm/topologies/ihm_bridge.py::cutaneous_nociceptive_routes` | A-delta and C delay per patch |
| gate | IBM-1 `scripts/gate_nociception_routes.py` → `out/gate_nociception_routes.json` | 9 checks |

---

## 1. The cutaneous receptor

### The law, and whose numbers it uses

```
x    = max(0, P - P_threshold) / (P_half - P_threshold)
rate = ceiling * x / (1 + x)                     (Hz; P in Pa, compressive positive)
```

| constant | value | derivation | source | kind |
|---|---|---|---|---|
| `P_threshold` | **133.3 kPa** | 4 N / 30 mm² | Adriaensen H, Gybels J, Handwerker HO, Van Hees J (1984). *Nociceptor discharges and sensations due to prolonged noxious mechanical stimulation — a paradox.* Hum Neurobiol 3:53–58. PMID 6330012. "forces greater than 4 N exerted on forceps faces of 30 mm2 elicited pain" | human pain-report threshold |
| `P_half` | **466.7 kPa** | 14 N / 30 mm² | Schmidt R, Schmelz M, Torebjörk HE, Handwerker HO (2000). *Mechano-insensitive nociceptors encode pain evoked by tonic pressure to human skin.* Neuroscience 98:793–800. doi:10.1016/s0306-4522(00)00189-5. "tonic pressure stimulation (14N at 30 mm(2); 120 s)" evoked augmenting pain and recruited silent C units | the one human suprathreshold pressure in the same geometry; half-saturation is **placed** here |
| `ceiling` | **10 Hz** | — | Van Hees J, Gybels J (1981). *C nociceptor activity in human nerve during painful and non painful skin stimulation.* J Neurol Neurosurg Psychiatry 44:600–607. doi:10.1136/jnnp.44.7.600. "Mechanically evoked C fibre discharge even up to more than 10 spikes/s" | an observed maximum, and "more than" makes it a **floor** on what was seen |
| the curve | rectified hyperbola | — | none | **modelling choice**: the sources give a threshold, a strongly painful level and a ceiling, not a curve |

**Why these two anchors and not the more familiar algometry.** Threshold and
half-saturation are both forceps squeezes on 30 mm² faces, so both are the same division
of force by the same area. No constant crosses a change of units or of probe geometry —
the trap this programme's ledger records as "a constant is not portable across a change
of units".

What the law gives, for orientation:

| pressure (kPa) | 100 | 133.3 | 165.9 | 200 | 270.9 | 300 | 466.7 | 1,000 | 2,000 | 5,000 |
|---|---|---|---|---|---|---|---|---|---|---|
| rate (Hz) | 0 | 0 | 0.89 | 1.67 | 2.92 | 3.33 | 5.00 | 7.22 | 8.49 | 9.36 |

### An independent human check, which disagrees by up to 2×

Algometer pressure pain thresholds over a 1 cm² probe (DFNS protocol: Rolke R et al.
2006, Pain 123:231–243, "probe area of 1 cm2 … ramp of 50 kPa/s"), in 130 healthy adults
(Pan LH et al. 2024, *The normative values of pain thresholds in healthy Taiwanese*,
Brain Behav 14:e3485, doi:10.1002/brb3.3485):

| site | men (n=55) | women (n=75) | pooled by n | ÷ model threshold |
|---|---|---|---|---|
| thenar eminence | 290.3 ± 91.3 kPa | 256.7 ± 94.8 | **270.9** | **2.03×** |
| masseter | 178.5 ± 56.7 | 156.6 ± 58.4 | **165.9** | **1.24×** |

These are not used in the law. The pooled figures are derived in code
(`nociception.pooled_ppt_kpa`), not typed. The direction of the disagreement matters:
spatial summation (Lautenbacher S et al. 2005, Pain 115:410–418) predicts the **larger**
probe gives the **lower** threshold, and the 1 cm² algometer reads **higher** than the
30 mm² forceps. So the difference is site and loading mode, not area. The honest statement
is that human mechanical pain thresholds span roughly **130–290 kPa** across these two
instruments and sites, and the model sits at the low end of that span. It applies one
threshold to all twelve patch regions. It was measured on dorsal hand skin, and no source
read here measured the others. DFNS trunk reference data (Pfau et al. 2014, Pain
155:1002–1015) would begin to fill that gap.

### What it does not model, each for a stated reason

1. **Sub-pain nociceptor firing.** The threshold is where humans *report* pain. Human
   polymodal C units have von Frey thresholds of 2.3–13.1 g (Van Hees and Gybels 1981)
   and < 160 mN (Schmidt et al. 2000), well below it. Those are forces on hair tips whose
   area neither abstract gives, and turning them into a pressure would mean inventing
   one. They are recorded in `EVIDENCE['sub_pain_von_frey']` as **NOT USED**. This
   channel is therefore silent over a range where real C-nociceptors are not. It
   understates sub-pain activity by construction.
2. **Adaptation and sensitisation.** Both are measured. Polymodal C units show "an initial
   high frequency dynamic discharge followed by adaptation" (Adriaensen 1984), and silent
   units were activated only "after more than 20s" and responded more strongly to a
   second stimulus (Schmidt 2000). Neither source gives a time constant. The transduction
   is static.
3. **A separate A-delta law.** No human A-delta mechanical stimulus-response function was
   read. Adriaensen H et al. (1983), J Neurophysiol 49:111–122,
   doi:10.1152/jn.1983.49.1.111, is the source, and its text was not available. The
   A-delta channel therefore uses the C law unchanged. **That is an unmeasured
   assumption** (`EVIDENCE['adelta_law']`), and that paper's high-threshold
   mechanoreceptor data would measure it. What does differ between the classes is
   conduction velocity, which is why they are routed separately.

---

## 2. Binding the stimulus to the patches

**The crude body's contact elements are not converted.** Each one (22 in the
supine-support run's `initial_native.json`) is a sphere with a force and no contact area, so a pressure cannot be computed without inventing an area.
Leaning on the scaffold for this is what ACTUATION_STAGES.md says not to do.

**The supine surface quadrature is bound exactly.** Each of its 21,382 samples records the
skin triangle it lies on. The 1,326 patches record a centroid, an area and a triangle
*count*, but not which triangles. Assigning samples to the nearest centroid is the coarse
proxy CLAUDE.md records failing three times on a mesh with a long facet tail, so the
membership was recovered instead. `scripts/build_dermatome_patches.py` is deterministic.
`build_nociceptor_patch_binding.py` re-runs it with the output redirected to a temporary
directory, records each cluster, and **refuses unless the rebuilt `dermatomes.json` is
byte-identical to the canonical one**. It was byte-identical: 109,183 exterior faces
bound, all 1,326 triangle counts equal, worst relative area mismatch 0.0. The canonical
file's mtime is unchanged. The binding is keyed by the sha256 of `dermatomes.json` and of
the skin geometry, and the loader refuses a mismatch.

**Pooling.** The rate is computed at each sample and then pooled per patch as a
skin-area-weighted mean. Unloaded skin counts as zero, so the result is the mean firing of
endings spread evenly over the patch. The peak sample rate is reported beside it. A
per-patch force (`patch_rates_from_patch_forces`) gives the patch-mean pressure, which is
a **lower bound** on the peak whenever the load is concentrated.

### Three corrections made during this work, recorded as they happened

**Run 1: the pooling mixed two measures and the known answer could not see it.** The
quadrature's `area_m2` is **projected** area (a 5 mm raster normal to the bed). The patch
area is **surface** area. Pooling one against the other put a uniform load on one patch
at **7.0985 Hz when every sample on it read 7.0588 Hz**, a mean above its own maximum.
Check B3 reported PASS, because its "expected" value used the same formula as the code.
It tested the arithmetic against itself. Record: `out/gate_nociception_run1_b3_tautology.json`.
Fix: B3m was added. It loads every sample and requires that no patch mean exceed the
sample rate. A mean above its maximum is impossible whatever the formula, so this check
can fail.

**Run 2: the first fix discarded skin.** Each sample's skin area is its projected area
divided by |cos θ| between the face normal and the ray. The ray is carried from the source
frame by the manifest's rigid `source_to_canonical_ground`, and **every sample, carried
back, lies on its own face to 4.5e-16 m**, so the map is checked rather than trusted. As a
guard against grazing cells, the samples on each face were then capped at that face's
area. That cut the posterior skin from **0.868 to 0.489 m², below the 0.535 m² projected
area it came from**, which a projection can never exceed. A 5 mm cell covers 25 mm², and 81%
of exterior skin faces are smaller (median face 7.7 mm², median edge 4.4 mm), so its
centre ray hits one face while its footprint spans many. Record: `out/gate_nociception_run2_face_cap.json`. Fix: no per-face cap. A
patch whose samples sum past its own area is normalised by that sum and **counted**. S2
was added: total skin area must be ≥ total projected area.

**Runs 2 and 3: an idempotence check was deleted by the edit that fixed run 1, and the
gate still said PASSED.** Rewriting the B3 block also removed I2 (pooling called twice),
and both runs reported "16 of 16" without noticing that the docstring declared 17.
Record: `out/gate_nociception_run3_missing_i2.json`. I2 is restored, I2s is added (the
skin-area conversion called twice), and a final check, ALL, now parses the declared
checks from the docstring and fails if any did not run. It was made to fail on purpose:
with I2 suppressed it prints `ALL FAIL … missing ['I2']`.

**What the final pooling still carries.** Under a uniform load on every sample, 717
patches are loaded and **234 of them are over-covered**. The median excess is 4.8%, the
90th percentile 23%, and the worst ×2.88. Samples at grazing incidence (|cos| < 0.1)
account for only 32 of the 234, so most of it is 5 mm raster excess at patch boundaries.
The quadrature's own manifest lists "boundary-cell area … unresolved". The weights still
sum to ≤ 1, so no patch rate can exceed its sample maximum (asserted in code).

---

## 3. Routing: fast and slow on the right fibre

Each patch's rate is delayed by `path_length_m / v`, separately for A-delta and C. The
typical velocities are 15 and 1.0 m/s, taken from IBM-1's `nerve.FIBRE_VELOCITY_M_S`
through IHM's `peripheral.json` snapshot, and IBM's gate C4 asserts the two still agree.
The human C conduction velocity measured by Van Hees and Gybels, 0.86 ± 0.17 m/s, sits
0.8 SD below the 1.0 used. IBM's join (`cutaneous_nociceptive_routes`) **raises** if any
patch's trunk does not carry both A-delta and C. All 21 trunks the patches use carry
both, and C2 confirms the raise works by removing `c` from one of them.

| class | delay min | median | max |
|---|---|---|---|
| A-delta | 4.0 ms | 30.9 ms | 80.2 ms |
| C | 60.0 ms | 463.3 ms | 1,202.3 ms |

**Every one of these is a schematic lower bound, not a measured conduction time.** Each
route is an authored polyline from the patch centroid, through a waypoint, to a relay. It
is not a dissected nerve. IBM-1's `ihm_bridge.py` carries the standing caveat, raised
when the spinal relays sat 89–337 mm below their segments, that spinal-route delays are a
lower bound. The relays have since moved onto the dura centreline and the caveat has not
been lifted, so it is carried here unchanged. Central and synaptic delays are excluded.

---

## 4. The visceral half (IBM-1)

`transduction.visceral_nociceptor` now exists. It sits on the `viscera` support, is
written by the `transduction` process and read by the afferent relay, so the registry
does not report it dead (a probe component with neither is reported dead, which shows the
check is live). The five splanchnic rows are bound to it instead of
`transduction.baroreceptor`, following Cervero F and Jänig W (1992), *Visceral
nociceptors: a new world order?*, Trends Neurosci 15:374–378: they describe high-threshold,
silent and intensity-encoding visceral receptors, and the first two belong on this
component. `interoception.check()` now refuses a row whose `nociceptive` flag and receptor
component disagree, in either direction. The gate corrupts one row each way and both
corruptions raise.

**Still open: a visceral transduction law.** IHM emits these channels in mL, mg/dL and
mL/min, and human visceral pain thresholds are published as barostat pressures in mmHg.
Converting one to the other needs an organ compliance the body does not declare, so no
threshold is applied. `ibm.interoception.ONTOLOGY_GAPS` records this.

---

## 5. Known answers, as printed

```
K1  PASS  rates [0.0, 0.0, 0.0, 0.0] Hz below/at threshold       (0, 0.5 P0, P0(1-1e-9), P0)
K2  PASS  rate at P_half = 5.0 Hz (ceiling/2 = 5.0)
K3  PASS  rate at 1e4 P0 = 9.997500 Hz of 10.0
K4  PASS  non-decreasing over 1 kPa..1 GPa
I1  PASS  firing_rate_hz called twice on 4,096 pressures
B1  PASS  1326 of 1326 patches match their recorded triangle count
B2  PASS  21382 of 21382 supine samples bound (supine-surface-contact-5jqy1juo)
S1  PASS  max sample off its face 4.52e-16 m; projected 0.5346 m2 -> skin 0.8681 m2
S2  PASS  skin area >= the projected area it came from
B3  PASS  one patch, uniform 2 P_half: 7.058823529 Hz = sample rate
B3m PASS  uniform 2 P_half everywhere: max patch rate 7.058824 <= 7.058824; 234 of 717 over-covered
I2  PASS  patch_rates_from_quadrature called twice on the same field
I2s PASS  supine_sample_surface_area called twice on the same quadrature
B4  PASS  0.999 P0 on all 21382 samples -> max patch rate 0.0 Hz
R3  PASS  zero stimulus: 0 events sent, no arrivals
I3  PASS  one delay-line state stepped twice with one input
R1  PASS  1326 of 1326 patches: A-delta first; lead min 56.0 ms, median 433.0 ms, max 1122.0 ms
R2  PASS  worst |first arrival - L/v| = 1.000 ms (dt 1 ms)
ALL PASS  every declared check ran: declared 18, ran 18
```

IBM-1 `scripts/gate_nociception_routes.py`: V1–V3 and C1–C6 all PASS. Every patch has
A-delta ahead of C at the typical velocities, and also in the worst case, with the
slowest A-delta (5 m/s) against the fastest C (2 m/s).

---

## 6. A prior nociceptor this does not reconcile

IBM-1 commit `8b85e1a` (9 September) added `nociceptor_polymodal` to
`ibm/processes/transduction.py`, and it covers the same cutaneous ground. It is not what
transduces at the patches, and it disagrees with the human data above on three points:

- **its threshold is a force with no area**, `mechanical_threshold_n` 8 N (provenance
  WEAK). As a pressure, that is 267 kPa over 30 mm², 80 kPa over 1 cm², and about 5 kPa
  over a 15 cm² patch, so its meaning depends entirely on an area it does not declare;
- **its ceiling is 100 Hz**, ten times the maximum mechanically evoked human C discharge
  cited here;
- its docstring describes its rate as "usable as a reward signal", which is the step the
  receptor should not take on its own.

It was not edited, because it is outside this work's territory. Either retire it or
reparameterise it from `nociception.EVIDENCE`, and that is a decision for whoever owns
`ibm/processes/transduction.py`.

---

## 7. What pain-guided RL still lacks after this

- **A reward decision.** Deliberately absent. Nothing here says how A-delta and C rates
  become a cost, and the human data say they are not pain.
- **A dynamic stimulus.** The binding is to the *reference* supine quadrature. There is
  no settled supine load in the repository to feed it, and the quadrature's own limits
  say large rotations need a rebuilt one. The live frame's contact, the crude
  body's spheres, carries no area. Until skin, not a sphere, meets the floor in the running
  body, there is no damage stimulus to transduce during motion.
- **Sub-pain firing, adaptation, sensitisation and a real A-delta law** (section 1).
- **Regional thresholds.** One threshold, measured on the hand, is applied everywhere.
- **Thermal and chemical nociception.** Only mechanical pressure is transduced.
- **A visceral law** (section 4).
- **A central pathway.** The rates arrive at a relay. Spinothalamic projection, dorsal
  horn gating and descending modulation are absent. IHM's routes end at a relay and
  exclude central delays.
- **Measured routes.** Every delay is schematic.

## Rerun

```
.venv/bin/python scripts/build_nociceptor_patch_binding.py   # once per dermatomes.json
.venv/bin/python scripts/gate_nociception.py
cd ../IBM-1 && PYTHONPATH=. .venv/bin/python scripts/gate_nociception_routes.py
```
