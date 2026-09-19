# the soft body: a deformable tissue layer over the rigid skeleton

`docs/ACTUATION_STAGES.md` defines two modes, `driven` and **fully present participant**,
and says the second is to exist now. Before this file, nothing in the stack could deform:
every compliant element was a 1-D law on a rigid carrier. The supine foundation is one
independent confined column under each quadrature point, and the segment contact meshes are
an elastic foundation with one spring per triangle. Neither can bulge, share load with a
neighbour, or tell a heel-sized footprint from a fingertip-sized one.

This is the first piece of the participant mode: **a 3-D soft-tissue layer per segment that
deforms under contact and returns the load it transmits to its bone.** It is real, verified
against analytic answers, and wired so a live body can ask for it. It was first built slow,
not mesh-converged at the contact surface, and apparently softer than every measurement it
was built from. **The 2026-09-18 section below revises all three.** The "softer"
claim was measured on a seam wedge, not the heel, and is withdrawn. The surface is now
body-fitted and its heel force moves monotonically with refinement, but is not converged.
The solve is 10–20x faster, and still 9–16x short of the plant step. It was also not in
the plant's integration loop; **the first section below puts it there, and measures what
that costs.** It is still not the whole participant mode: it is the part that was buildable
from what was already on disk, with the reason each remaining part is blocked.

**And the coupling is not how the layer should reach the plant.** The first section below
measures that an explicit reaction is 10x outside its own accuracy bar even at one solve per
step; the section above it takes the consequence — the layer is used OFFLINE to calibrate the
engine's own implicit contact element, which is solved inside the integrator at no per-step
cost. That works on the heel to 4.8–6.1% and fails on the forearm at 39.8%, and both halves
are results.

Code: `ihm/assembly/soft_tissue_layer.py` (the layer) and `ihm/assembly/contact_law.py` (the
fitted law). Selection: `ihm/assembly/plant_options.py` — `soft_tissue`: `layer_map_*` is the
layer as first built, bit for bit; `layer_fitted_local_*` is the fitted, local-depth layer;
`*_coupled` is the same layer, in the plant's loop; and `segment_contact: skin_layer_fitted`
is the implicit law the layer implies.
Batteries: `scripts/verify_soft_tissue.py` → `data/derived/soft-tissue-layer-v1/report.json`,
`scripts/verify_soft_tissue_coupled.py` → `data/derived/soft-tissue-coupled-v1/report.json`,
and `scripts/verify_contact_law_fit.py` → `data/derived/contact-law-fit-v1/verification.json`.
Also: `scripts/build_soft_tissue_local_depth.py`, `scripts/census_soft_tissue_fitted.py`,
`scripts/measure_soft_tissue_speed.py`, `scripts/measure_contact_law_fit.py`.

---

## 2026-09-18 (later still): the contact law the layer IMPLIES, fitted offline and solved implicitly

The section below this one ends on a negative that is really an instruction: the reaction
grows at **1,472.86 N/s**, so a 10% bar allows a hold of **1.0 ms**, and *even one solve per
plant step is 10x outside that bar* — an explicit coupling cannot be repaired by running it
more often, it has to be **implicit**. It also costs 53.9 ms per 10 ms step on one segment.

The engine already carries an implicit contact element: `ElasticFoundationForce` over a
`ContactMesh`, solved inside the error-controlled integrator, which is what the segment
contact bundles use. So the layer is used **offline as the ground truth** and that element's
one elastic parameter is fitted to it. The answer is not uniform, and the shape of the split
is the result:

| | |
|---|---|
| **on the heel it works, on both heels** | a single stiffness reproduces the 3-D layer to **4.84%** (`calcn_l`) and **6.11%** (`calcn_r`) median on penetrations the fit never saw, inside the layer's own 10% convergence bar |
| **and the shipped rule was already right there** | the bundle's own **unfitted** `k = E_app/h` reads **4.47%** on `calcn_l` — very slightly BETTER than the fit — and 16.00% on `calcn_r`. F2 FAILED on both: neither fit beats the trivial baseline by the √2·bar a measured baseline demands |
| **on the forearm it does not work** | on `ulna_l`'s own collapse poses the fit is **39.84%** median and **470.71%** at worst, and `F_layer/G` spans **12.0x**. F1 FAILED, and by this measurement's own pre-registered outcome (c) `ulna_l` carries **no fitted stiffness at all** in the delivered bundle |
| **and a single k is pose-dependent even within one segment** | the `calcn_l` stiffness fitted at the 17.26° heel pose is **22.76%** off at the 12.59° one. The two poses' effective stiffnesses are 1.031e7 and 1.352e7 Pa/m — **31% apart** |

Code: `ihm/assembly/contact_law.py`. Measurement: `scripts/measure_contact_law_fit.py` →
`data/derived/contact-law-fit-v1/report.json`. Battery:
`scripts/verify_contact_law_fit.py` → `.../verification.json`, 22 gates, 1 FAILED and
recorded. Selection: `plant_options.SEGMENT_CONTACT_BUNDLES['skin_layer_fitted']` and
`['skin_layer_fitted_carried']`. Logs: `logs/measure_contact_law_fit.run{1-F2H-F1U-FAILED,2}.log`,
`logs/verify_contact_law_fit.run{1-V1c-FAILED,2}.log`.

### The law was read out of the engine's source, and the engine confirms it to 2.7e-15

`data/raw/mechanics/simbody/Simbody/src/ElasticFoundationForce.cpp`: `springPosition[i]` is
the face **centroid**, `springArea[i]` its area, and `processContact` skips any face whose
centroid is not inside the other object and otherwise applies
`f = k · area · distance · (1 + c·v_normal)` at the **nearest point on the support**.
`areaScale` is 1 because OpenSim registers parameters only for a `ContactMesh`
(`ElasticFoundationForce.cpp:86-92`) and the floor is a `ContactHalfSpace`. Against an
axis-aligned half-space the nearest point is the perpendicular projection, so

    F(pose, k) = k · Σ_{centroid inside} area_f · depth_f · n̂

— **linear in k at a fixed pose**, which makes every fit here a one-parameter closed-form
fit against a purely geometric integral `G(pose)`. No seed, no iteration, no generator.

**The cross-implementation gate (V2).** A plant is built from the fitted bundle and read
**at rest** — the run's own reported speeds are checked to be exactly zero rather than
inferred from the fact that nothing has been advanced — because at rest `v_normal` and the
slip velocity are identically zero and only the elastic term survives. Over ten contacting
skin meshes at three heights and 72 to 5,519 faces, the Python law and OpenSim's own
reported force agree to **2.72e-15** relative, and the engine reports back exactly the
stiffness the artefact fitted. The control that can fail, a 1.5x stiffness, differs by
0.5000.

**In MOTION the same comparison is out by up to 3.2x, and that is not an error.** The
engine's force there also carries friction and the dissipation factor. A quasi-static layer
contains no information about either, so **only the elastic term is fitted** and the
dissipation and friction coefficients stay the shipped engineering constants. That is a
limit of what the layer can say, not a choice.

### The fixtures, the split and the baselines

Every fixture is one this repo had already chosen by measurement.

| | fixture | penetrations |
|---|---|---|
| **H** | `calcn_l` **heel fixture** — least rotation about z keeping the forefoot half clear at 6 mm, by bisection on the skin mesh: **17.2552236013°** | 0.5 → 4.0 mm by 0.5 |
| **H2** | the same derivation at **3 mm** forefoot clearance: **12.5920935886°**. Cross-fixture only, never trained on | 0.5 → 4.0 mm |
| **HR** | `calcn_r` heel fixture: **17.8740180872°** | 0.5 → 4.0 mm |
| **S** | `calcn_l` at identity rotation — **the seam wedge**, see the correction below | 0.5 → 2.5 mm |
| **U** | `ulna_l` at **the plant's own poses**: a 140-step uncoupled collapse from `pelvis_ty = 1.03`, every step whose skin is below the floor, solved COLD at each so the measurement is a function of its own pose | steps 81–103, 1.65 → 10.37 mm |

TRAIN on the odd-numbered penetrations, SCORE on the even ones — interleaved, so the
held-out points are interpolation; `U` trains on even-indexed loaded poses and scores on the
odd ones. A pose where the layer carries exactly zero, or where the foundation touches no
face, is dropped and reported rather than weighted as zero. Baselines on every scored set:
predict zero, predict the mean, **the bundle's own unfitted `k = E_app/h`** (the trivial one
the brief asks for), and the uniform `skin_material` triple E = 3000 Pa, ν = 0.45,
h = 6.6 mm → k = 1.7241e6 Pa/m.

Bars fixed before any fit existed: **F1**, held-out median relative error ≤ **0.10**, which
is the layer's own convergence bar (CV7–CV10) and the only place the number comes from;
**F2**, beat the unfitted baseline by more than **√2 · 0.10 = 0.1414**, because that baseline
was measured against the same uncertain layer and a bare "beats the baseline" compares one
uncertain number with another.

### What each segment came out at

| | `calcn_l` (H) | `calcn_r` (HR) | `ulna_l` (U) |
|---|---:|---:|---:|
| fitted k (Pa/m) | **1.0517e7** | **9.5056e6** | 4.2623e5 |
| as a Young modulus at the declared ν, h | 51.60 kPa, h 18.61 mm | 46.43 kPa, h 18.53 mm | 1.34 kPa, h 11.92 mm |
| the bundle's unfitted k = E_app/h | 1.0347e7 | 1.0392e7 | 3.2826e6 |
| fitted / declared | **1.016** | **0.915** | **0.130** |
| held-out median, **fit** | **4.84%** | **6.11%** | 39.84% |
| held-out max, fit | 8.84% | 12.10% | 470.71% |
| held-out median, **unfitted k** | 4.47% | 16.00% | 374.93% |
| held-out median, uniform E/ν/h | 82.81% | 80.75% | 149.45% |
| skill vs predict-zero / predict-mean (fit) | 0.9985 / 0.9965 | 0.9980 / 0.9953 | 0.8478 / 0.7004 |
| `F_layer/G` spread over the fixture | 1.37x | 1.45x | **12.0x** |
| F1 / F2 | PASS / **FAILED** | PASS / **FAILED** | **FAILED** / PASS |
| pre-registered outcome | (b) | (b) | **(c)** |

Two heels, two segments, and they agree: the elastic foundation **can** carry the layer's
heel response with one number, and the layer's own two heels differ from each other by 10%
in fitted stiffness while the declared rule says they are the same to 0.4%. Both differences
are inside the layer's 10% bar, which is the honest reading of all of it.

**The uniform `skin_material` triple is 81–83% wrong on both heels** and is the baseline the
per-segment layer map replaced; this measurement is the first independent confirmation that
replacing it was right.

### Where it fails, in newtons, and what the failure is NOT

On `ulna_l` the effective stiffness `F_layer/G` reads, pose by pose:

| step | penetration | layer | `G` | `F/G` | layer contact nodes |
|---:|---:|---:|---:|---:|---:|
| 81 | 1.654 mm | 0.0142 N | 1.844e-07 | **7.69e4** | **4** |
| 82 | 1.732 mm | 0.0151 N | 2.017e-07 | **7.47e4** | **4** |
| 98 | 1.999 mm | 0.4472 N | 4.991e-07 | 8.96e5 | 19 |
| 99 | 3.921 mm | 1.6869 N | 2.359e-06 | 7.15e5 | 67 |
| 100 | 5.741 mm | 4.0092 N | 5.818e-06 | 6.89e5 | 121 |
| 101 | 7.436 mm | 7.3764 N | 1.084e-05 | 6.81e5 | 182 |
| 102 | 8.984 mm | 11.7449 N | 1.689e-05 | 6.95e5 | 231 |
| 103 | 10.371 mm | 16.9001 N | 2.407e-05 | 7.02e5 | 268 |
| 104 | 11.604 mm | — | — | — | `SoftTissueLeftDomain`, min J 0.135 |

From step 99 on, `F/G` is constant to **±2.5%** — a single k follows the layer there
perfectly well. The 12x is two poses, and they are the two where the layer's own contact
patch is **four nodes**. In absolute terms the whole disagreement at steps 81–82 is
**0.014 N**; the fit's largest held-out error in newtons is **6.64 N** at step 103 against a
16.90 N layer force, and the unfitted stiffness's is **62.11 N**.

**A6, the arm built to test that, and what it actually established.** `ulna_l` was rebuilt
at 4 mm and 3 mm spacing (14,874 → 27,555 → 58,071 DOF) and re-solved at the same poses.
At step 81 `F/G` reads **7.69e4 / 6.30e4 / 6.43e4** — it does not move toward the deep-pose
value. But the contact node count went **4 → 2 → 4**: refining the mesh globally did not
refine the *contact patch*, so the arm **excludes a global mesh-coarseness effect and does
not test the node-count hypothesis at all**. It is recorded as that, neither confirming nor
excluding. The instrument that would settle it is a locally refined patch, and it does not
exist. (The finer meshes also leave the constitutive domain at the deeper poses — min J
0.018 at step 103, 4 mm — so the comparison cannot be carried there.)

**Two mechanisms proposed for the forearm's 7.7x, one measured down and one excluded.**

* *The declared thickness is a segment MEDIAN and the contact is not at the median.* Measured
  from the depth map's own values at the skin vertices under the contact patch: `ulna_l`
  reads **23.71 mm** (range 15.16–26.61) against its declared median of 11.92 mm — **2.0x**.
  `k = E_app/h_local` is then 1.650e6 Pa/m, which is exactly **half** the declared 3.283e6
  and accounts for a factor of 2.0 of the 7.7x. **The remaining 3.9x is not thickness.** The
  same correction makes the heels slightly WORSE: their local depths under the contact are
  15.79 and 16.57 mm against declared 18.61 and 18.53, giving 1.219e7 and 1.162e7 Pa/m
  against fitted 1.052e7 and 0.951e7.
* *The footprint is too narrow for the confined reading.* **Excluded.** The confined reading
  of the layer modulus is a claim about a load much wider than the layer is thick, so the
  contact footprint's minor extent was measured at every pose: `calcn_l` 8.9–41.1 mm over
  h 18.61 (ratio **0.48–2.21**), `calcn_r` 10.9–38.5 over 18.53 (**0.59–2.08**), `ulna_l`
  14.4–43.9 over 11.92 (**1.21–3.69**). The forearm's footprint is *wider* relative to its
  thickness than the heel's and it is 7.7x softer anyway, so the footprint-to-thickness
  reading does not explain it. What is left is the shape of the rigid core — a rounded
  forearm lets the tissue escape around it in a way a broad flat heel does not — and that is
  **a hypothesis with no instrument here**, stated as one.

### A correction, and a fixture that landed on the wrong thing

Fixture **S** was pre-registered as "`calcn_l`, the flat stance fixture (identity rotation)",
citing the coupled battery's P2. **That citation was wrong.** P2's flat stance pose is the
*plant's own* `calcn_l` transform; identity rotation with the plane taken from the skin
mesh's lowest vertex is vertex 780 on the rim of the cut between `calcn_l` and the toes —
**the 1.7 mm seam wedge this file already withdrew every earlier "heel" number for**. Run 1
duly measured the layer carrying exactly **0 N at 0.5, 1.0, 1.5 and 2.0 mm** and 0.0909 N at
2.5 mm, which reproduces the withdrawn fixture's own CV5/CV6 row (0 / 0 / 0 / 0.094 /
0.096 N). Its cross-fixture figure is **948%** on the one loaded point.

It is kept and recorded exactly as measured, and it is **not used to judge the contact law**:
one loaded point on a seam wedge is a statement about the wedge. The cross-fixture test that
does land on the heel was added beside it (H2, above) and reads 22.76%. The failure mode is
worth the line it costs: *a cross-fixture test has to land on the part of the object the
question is about* — the same correction this file made when it derived the heel fixture in
the first place, made again, in the same place, by someone who had read it.

### The per-step cost, and the one number that is NOT the substitution's cost

| | ms per 10 ms plant step |
|---|---:|
| historical plant, inertia-ellipsoid COM spheres (`None`) | **132.9 / 138.2** — two bitwise-identical runs, so the noise floor here is ±4% |
| the layer-map segment-contact bundle, 20 skin meshes, **unfitted** | **459.8** |
| the same bundle, **fitted** | **455.1** |
| the coupled LAYER, on ONE segment, ON TOP of its plant | **53.9** (this file, previous section) |

**`fitted − unfitted` read +29.4 ms in run 1 of the battery and −4.7 ms in run 2. The sign
flips, and that is the finding.** The fitted bundle adds no force element and no per-step
evaluation: its meshes are byte-identical to the source bundle's by sha256 and its record
set is the same, so the two plants differ by **exactly one constant**. A different stiffness
gives a different trajectory and the error controller does different work on it — a
numerator holding work the treatment never did, the same shape this repo has recorded twice
before. **The difference is not attributable and is not quoted as a cost.**

What IS attributable is structural and is the whole point: **the law is solved inside the
integrator, so it costs nothing per step that the segment-contact bundle did not already
cost.** Against the coupling it replaces — 53.9 ms per 10 ms step, for one segment, and a
step behind by construction — that is the difference between 5.4x short of real time and not
short at all. The bundle itself is not free (460 ms against the COM spheres' 133), but that
is the price of contacting 20 real surfaces instead of 20 balls, and it is paid by
`skin_layer_map` already.

### A bug this work uncovered: the per-segment stiffness had never been selectable

`{'segment_contact': 'skin_layer_map'}` **raised** at plant construction:
`Bundle carries a per-segment layer map; a caller E, p or h would override it`.
`plant_options.resolve_fidelity` decided "does this bundle carry a per-record layer map?"
from a top-level manifest key that `skin-layer-map-v1` does not carry, while
`NativeMechanicalStream` decides it from the **records**. So the resolver passed the uniform
triple with every layer-map bundle and every such selection died. **This body's own measured
per-segment stiffness — the bundle's entire reason to exist — had never been reachable
through the resolver.** Both sides now read the same records, which is the only form of the
fix that cannot drift apart again, and it is gated from both sides (V7: a layer-map bundle
resolves without a triple, *and* a bundle without a layer map still receives one — fixing the
first by dropping the material everywhere would have broken the second).

### What the delivered bundle is, and what it is not

`skin_layer_fitted` (feet replaced) and `skin_layer_fitted_carried` (source feet kept) are
`ihm.segment-contact-meshes.v1` bundles whose meshes are byte-identical to
`skin-layer-map-v1` and whose per-segment stiffness is:

* **fitted**, for `calcn_l` and `calcn_r`, with each record carrying its held-out median and
  max error and its ratio to the layer map's own number;
* **refused**, for `ulna_l` — F1 failed, so by the pre-registered outcome (c) the deliverable
  is the measurement and the record keeps the layer map's k and says so in words;
* **unfitted**, for the other 17 records, marked as such, because no layer measurement exists
  for them.

Fitting the rest is compute, not a blocker: a segment's layer must be built (0.3 s for a
heel, 2.7 s for `ulna_l`, 49.8 s for the torso) and solved at 8 poses, and the largest
segments — pelvis 192,870 DOF, torso 183,729 — would not fit the 4 GiB budget at the
spacings that matter.

**Three things it is not**, and they are in the selection's own disclosure:

1. It is a **fit to a layer converged to no better than ~10%**. No force from it is good to
   better than that, and the two heels' 10% disagreement with each other is that bar being
   visible.
2. It **does not deform**. Each spring is independent; there is no lateral bulge and no load
   sharing between neighbours. That is the whole physical difference between a foundation and
   a continuum, and it is exactly what the `ulna_l` failure is made of.
3. Only the **elastic** term is fitted. The dissipation and the friction coefficients are the
   shipped engineering constants and this measurement says nothing about them.

**And one gate is recorded FAILED.** V1c asked that re-ordering and re-winding the mesh's
faces leave the foundation integral **bitwise** unchanged. It does not: 9.9e-23 m³ on
3.55e-08, **2.79e-15 relative** — the summation order of four floats. A bitwise bar on a sum
is a bar on the order of the sum. The bar was not moved; V1d was added beside it at 1e-14 and
passes, and V1c still prints FAILED.

---

## 2026-09-18 (later): the layer is IN the plant's loop, and no cadence is acceptable

`docs/ACTUATION_STAGES.md`'s **fully present participant** mode is this coupling, and
`docs/WORKBENCH_AUTHENTICITY.md` 2.1 was the gap: the layer deformed, pushed back on its bone,
and nothing it computed reached the integrator. It reaches it now. What that is worth, and
what it costs, is below, and the headline is a negative: **on a contact that develops as fast
as this one, no sub-cycling cadence is acceptable, so the coupling runs at one solve per plant
step and is 5.4x short of real time.**

Code: `SoftTissueCoupling` in `ihm/assembly/soft_tissue_layer.py`; call site
`ArticulatedBodyPlant.__init__`/`.advance`; selection
`plant_options.SOFT_TISSUE_COUPLED`. Battery: `scripts/verify_soft_tissue_coupled.py`,
pre-registered in its own header in `5a511f4`, before any gated run. Runs:
`logs/verify_soft_tissue_coupled.run1-PASS.log` and `.run2.log` — **32 gates, 0 FAILED, both
runs, and run 2 reproduces run 1 gate for gate and number for number** except the wall clocks
and two recorded figures run 2 corrected (below).

### The coupling, in three parts

| | |
|---|---|
| **only where it matters** | every plant step, each coupled segment's layer gets ONE matrix-vector product over its own nodes, giving its minimum gap and its **core's** minimum gap against the support. Skin clear → exactly zero, **no force port at all**, no solve. Measured cost: **0.0000 ms median over 96 unloaded steps** (below the timer). This is not an approximation — it is the same predicate `solve` uses to decide there is nothing to do |
| **sub-cycling** | a loaded segment is re-solved when `resolve_interval_s` has elapsed and the reaction is HELD in between: the world force vector (the support's normal is world-fixed) at a centre of pressure carried as a station in the **segment** frame, so the station follows the segment and the moment it delivers changes as the segment moves |
| **warm starting** | each solve begins from the last solve's shape in the segment frame, which `solve` already takes as an argument like any other |

**Two gates are NOT sub-cycled**, because both are free: contact onset, and **bottoming out**.
A held reaction that sailed past the step where the rigid core entered the support would be a
coupling that hid its own failure.

**Nothing is swallowed.** `solve` raises when the core would enter the support and when an
element leaves the constitutive domain; both are re-raised with the segment named
(`SoftTissueBottomedOut`, `SoftTissueLeftDomain`), `articulated.py` catches neither, and the
step simply does not happen. Gate X9: the reference run terminates that way and **the plant's
coordinates after the raise equal the coordinates before it, bitwise.**

**The one wrench component a point force cannot carry.** The engine's port is a point force on
a body, so the layer's (force, moment) is delivered at `p = (F × M)/|F|²`, which gives
`p × F = M` exactly *when M ⊥ F*. Frictionless contact against a half-space puts every nodal
force along the normal, so `M·F` must be zero — **checked on every emit, not assumed** (X3:
2.10e-10 of the wrench scale over 8 solves). A wrench with a real axial moment raises rather
than losing it silently.

### The fixture, and why it is this one

**This is the 22-segment scaffold collapsing.** Not the body. Nothing here is a statement
about a human forearm.

Upright plant, `pelvis_ty = 1.03` (the scaffold's own near-stance height: 616 N under the left
foot at t=0 against 761 N of weight; the model's default 0.93 starts the foot 63.7 mm *inside*
the floor), no actuation, 10 ms steps, one coupled segment: **`ulna_l`**.

`ulna_l` was chosen by measurement, over a 200-step uncoupled collapse that counted, for every
anchored segment, the steps in which its skin is below the floor **while its own engine contact
element is not**:

| segment | steps skin below floor | of which the layer is the ONLY contact |
|---|---:|---:|
| `ulna_r` | 99 | **99** |
| `ulna_l` | 69 | **69** |
| `tibia_r` | 33 | 33 |
| `torso` | 129 | 93 |
| `hand_r` | 137 | 43 |
| `calcn_l`, `calcn_r` | 0 | **0** |

In that window the layer is the only thing touching that segment and there is nothing to
double-count against. Measured in the run: the engine's own force on `ulna_l` peaked at
**1.67e-04 N** while the layer carried up to 14.73 N — five orders down, so the double count
is not what any number here is measuring.

### **The foot cannot be the fixture, and that is a result** — measured

| | |
|---|---|
| at the scaffold's own stance the foot's **skin** sits **36.3 mm above the floor** while the source foot spheres carry 616 N | the spheres reach ~36 mm below the skin |
| in a free drop from `pelvis_ty = 1.06`, peak foot load **770 N**, the skin still clears the floor by **32.4 mm** | |
| → **a coupled layer on `calcn_l` returns exactly 0 N in every upright trajectory the scaffold can produce** | this is the skin bundle's own documented caveat, now measured for the deformable layer |

**What the 36.3 mm is made of — measured 2026-09-18, `docs/SEGMENT_CONTACT_SURFACES.md`,
`scripts/verify_skin_contact.py`.** The right-hand column above, *the spheres reach ~36 mm
below the skin*, is true and reads as a fact about the spheres. It is not: **35.700 =
15.407 + 20.293** at the stance pose. The spheres hold the calcaneus **15.407 mm** off the
floor under the body's whole weight, and this body's own soft-tissue depth map reads
**14.91 mm** median over the plantar band — the sphere proxy's loaded stand-off is *inside*
the specimen's own heel pad. The remaining **20.293 mm** is the skin seated **above** the
bone by the global similarity, which is not a thing a body can do. A skin at the declared
pad would sit at **+0.5 mm**. So the conclusion here stands — a coupled layer on `calcn_l`
returns 0 N — but the reason is the skin's placement, not the spheres, and it is repairable:
`segment_contact: skin_per_segment` seats the heel at **−10.360 mm**, which would put the
`calcn_l` skin **+5.079 mm** above the floor at this pose instead of +35.700. Still clear,
so this fixture result does not change; the margin does.

And if it did touch: at the flat stance pose `calcn_l`'s layer carries **0.119 / 0.430 / 1.278
/ 2.866 / 5.574 N** at 0.5 / 1 / 1.5 / 2 / 2.5 mm, and **leaves the constitutive domain at
3.0 mm** (min J 0.180) — its rigid core is 3.76 mm under its lowest skin point there. Body
weight is 761 N. **This layer cannot carry a standing body, and no cadence changes that.**

### Known answers, as printed (run 2; run 1 identical)

| gate | verdict | what it says |
|---|---|---|
| X0 | PASS | two identical uncoupled runs agree **bitwise** in every coordinate over 140 steps. Without this every comparison below is void, so the battery stops here if it fails |
| X1 | PASS | the canonical force round trip: `registration.force` returns the source body, station and force the coupling handed it — **4.16e-16 m, 2.49e-14 N** |
| X2 | PASS | **the coupling invents no force path**: a fixed wrench through the coupling and the same force passed by a CALLER through `forces=` give **bitwise identical** plants over 20 steps |
| X2b | PASS | CONTROL THAT CAN FAIL: the same comparison with the force **negated** DIFFERS, at step 1. Without it X2 could be comparing two plants that never felt the force |
| X3 | PASS | the delivered point force reproduces the layer's wrench, axial moment included: **2.10e-10** over 8 solves |
| X4 | PASS | the layer's own force balance on every coupled solve: worst **7.10e-09** |
| **X5** | **PASS** | **OFF IS OFF, BIT FOR BIT.** `None`, the uncoupled identity and the **coupled** plant agree bitwise in every coordinate for **81 steps**, up to the step the layer first emits a port. An unloaded coupled plant **is** the historical plant |
| X6 | PASS | **call it twice at the same input**: `advance_state` is a pure function of `(state, transforms, dt)` — bitwise, unloaded and loaded |
| X7 | PASS | two identical coupled runs agree bitwise over 105 steps and stop for the same reason |
| X8 | PASS | checkpoint from a **loaded** step, two steps, restore, two steps: bitwise. A coupling whose held reaction is not checkpointed fails this |
| X9 | PASS | the run raises `SoftTissueLeftDomain` at step 104, re-raises on the next call, and **the plant is left exactly where it was** |
| X9b | PASS | a pose with the rigid core 10.00 mm inside the support raises `SoftTissueBottomedOut` naming the segment and the depth — from the cheap per-step gate, not from a solve |
| X10z | PASS | KNOWN ANSWER on the staleness instrument: the replay reads **exactly 0.000000** at the per-step interval |
| X10b | PASS | CONTROL THAT CAN FAIL: holding the first reaction forever reads **0.9697**, which exceeds the 0.10 bar. If it did not, the metric could not see staleness and X10 would mean nothing |

### What staleness costs, and the cadence

`e_N` is the held reaction against the per-step reference, **on the reference run's own poses**
(so the comparison is not confounded by the trajectories diverging), peak-normalised. `d_N` is
the closed loop: each cadence's own run against the reference, as the coupled segment's
greatest displacement difference.

**The bar for `d_N` is measured, not chosen**: `d_u` is what the layer's OWN force uncertainty
does to the same plant — the same run with the reaction scaled by 0.85, the 15% the layer's
convergence study (CV7–CV10) leaves on the table. Staleness must not be worse than the model's
own error.

| interval | steps held | `e_N` (bar 0.10) | `d_N` (bar `d_u` = 1.913e-04 m) | solves | terminates |
|---:|---:|---:|---:|---:|---|
| **10 ms** | 1 | **0.0000** | **0.000e+00** | 8 | LeftDomain @ 104 |
| 20 ms | 2 | 0.2606 | 2.776e-04 | 4 | LeftDomain @ 104 |
| 50 ms | 5 | 0.7091 | 9.605e-04 | 3 | LeftDomain @ **108** |
| 100 ms | 10 | 0.9697 | 1.160e-03 | 2 | LeftDomain @ **108** |
| 200 ms | 20 | 0.9697 | 1.160e-03 | 2 | LeftDomain @ **118** |
| 500 ms | 50 | 0.9697 | 1.160e-03 | 2 | **BottomedOut** @ 119 |

**Every interval above 10 ms fails, and the first one fails by 1.45x.** At and past 100 ms the
error equals "hold forever" — the interval is longer than the 80 ms loaded window, so the
cadence never re-solves inside it, and by this repo's own rule those rows are VOID and count as
FAILED. The termination step and even the termination *reason* move with the cadence, which is
a qualitative change, not a loss of resolution.

**OUTCOME (c) of the three pre-registered outcomes: no sub-cycling cadence is acceptable at
these costs.** The coupling ships at one solve per plant step.

### Why, in one number that carries to a different contact

The reaction grew from 0.45 N to 14.73 N over five steps: a peak loading rate of
**1,472.86 N/s**.
A 10% bar on a 14.73 N peak therefore allows a hold of **1.0 ms** — a *tenth* of the plant
step. That is the transferable statement: **a cadence is set by how fast the contact develops,
not by how much the solve costs**, and
`τ_max = 0.10 · max|F| / max|dF/dt|`. A contact loading 50x more slowly would tolerate a 50 ms
hold; this one does not tolerate one step. Note what that also says about the coupling at
N = 1: the force is computed at the pose the step starts from and held through it, so even the
per-step coupling is 10x outside its own bar on this trajectory. **An explicit coupling cannot
be made accurate on a contact this fast by choosing a cadence; it needs an implicit one.**

### Cost, as measured

| | |
|---|---|
| unloaded segment | **0.0000 ms** median over 96 steps (below the timer) |
| one loaded solve, `ulna_l` (14,874 DOF, penetration to ~10 mm) | **676.9 ms** median, 1,393 ms max (run 2); 2,012 / 3,869 ms in run 1 on a 3x busier machine |
| the coupling, amortised over the whole 105-step run | **53.9 ms per 10 ms plant step → 5.4x short of real time** |
| the arm's whole wall clock per step, coupling AND plant | **221.7 ms** (run 2) — mostly the collapsing scaffold's own contact solve, **not** the coupling. Run 1 reported this number as the coupling's; see the correction below |
| what the reaction does to the scaffold at all | the coupled run's forearm ends **1.310e-03 m** from the uncoupled run's |
| peak reaction | **14.73 N**, against 761 N of body weight |
| peak RSS | 1,054 MB |

`ulna_l`'s solve is far more expensive than the 136–159 ms heel figure for two measured
reasons and no mysterious one: it is **14,874 DOF against the heel's 7,296**, and the collapse
drives it to ~10 mm of penetration where the load stepping multiplies, against the heel ramp's
4 mm.

### Two recorded figures run 1 got wrong, and the shape of it

Run 1 passed all 32 gates. Two of its RECORDED numbers were still wrong, and both were
corrected before run 2. **No bar moved, and no gate was rescored.**

1. **The amortised cost counted work the coupling never did.** Run 1 divided each arm's WHOLE
   wall clock by its steps and called it the coupling's per-step cost: **552 ms**. But a
   collapsing 22-segment scaffold with 28 contact elements costs most of that by itself (run 2:
   221.7 ms a step in total against the coupling's 53.9 ms). Same family as this repo's *"a ratio whose denominator the
   treatment also changes"* — here the **numerator** held work the treatment never did, and
   the error was **10x** and in the direction that made the coupling look worse than it is.
2. **A typed constant.** Run 1 typed `ulna_l`'s bounding radius into X3's scale. It is read off
   the coupling now. The value is identical; a number nobody transcribes cannot be
   mistranscribed.

### What the coupling does NOT do

- **It does not replace the engine's own contact.** A coupled segment that also carries an
  engine contact element is a **second path to the floor**. Nothing in this repo's option set
  removes a segment's engine contact per segment, so the double count is inherent and is
  *reported*: every frame carries the engine's own resultant on each coupled segment beside
  the layer's. On this fixture it is 1.67e-04 N against 14.73 N; on the foot it would be the
  whole of body weight.
- **It does not make the layer converged.** No force from this layer is good to better than
  about 10% (CV7–CV10). A cadence cannot fix a discretisation error, and `d_u` above is that
  error expressed as plant motion.
- **It is not real time**, and the cadence does not make it so.
- **It is not the participant mode.** One segment of one scaffold, against a frictionless
  rigid half-space, with the muscle still inside a rigid core.

---

## 2026-09-18: the heel was a seam, the depth was a median, and the solve is 10-20x faster

Three problems were stated for this layer: it is not mesh-converged, its core sits 30 mm
above the sole against a declared 18.6 mm, and one heel solve costs 360-3,600 plant steps.
Each was worked. The main finding came first, and it changes what the other two mean.

### The "heel" in every table below was a seam wedge, not the heel — WITHDRAWN

Every "heel" force in this file (M1/M2/M3, the budget table, the 17x/4.4x comparison, the
2/6/12 mm rows) pressed `calcn_l` in its bundle frame, with the plane measured from the
skin's lowest point. **That point is vertex 780 of `skin_calcn_l.obj`, on the joint-cap rim
where `calcn_l` is cut from the toes (x = 92.7 mm, the forefoot).** Within 5 mm of it the
tissue is at most 1.7 mm thick: exact distances, sampled every 1 mm. The contact nodes of
the old 6 mm solve sit at x = 83-99 mm. So those numbers measure a thin wedge of skin at a
partition seam, cantilevered off a midfoot core, not a heel pad under a calcaneus. They are
kept below as the record of what was run. **None of them is a heel force.**

The battery now also has a **heel fixture**: `calcn_l` rotated about its frame z axis by
the least angle that keeps every vertex of the forefoot half clear of the support at 6 mm.
It is found by bisection on the skin mesh, and no angle is typed; it comes out at
17.2552°. Contact then runs from x = −4.6 to 42.7 mm, which is the heel.

### Depth: the 30 mm was two different places; the median was the wrong number

**The 30 mm column is withdrawn.** It was the lowest CORE node, at the midfoot (x ≈ 38 mm),
minus the lowest SKIN point, at the forefoot seam (x = 92.7 mm): two places 55 mm apart,
where the sole's own height differs by about 10 mm. Measured properly (exact distances to
the triangles, 1 mm samples, the core's lowest point above the LOCAL sole), the three
candidates come out as follows.

| | core above the local sole, heel (x 10–40 mm) |
|---|---:|
| exact offset of the declared 18.61 mm (1 mm samples), by 10 mm bin | 18.6–21.4 mm |
| voxel rule, 5 mm (M4) | 21.1 mm |
| fitted surface, same median rule (M4) | 18.5 mm |

- **The voxeliser** adds up to one cell of quantisation: 30.85 against the exact 28.35 mm on
  the old two-place measure, and 21.1 against 18.5 mm locally. The fitted surface removes it.
- **The core-extraction rule** does what it declares: under the heel the core starts at h.
- **The declared depth is the wrong one.** It is ONE number per segment, the median of the
  depth map over the segment's 628 skin points. The same depth map, read locally and mapped
  into the segment frame, gives:

| site (calcn_l, plantar) | depth map, local | in vivo (source card) |
|---|---|---|
| heel, x 12–37 mm | 6.4–14.1 mm; the nearest structure is the calcaneus for 26 of the 30 points | 15.99 mm median, 9.60–17.74 (Teng 2022); 14.85 ± 2.81 (Yang 2022) |
| forefoot, x 60–100 mm | 3.9–8.9 mm | fat pad 5.8–7.1 mm at MTH1-5 (Sanchez-Rodriguez 2025); all plantar soft tissue 8.9–13.5 mm (Garcia 2008) |
| segment median | 18.61 mm (set by the sides and back of the heel) | — |

Applied everywhere, 18.6 mm leaves the forefoot with **no core at all**: it is thinner than
2 × 18.6 mm everywhere (no point of it is more than 15.1 mm from the skin), so the
metatarsals are deleted and the tissue under them hangs off a midfoot core. That, and not
the voxeliser, is why the "heel" was so soft.

**Fix: the local depth rule.** `scripts/build_soft_tissue_local_depth.py` maps the depth
map into each segment frame through the bundle builder's own
`binding_registration`, and gives every skin vertex its depth. Two known answers are checked:
- **A1:** every depth point lands on a bundle vertex. It passes for 19 of 20 segments at
  ≤ 8.0e-10 m.
- **A2:** each segment's median is exactly its declared thickness, as a float equality. It
  passes for all 20.

A cell is core when it is deeper than the depth at its nearest skin point, measured to the
SKIN faces only: a joint cap is a cut, not skin, and the depth map has no point on it. With
the local rule the core sits 6.9 mm above the local sole under the heel and 5.9 mm under
the forefoot (M4). The forefoot figure agrees with the fat-pad literature.

The local heel depth is shallower than the in vivo pad by about 2x. The depth map measures
the NEAREST structure (it says so: "underestimates in folds"), and the in vivo value is
along the load axis. The modulus E_app was measured over that in vivo thickness, so pairing
it with a 7–13 mm layer makes this heel stiffer than the in vivo layer stiffness E/h, by
the thickness ratio. That is recorded, not fitted.

**A1 FAILED for `radius_l`, and the failure is real.** One of its 306 depth points lies on
canonical vertex 89206. That vertex's binding weight is 0.401 `radius_l` against 0.393
`ulna_l`, so the layer map (vertex argmax) counted it for the radius, while the bundle
(triangle argmax of the mean weight) gave all three of its triangles to the ulna. The two
partitions disagree on that vertex. The bar was not moved: the local rule is refused for
`radius_l` (D1 recorded FAILED, D1b checks the refusal), and the plant identities say so.

### The body-fitted surface, and two versions that failed first

`fitted_layer_mesh` cuts the voxel lattice by two level sets: the skin (signed distance to
the capped mesh), and the declared depth. It uses isosurface stuffing (Labelle & Shewchuk
2007): lattice nodes near a cut are warped onto the surface, every tet is cut, the tissue
side is kept, and every surface vertex is then projected EXACTLY onto the skin or onto the
depth level set. Its alpha, 0.24969, is their BCC value and is **borrowed**: this is a Kuhn
lattice, for which they prove nothing, so the mesh quality it gives is measured on every
build instead.

Two earlier versions moved the voxel surface nodes instead of cutting, and both failed on
the same object. A tet whose four nodes all land on one smooth surface is flat, and no
interior motion can fix it. On calcn_l at 5 mm those tets' volume ratios form a continuum
(−0.17 to 0.2 at the median), with no gap to cut at.
- **A penalty on distance**, raised in decades over a mesh-quality energy, stalled at
  0.14 mm on the skin and 1.07 mm on the interface after 353 s at 6 mm.
- **Exact projection with the interior relaxed** reached both surfaces exactly. But it
  kept near-flat slivers (volume ratio < 5e-7), which would stiffen the contact, and it
  carried only 43–49% of the local rule's thin-layer load.

| gate | result |
|---|---|
| G1 exact closest point vs brute force over every face, and a function | 6.9e-18 m, bitwise twice |
| G2 fitted heel (median rule, 5 mm): no flat tet, skin and interface nodes on their surfaces | 3.5e-17 m / 1.0e-17 m; min quality 6.0e-3 (1st percentile 0.19) |
| G3 built twice | bitwise |
| F1 every segment of `layer_fitted_local_*` builds (`scripts/census_soft_tissue_fitted.py`) | 19 of 19 named, every residual 0.000 µm |
| F2 radius_l refused by name and by `bodies=None`; calcn_l built fitted/local | pass |

Whole body under the local rule at 5 mm: **936,669 DOF, 1,668,019 tets, 23.58 L of
layer**, against 31.5 L for the voxel/median layer and 24.86 L of measured skin-to-bone/
muscle space. Peak memory for one segment is 922 MB (the torso), and the build takes 484 s
in all. The census's first run spanned the machine-wide OOM of 19:58 and read 800 s; it was
re-run, and only the re-run is quoted. Gates and counts were identical. The worst elements are poor: minimum quality 4e-8 (pelvis), with the 1st
percentile at 0.02–0.13 per segment. One consequence: the fast gradient agrees with the
shipped one only to 2.3e-11, not the pre-registered 1e-12. The worst element's Dm has
condition number 1.08e4, so two mathematically identical ways of forming Dm^-1 differ by
1e-11. S1 is recorded FAILED for that reason; with an identical inverse it reads 1.0e-12.

### Convergence

| fixture | spacing 6 / 5 / 4 / 3 / 2.5 mm | verdict |
|---|---|---|
| voxel, M3 (seam wedge), 6 mm | 1.113 / 0.457 / 0.657 N | non-monotone (unchanged) |
| fitted median, M3, 6 mm | 0.3050 / 0.2871 / 0.2926 / 0.2811 / not run (memory) | CV1, CV2 FAILED: the wedge is thinner than a cell, and the mesh's lowest node stays 1.5–2.4 mm off the skin |
| fitted local, M3, 6 mm | uncomputable at every spacing | CV3, CV4 FAILED: the local core at the seam is shallower than 6 mm (bottomed out, or J < 0.2) |
| fitted local, M3, 2 mm | 0 / 0 / 0 / 0.0940 / 0.0964 N | CV5, CV6 VOID (the plane does not reach the mesh until 3 mm) |
| **fitted median, HEEL, 6 mm** | **27.595 / 27.194 / 26.409 / 25.653** / not run (memory) | **CV7 PASS** (monotone); **CV8 FAILED** (differences 0.400 / 0.785 / 0.756 N) |
| **fitted local, HEEL, 2 mm** | **6.283 / 5.974 / 5.536 / 5.290 / 5.165** | **CV9 PASS** (monotone); **CV10 FAILED** (differences 0.309 / 0.438 / 0.246 / 0.125 N) |

On the heel, the body-fitted surface did what it was for: the force now moves in one
direction as the mesh refines, where the voxel force jumped with the staircase. It is **not
converged**. The differences shrink only after the first step, so the strictly-shrinking
gates fail as written. From the last three spacings the observed order is 1.15 (median)
and 0.94 (local): first order. The Richardson estimates are 23.7 N and 4.50 N, so the
finest forces are still about 8% and 15% above their extrapolated limits.

(Run 6 printed these estimates as 27.58 and 5.83 N. The sign of the correction was
reversed in the instrument; the observed order was right. Both are reported, never gated,
and run 7, with the corrected formula, prints 23.72 and 4.50 N.)

The force falls with refinement, which is what over-stiff linear tetrahedra near
incompressibility (ν = 0.45, volumetric locking) would do. That is a hypothesis, untested
here. **No contact force from this layer is a converged number to better than ~10%.** The
2.5 mm median heel did not fit the 4 GiB budget; the 2.5 mm local heel did (47,340 DOF,
1.19 GB).

### Every gate added 2026-09-18, as printed (`logs/verify_soft_tissue.run6.log`; run 7 reproduces it gate for gate and force for force)

All were pre-registered in the battery's header, in four commits made before the results
they judge existed (e140edf, 48af5cc, a282142, 8a5cf2f; the census in f615ab6). **Every gate
that existed before still passes with the same printed numbers**: K1–K12d and C1–C3, with R1
still recorded FAILED at 67.9%. Run 5 crashed at G2 on a key the new mesh no longer reports;
its log is kept as `run5-G2-KeyError`.

| gate | verdict | what it says |
|---|---|---|
| G1 G2 G3 D1b F1 F2 | PASS | exact geometry; fitted surfaces exact; refusal works |
| D1 | FAILED (recorded) | radius_l's local depth: the two partitions disagree on one vertex |
| S1 | **FAILED** | fast gradient vs shipped: 2.33e-11 against a 1e-12 bar; cause, element conditioning (1.08e4), shown above |
| S2 | PASS | H·v vs assembled Hessian, 1.31e-11 |
| S3, S6, S3m, S6m, S4, S5 | **FAILED (VOID)** | on M3's fixture the fitted meshes are unloaded at 2 mm and the local core is reached at 6 mm: nothing to compare |
| S3h, S4h, S5h, S6h, S6h-median | PASS | on the heel: fast = Newton to ≤ 1.4e-7 N (bound 6e-5), bitwise through a cache clear, warm = cold, faster at every depth |
| S3h-median | **FAILED** | at 12 mm the reference Newton does not converge (stagnates at 1.3e-6 N; 5.7e-7 with COLAMD), so "both converged" is false; the forces agree to 1.4e-7 N |
| CV1–CV6 | **FAILED** (CV5/6 VOID) | M3's seam-wedge fixture cannot be refined to convergence at ≥ 2.5 mm |
| CV7, CV9 | PASS | heel force monotone in the spacing |
| CV8, CV10 | **FAILED** | differences do not shrink from the first step |
| RT1, RT1m | FAILED (recorded, VOID) | M3's fixture ramps start unloaded |
| RT1h, RT1mh | FAILED (recorded) | heel ramps: 147 / 135 ms and 330 / 295 ms median per step (runs 6 / 7) |

### Speed, change by change (fitted-local heel, heel fixture, 7,296 DOF; `scripts/measure_soft_tissue_speed.py`)

Machine: shared GB10, 20 cores. The speed runs ended 19:50 and 20:06; the second (quoted
here) ran after the 19:58 machine-wide OOM had cleared, and agrees with the first to within
load noise. Each configuration ran in its own process; peak RSS was 190 MB (261 MB for the
component timings).

**One Newton iteration, taken apart at the 2 mm solution:**

| change | cost before | cost after |
|---|---:|---:|
| gradient: one sparse scatter product instead of four `np.add.at` | 8.0 ms | 0.86 ms |
| ...of which the BLAS vector dot `volumes @ density` (OpenBLAS thread start-up on this machine; `np.sum` does it in 7 µs) | ~5 ms | 0 |
| Hessian: matrix-free H·v from F instead of element Hessians + COO assembly | 27 + 23 ms (116 ms projected) | 0.37 ms prepare + 0.94 ms per product |
| factorisation order: symmetric-mode MMD instead of COLAMD (also now in `method='newton'`) | 110 ms, 3.01 M entries | 45 ms, 1.56 M entries |
| **factorisation reused**: the rest stiffness is factored once per layer and free set, and reused as the CG preconditioner across iterations, load steps, calls and poses (a rotation only rotates it) | one factorisation per iteration | 1.05 ms per triangular solve |
| batched 3×3 products held component-major, and one sparse gather for F | 1.33 ms per product | 0.14 ms |
| line search: Φ's change computed from dF itself (exact `det(F+dF)-det F` expansion, `log1p`), so it resolves below the floating-point floor of the total | one ramp step: 16 it, 1.9 s, "stagnated" | 7 it, 0.37 s |

The symbolic factorisation proper is not reusable through scipy's SuperLU, and CHOLMOD is not
installed. What is reused is the whole NUMERIC factor: it is computed once, and every
iteration costs a triangular solve instead of a factorisation.

**Whole solves, cold** (fitted-local heel, heel fixture, method='newton' vs method='fast',
same process, same answer to 1e-9 N: S3h):

| penetration | force | newton | fast |
|---:|---:|---:|---:|
| 1 mm | 1.3397 N | 2.16 s (10 it) | 0.21 s (7 it) |
| 2 mm | 5.9739 N | 1.96 s (15 it) | 0.17 s (10 it) |
| 4 mm | 31.467 N | 3.34 s (25 it) | 0.44 s (22 it) |
| fitted-median heel, 12 mm (battery) | 100.99 N | 35.1 s, and NOT converged (stagnated at 1.3e-6 N) | 3.33 s, converged to 1.4e-10 N |

**The plant loop.** The call a caller makes each 10 ms step is the next section. Here it is
warm-started from the previous step's shape along a 0 → 4 mm ramp at 0.5 mm per step:

| force tolerance per DOF | median per step | max | force change vs the default |
|---|---:|---:|---:|
| 4.4e-9 N (default, 1e-8 μ h²) | 159 ms (warm steps 136 ms) | 212 ms | — |
| 1e-6 N | 94 ms | 153 ms | ≤ 1.5e-7 relative |
| 1e-4 N | 93 ms | 111 ms | ≤ 2.4e-6 relative |

At 2 mm per step: 176 ms, then 308 ms. **The plant loop is 9–16x away, not real-time.**
RT1h (the battery's own ramp) is recorded FAILED at a 147 ms median (run 6), 135 ms (run 7). Loosening the
tolerance buys little: most of the 5–8 Newton iterations per step go to the contact set and
the nonlinearity, not to the last digits.

**The reduction that would be real-time, and why it is not offered.** The rest stiffness,
condensed onto the 154 candidate contact nodes, with contact as a bound-constrained QP (the
linear Signorini problem), costs 1.3–3.2 ms per step to 3.5 mm, plus a 0.2 s setup. Its
naive active set took 1.17 s at 4 mm, where it adds contacts one at a time. It is linear,
so it cannot pass K2/K3 at finite strain, and on this heel its force is **−1.8% at 0.5 mm,
−12% at 2 mm and −27% at 4 mm** against the nonlinear solve. A reduced model that fast and
that wrong at stance-sized strains is not a layer.

### The call, for whoever integrates it — **DONE (see the first section)**

`ArticulatedBodyPlant` now makes this call itself when the selection asks for a coupled
identity. What follows is the raw call, kept because it is what the coupling does internally
and what an offline caller still uses.

```python
from ihm.assembly.plant_options import resolve_fidelity
from ihm.assembly.soft_tissue_layer import build_selected_layers
_, sel = resolve_fidelity(root, {'soft_tissue': 'layer_fitted_local_confined'}, environment='upright')
heel = build_selected_layers(root, sel, ['calcn_l'])['calcn_l']       # ~1 s, once
state = None
# every plant step, with the segment's world rotation R (3x3) and origin t (3,):
r = heel.solve(rotation=R, translation=t, plane_axis=1, plane_value_m=floor_y, plane_sign=1.0,
               method='fast', warm_start_local_m=None if state is None else state['positions_local_m'])
state = r['state']
# apply r['segment_force_n'] at the segment origin, plus r['segment_moment_nm'] (about t)
```

`solve` raises `ValueError('Bottomed out ...')` when the rigid core itself would enter the
support. That is the tissue saying it cannot carry the pose, and the caller has to handle
it. An unloaded segment returns exact zeros in microseconds (K4). Measured cost: 136–159 ms
per loaded step at 7,296 DOF and 190 MB.

---


The sections from here on describe the layer as first built (commit a1b2812), and are kept
as written except where marked.

---

## What it is

| | |
|---|---|
| **geometry** | the segment's own closed skin surface (`segment-contact-meshes/skin-layer-map-v1`), voxelised by the repo's own mesher `material_domains.voxel_partition`: six tets per cube, so no element is born inverted |
| **layer / core split** | a cell whose centre is deeper than the segment's **measured** soft-tissue depth `h` (the layer map's own per-segment median, `soft-tissue-depth-v1`) is rigid core; the rest is layer. Every node of a core cell is carried by the segment |
| **material** | `DeformableRegion`'s compressible neo-Hookean `W = mu/2 (I1-3) - mu ln J + lam/2 (ln J)^2`, the same law as the shipped foundation column. `E` is the layer map's in-vivo apparent modulus, `nu` = 0.45 is the bundle's declared value, and `rho` = 950 kg/m3 is adipose from `tissue_materials.DENSITY` |
| **contact** | a frictionless rigid half-space as an **exact bound** on every free node, with no penalty stiffness to choose |
| **coupling back** | the reactions at the carried nodes, summed, are the force and moment on the segment. A frame-indifferent energy's internal forces sum to zero and carry no moment, so `segment + crop - contact + body = residual` is an identity, and it is reported on every call |
| **solve** | projected Newton on the analytic Hessian of `DeformableRegion`'s own energy (checked against finite differences of the shipped gradient), with load-stepped bounds. L-BFGS-B (`DeformableRegion`'s own minimiser) is kept as the cross-check |
| **time** | quasi-static by default (the layer is massless and the segment carries all mass, as today), or backward-Euler incremental potential with lumped mass (`dt_s`) |
| **purity** | `solve` is a function of its arguments: a cold start every call, never the last answer. The battery checks this bitwise |

### The modulus: a measurement that brackets a material rather than fixing it

The in-vivo source card (`data/sources/in-vivo-soft-tissue-compression.json`) says the heel
modulus is *pressure over thickness strain*: a **layer** modulus. A 3-D material needs a
Young's modulus, and a layer modulus only brackets one:

- `confined` (default): `E_young = E_app (1+nu)(1-2nu)/(1-nu)`, so `lam + 2mu = E_app`. Under
  a load much wider than the layer is thick, the 3-D layer has exactly the measured stiffness
  `E_app/h`, and it reproduces the shipped foundation column exactly in the limit where that
  column is right. This is the softest reading.
- `unconfined`: `E_young = E_app`. At nu = 0.45 this is 3.79x stiffer under a wide load. This
  is the stiffest reading.

Which reading is right depends on the in-vivo footprint-to-thickness ratio, which the card does
not report. Both are offered as separate identities, and neither is fitted.

---

## Known answers: what the battery prints

Every bar was fixed in the script header before the first run, except the ones marked there as
MEASURED or RECORDED.

| | case | result |
|---|---|---|
| K1 | Hessian vs central differences of the **shipped** `energy_gradient` | rel **2.79e-11** |
| K2 | confined uniaxial compression vs **`supine_contact.foundation` called directly**, 10/30/50% | worst rel **6.99e-12** |
| K3 | unconfined uniaxial vs compressible neo-Hookean with zero lateral stress (force **and** lateral stretch) | worst rel **8.49e-10** |
| K4 | null case: a support that does not reach the layer, and one exactly touching | force, moment and displacement **exactly 0.0**, no solve |
| K5 | heel 2 mm into the floor: segment force = contact force; segment moment = contact moment | rel **2.26e-8** / **2.41e-8** |
| K6 | call it twice: box solve, heel solve, heel mesh build | **bitwise identical** ×3 |
| K7 | released box, backward Euler, 40 steps: kinetic + elastic never rises | worst step **-8.92e-10** of E0 (monotone) |
| K9 | projected Newton vs L-BFGS-B, unconfined box | rel **3.23e-8**; Newton 37 it / 0.12 s, L-BFGS-B 793 it / 3.64 s |
| K11 | the unanchored list in `plant_options` re-derived | exactly `patella_l`, `patella_r` |
| K12 | `None` is the historical plant; `soft_tissue` never changes the plant's kwargs; path-as-id, boolean, free environment, tampered mapping and unanchored segment are all refused; the built layer carries the bundle's own `h` and `E` | all pass (5/5 refused) |

**Controls that must fail, and did:**

| | control | result |
|---|---|---|
| C1 | Lamé parameters swapped | K2 fails at rel **0.595** |
| C2 | solve starved to one iteration | K5's balance gate fails (relative residual **1.00**) |
| C3 | explicit symplectic Euler, same state and dt as K7 | K7's gate fails: an element inverts and the energy diverges |

At 30% strain, tissue free to bulge carries **0.237** of the confined force (K3). That one number is
the whole physical difference between this layer and the column it generalises.

### Three gates that FAILED on the way, and what each one caught

The battery ran four times; run 4 (`logs/verify_soft_tissue.run4-PASS.log`) passes every gate, with R1 recorded as FAILED. **The bars did not move. Instruments did, and each change is
recorded where it was made.** The logs are kept: `logs/verify_soft_tissue.run1-K5-FAILED.log`,
`run2-K7-crashed.log`, `run3-C2-FAILED.log`.

1. **Run 1, K5 FAILED: the heel transmitted 0.000000 N.** Newton had converged to 5.17e-9 N
   against a 4.38e-9 N tolerance and then stalled at the floating-point floor, because Armijo
   accepts a zero-length step when Φ comes out exactly equal. It ran to the iteration limit.
   The load-step loop then **reported its start**, which was the rigid rest state sitting inside
   the floor, carrying nothing. A receipt described a state no solver had produced. Fixed: a
   failed solve now reports its last iterate, flagged `converged=False`. Stagnation is detected
   (no strict decrease), and it counts as convergence only within 100× the tolerance. The
   achieved residual is always in `free_residual_n`.
2. **Run 2 crashed at K7.** Dynamic mode had never been run through Newton, and the lumped mass
   ((n,1)) was added to a 3n×3n matrix without being expanded.
3. **Run 3, C2 FAILED: a starved, unconverged solve read "balanced" at 0.00e+00.**
   `balance_force_relative` divided by the contact force and defaulted to 0 when that force was
   zero, so a solve that never reached contact certified itself. K5 was protected only because
   it also demands contact > 0. **The control caught an instrument that could pass for a reason
   unrelated to what it tests.** The metric now divides by the largest force in the balance, and
   returns inf when there is a residual and no force.

---

## The measured budget

Machine: shared GB10, 20 cores, about 18–36 GB available during the runs. Every run was capped
at `prlimit --as` 4 GiB.

**Whole body at 5 mm** (K11, from the loop that builds each layer):

| | |
|---|---|
| free DOF | **855,468** |
| tetrahedra | **1,514,082** |
| layer volume | 31.5 L (voxel staircase; the declared skin layers are 11.75 L, and the measured skin-to-bone/muscle space is 24.86 L) |
| largest segments | pelvis 192,870 DOF, torso 183,729, femurs 92–94k |
| build (voxelise + split), one-time | torso 49.8 s, pelvis 17.0 s, heel 0.3 s |
| peak RSS of the whole battery | **542 MB** (run 4; 536 MB in run 3; the torso voxelisation dominates) |
| unanchored | `patella_l`, `patella_r`: the skin patch is thinner than its own measured depth everywhere, so there is no core and they are refused |
| cells through the layer | 1.5–1.6 for hands and toes (h 7.7 mm): **under-resolved**; 3–4.7 elsewhere |

**Per-solve cost, `calcn_l` pressed at its SEAM WEDGE (not the heel: see the first section), 8,256 DOF,
11,922 tets, quasi-static, confined reading:**

| physical penetration | force on segment | Newton it | wall, run 4 (run 3) |
|---:|---:|---:|---:|
| 2 mm | 0.117 N | 8 | 3.6 s (5.3) |
| 6 mm | 0.457 N | 28 | 9.6 s (10.3) |
| 12 mm | 7.783 N | 95 | 36.2 s (34.7) |

Wall time moves by up to ~1.7 s between runs on this shared machine; forces and iteration
counts are identical to the last printed digit.

**The cost per Newton iteration is the sparse LU factorisation:** 0.21 s with COLAMD, of which
energy, gradient, Hessian and assembly are about 0.06 s. L-BFGS-B on the same heel took 45–50 s
per solve (940–1,037 iterations). Memory for one heel solve is about 100 MB.

**Against the budget that matters:** the plant advances in 10 ms steps. One heel solve cost
360–3,600 plant steps of wall time. **Superseded (2026-09-18):** the fast path reuses one
factorisation and brings a warm-started heel step to 136–159 ms, still 9–16x the plant step.
The condensed small-strain compliance was built and measured: 1–3 ms, but 2–27% low. See
the first section.

---

## What it says about the heel, against the baselines — WITHDRAWN AS A HEEL RESULT (2026-09-18)

**Every row and cause below was measured on the seam wedge at the calcn/toes cut, not on the
heel** (first section). The table is kept as the record of what was run. Its three causes
are superseded: (1) the "30 mm column" was two places 55 mm apart, and under the heel the
core sat at the declared depth; (2) lateral escape was never separated from the missing
forefoot core; (3) the staircase is gone from the fitted layer. On the actual heel (heel
fixture, fitted, 5 mm) the layer carries 3.31 N at 2 mm and 27.2 N at 6 mm under the median
depth, and 5.97 N at 2 mm under the local depth. That is 28x and 51x what the table below
says for "2 mm".

Same physical penetration, measured from the **skin mesh's** lowest point (not the voxel
surface), 5 mm cells:

| penetration | 3-D confined | 3-D unconfined | shipped linear `k=E/h` | NH column, h=18.6 declared | NH column, h=30.0 (the voxel model's own column) |
|---:|---:|---:|---:|---:|---:|
| 2 mm | 0.117 N | 0.444 N | 0.37 N | 0.38 N | 0.23 N |
| 6 mm | 0.457 N | 1.733 N | 13.74 N | 17.28 N | 9.75 N |
| 12 mm | 7.783 N | 29.521 N | 129.94 N | 225.42 N | 108.06 N |

**The 3-D heel is 17× softer than what the plant ships (confined reading) and 4.4× softer
(unconfined reading).** Three causes, and they can be separated:

1. **The column is 30.0 mm, not 18.6 mm.** On a segment this small (bbox 115 × 105 × 88 mm), the
   depth-defined core is set by the nearest of *all* surfaces, including sides and joint caps, so
   the lowest core node sits 30 mm above the plantar skin. The neo-Hookean column at 30 mm
   carries 108 N against 225 N at 18.6 mm: depth alone is a factor of 2.1.
2. **Lateral escape.** Tissue in an isotropic continuum bulges. At matched column depth the
   stiffest reading still carries 3.7× less than the column.
3. **The voxel staircase.** The voxel surface's lowest point is 0.85 mm above the mesh's, and at
   12 mm only 76 nodes (≈19 cm²) touch, against 32.6 cm² of mesh triangles below the plane.

**And cause 3 dominates the mesh study, which does NOT converge (M3; the fitted surface's own study is in the first section):** at 6 mm physical
penetration the force reads **1.113 / 0.457 / 0.657 N at 6 / 5 / 4 mm spacing**, non-monotone,
while the voxel surface's lowest point jumps 12.0 / 10.0 / 12.0 mm against the mesh's 9.15. Where
the staircase lands decides how far the floor actually reaches into the tissue. **No contact
force from this layer is a converged number.** The fix is a body-fitted boundary: snap the
voxel surface nodes onto the skin triangles. It is named here and was not done. **Done
2026-09-18**: snapping failed twice, and cutting (isosurface stuffing) worked. See the first
section.

**The anatomical reading, as a hypothesis with its test:** a real heel pad is fat in closed
chambers walled by fibrous septa, and the septa resist the lateral escape that an isotropic
continuum allows freely. That would explain why a homogeneous 3-D model built from the in-vivo
layer modulus reproduces the measurement only in its confined limit, and why the 1-D column,
which forbids bulging, is closer to the in-vivo number for the heel. The test: with a
body-fitted surface, measure the 3-D heel's apparent layer modulus (pressure over thickness
strain) under a heel-sized footprint and compare it with 192.55 kPa. If the gap survives, the
constitutive model needs anisotropy (septal confinement), not a better modulus.

---

## The reduction that failed: the Saint-Venant crop

`crop_radius_m` holds every free node farther than R from the contact region at its rigid rest
state. On the heel it **fails (R1, recorded):** at R = 3h, 6 mm, **67.9% of the load lands on
the truncation, not the segment, and the segment force is 58.1% wrong**. At 2h the segment
receives essentially nothing. Cause: the contact is the seam wedge (first section), whose
nearest core is at the midfoot, more than 30 mm away, so a crop around the contact never
reaches the thing that carries the load. R1 still runs on that fixture and is unchanged. The crop is kept because it may hold on
large segments where the core is near the skin (the supine pelvis and torso, which is where it
was meant to pay), but **nobody has measured that yet**. Every cropped call reports
`crop_load_fraction`, and a cropped force should not be believed until that fraction is small.

The bar (≤ 1% at 3h) was written after development runs had already shown 17–100%, so R1 records
an observed failure, not a tested prediction. The script says so.

---

## What it is NOT

- **A LAW fitted to it IS in the native integrator; the layer itself is not.** `skin_layer_fitted`
  puts the elastic foundation's stiffness, calibrated offline against this layer, inside the
  engine's own error-controlled step. That is a fit to a ~10%-converged layer, it does not
  deform, and it is refused on the one segment where it failed its bar. See the third section.
- **Still not in the native integrator, but now in the plant's LOOP.** The engine integrates
  the rigid scaffold and its own contact, and the layer is not one of its force elements. What
  changed 2026-09-18 (first section): a coupled identity (`*_coupled`) makes
  `ArticulatedBodyPlant.advance` pose each coupled segment's layer at that segment's own
  transform every step and apply what it transmits through the same `forces=` port a caller
  uses. An uncoupled identity still flows back **nothing**, and is bit-for-bit the historical
  plant. The coupled one costs 53.9 ms per 10 ms step on one 14,874-DOF segment, and **no
  sub-cycling cadence was acceptable** on the contact it was measured against.
- **Not muscle.** Everything deeper than `h` is rigid, which is the shipped layer map's own
  assumption.
- **Not bone-anchored.** The bone meshes would be the right core, but the skin bundle records
  `bone_vertices_inside_skin` **0.444** and `segments_enclosing_their_bone` **0** (0.084 for
  `calcn_l`). A bone core would be built on a known registration mismatch.
- **Not continuous across segments.** Each layer is independent, with a seam at every segment
  boundary, as in the skin bundle.
- **Not converged at the contact surface.** The fitted heel force is monotone in the spacing,
  at first order, and ~8–15% above its extrapolated limit at the finest spacing affordable
  (CV7–CV10). It is also **not frictional**.
- **The dynamic mode is not energy-accurate.** Backward Euler at 1 ms takes the released box from
  1.55e-2 J to 3.4e-11 J in 40 steps. That is monotone, as K7 demands, and very dissipative. Its
  lumped mass is also *shadowed*: it duplicates mass the rigid segment already carries.
  Quasi-static is the default for that reason.
- **Not stance.** At the confined reading the heel carries 7.8 N at 12 mm. Body weight is 761 N.

---

## Why each remaining piece is blocked, specifically

| wanted | blocked by | kind |
|---|---|---|
| ~~the layer in the plant's integration loop~~ | **DONE 2026-09-18**, first section. What remains is cost (53.9 ms per 10 ms step on one segment) and the fact that no cadence amortises it: the contact develops at 1,473 N/s and a 10% bar allows a 1.0 ms hold | — |
| ~~an **implicit** coupling (the layer inside the plant's own step, not held across it)~~ | **ANSWERED SIDEWAYS 2026-09-18**, and the answer is not to put the layer in the step. The engine's force port takes a force and not a stiffness, so the layer itself cannot be implicit through it; what CAN be is a law CALIBRATED to the layer. `ElasticFoundationForce` fitted offline reproduces the heel layer to 4.8% (`calcn_l`) and 6.1% (`calcn_r`) held out, inside the layer's own 10% bar, and is solved inside the integrator at no per-step cost. It FAILS on `ulna_l` at 39.8% and that segment is refused a fitted number. What remains wanted is the layer itself implicit, which still needs a stiffness port | engine interface (for the layer); **done, with a measured domain, for the law** |
| a coupled foot | the scaffold's source foot spheres reach 36 mm below the skin, so the layer never touches in any upright trajectory; and at the stance pose it carries 5.57 N at 2.5 mm and leaves the constitutive domain at 3.0 mm against 761 N of weight | geometry + model |
| coupling without a double count | no option in this repo removes ONE segment's engine contact, so a coupled segment that has an engine contact element carries both. Visible in every frame, not removable from here | plant options |
| real-time cost | Python/numpy at 7,296 DOF: triangular solves and element kernels take ~1 ms each, times 5–8 Newton iterations and 30–60 CG steps. The linear condensed reduction is fast but 2–27% wrong. Would need a compiled solver, or a nonlinear reduced basis that passes K2/K3 | compute (buildable) |
| a converged contact force | staircase removed (fitted surface); what remains is first-order convergence of linear tets, ~8–15% at the finest affordable spacing. Would need quadratic or mixed elements, or a finer mesh than 4 GiB allows (2.5 mm median heel did not fit) | discretisation (buildable) |
| local depth for `radius_l` | the layer map (vertex argmax) and the bundle (triangle argmax) disagree on vertex 89206 | data / partition |
| a bone-anchored core | skin/bone registration: 0.444 enclosure, 0 of 20 segments enclose their bone | data / registration |
| patellae | skin patch thinner than its own measured depth at 5 mm | geometry at this resolution |
| hands and toes resolved through their thickness | 1.5–1.6 cells at 5 mm; 2.5 mm cells would multiply their DOF about 8× | compute |
| muscle as a deformable, actively stiffening layer | `h` stops at the nearest bone **or muscle**; there is no muscle-volume constitutive model or activation coupling here | model not built |
| septal confinement of the heel pad | an isotropic continuum cannot express it; it needs anisotropy with a measurement to set it | model + data |
| fascia as a glide plane (TISSUE_MECHANICS.md) | a layer exists now, but a glide plane is **contact between two soft surfaces**, and this layer only contacts a rigid half-space | solver (sliding soft-soft contact) |
| cross-segment continuity | one independent layer per rigid segment; the skin bundle is cut at every joint | geometry |

**Pain and afference are downstream of all of this:** the layer's per-node strain and pressure
are exactly what a cutaneous or nociceptive transducer would read (ACTUATION_STAGES.md,
"blocked on a receptor"), and nothing reads them yet.
