# the soft body: a deformable tissue layer over the rigid skeleton

`docs/ACTUATION_STAGES.md` defines two modes, `driven` and **fully present participant**,
and says the second is to exist now. Before this file, nothing in the stack could deform:
every compliant element was a 1-D law on a rigid carrier. The supine foundation is one
independent confined column under each quadrature point, and the segment contact meshes are
an elastic foundation with one spring per triangle. Neither can bulge, share load with a
neighbour, or tell a heel-sized footprint from a fingertip-sized one.

This is the first piece of the participant mode: **a 3-D soft-tissue layer per segment that
deforms under contact and returns the load it transmits to its bone.** It is real, verified
against analytic answers, and wired so a live body can ask for it. It is also slow, not
mesh-converged at the contact surface, and softer than every measurement it was built from.
Each of those is measured below. **It is not the participant mode.** It is the part that was
buildable from what was already on disk, with the reason each remaining part is blocked.

Code: `ihm/assembly/soft_tissue_layer.py`. Selection: `ihm/assembly/plant_options.py`
(`soft_tissue`). Battery: `scripts/verify_soft_tissue.py`, which writes
`data/derived/soft-tissue-layer-v1/report.json`.

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

**Per-solve cost, heel (`calcn_l`, 8,256 DOF, 11,922 tets), quasi-static, confined reading:**

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

**Against the budget that matters:** the plant advances in 10 ms steps. One heel solve costs
360–3,600 plant steps of wall time. **This is not real-time and is not close.** A
participant-mode heel strike at 100 Hz would need two to three orders of magnitude, which means
a compiled solver with a reused symbolic factorisation (CHOLMOD is not installed here; scipy has
only SuperLU), a small-strain condensed compliance on the surface nodes, or a reduced basis.
None of these was built.

---

## What it says about the heel, against the baselines

Same physical penetration, measured from the **skin mesh's** lowest point (not the voxel
surface), heel, 5 mm cells:

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

**And cause 3 dominates the mesh study, which does NOT converge (M3):** at 6 mm physical
penetration the force reads **1.113 / 0.457 / 0.657 N at 6 / 5 / 4 mm spacing**, non-monotone,
while the voxel surface's lowest point jumps 12.0 / 10.0 / 12.0 mm against the mesh's 9.15. Where
the staircase lands decides how far the floor actually reaches into the tissue. **No contact
force from this layer is a converged number.** The fix is a body-fitted boundary: snap the
voxel surface nodes onto the skin triangles. It is named here and was not done.

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
receives essentially nothing. Cause: the core is 30 mm above the contact, and a crop around the
contact never reaches the thing that carries the load. The crop is kept because it may hold on
large segments where the core is near the skin (the supine pelvis and torso, which is where it
was meant to pay), but **nobody has measured that yet**. Every cropped call reports
`crop_load_fraction`, and a cropped force should not be believed until that fraction is small.

The bar (≤ 1% at 3h) was written after development runs had already shown 17–100%, so R1 records
an observed failure, not a tested prediction. The script says so.

---

## What it is NOT

- **Not in the native integrator.** Selecting `{'mechanical_fidelity': {'soft_tissue':
  'layer_map_confined'}}` gives the live body a disclosed identity: `articulated.py` already
  writes the whole selection to `mechanical_fidelity.json`. `build_selected_layers` builds the
  layers from it. The engine still integrates the rigid scaffold and its own contact.
  **Nothing flows back into the plant** unless a caller poses the layer each step and applies
  `segment_force_n` / `segment_moment_nm` as external loads. That caller would live in
  `articulated.py` or the native stream, which were outside this work's territory. At the costs
  above, closing that loop would cost tens of seconds per 10 ms step.
- **Not muscle.** Everything deeper than `h` is rigid, which is the shipped layer map's own
  assumption.
- **Not bone-anchored.** The bone meshes would be the right core, but the skin bundle records
  `bone_vertices_inside_skin` **0.444** and `segments_enclosing_their_bone` **0** (0.084 for
  `calcn_l`). A bone core would be built on a known registration mismatch.
- **Not continuous across segments.** Each layer is independent, with a seam at every segment
  boundary, as in the skin bundle.
- **Not converged at the contact surface** (M3), and **not frictional**.
- **The dynamic mode is not energy-accurate.** Backward Euler at 1 ms takes the released box from
  1.55e-2 J to 3.4e-11 J in 40 steps. That is monotone, as K7 demands, and very dissipative. Its
  lumped mass is also *shadowed*: it duplicates mass the rigid segment already carries.
  Quasi-static is the default for that reason.
- **Not stance.** At the confined reading the heel carries 7.8 N at 12 mm. Body weight is 761 N.

---

## Why each remaining piece is blocked, specifically

| wanted | blocked by | kind |
|---|---|---|
| the layer in the plant's integration loop | a caller in `articulated.py` / the native stream (not in this work's territory), and **cost**: 3.6–36 s per solve against a 10 ms step | territory + compute |
| real-time cost | no compiled sparse Cholesky here (SuperLU only); no reused symbolic factorisation; no condensed or reduced basis | compute (buildable) |
| a converged contact force | voxel staircase; needs a body-fitted boundary | discretisation (buildable) |
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
