# Skin rendering follows the contact material support

`SegmentSurfaceBinding` exposes the same hard material attachment used by the
world's sampled skin contacts. The whole skin no longer needs to inherit one
segment transform merely because it is one source mesh.

For every **canonical rest vertex**, choose the nearest named bone group's
axis-aligned rest envelope by Euclidean distance; ties choose the lexicographically
first native segment name. The first retained canonical bone in that group owns
the attachment. Subtract its reference centroid, rotate the resulting offset by
that bone's current rotation, then add its current centroid. Apply the common
canonical-to-world/view transform after this operation, once.

The mechanical frame carries frozen `surface_binding` metadata (including source
hashes, reference envelopes/centroids, segment and bone IDs, a stable
`binding_identity`, and `surface_entity_ids`) plus small per-frame
`surface_transforms` containing the 22 current segment centroid/rotation pairs.
The binding is fixed in rest coordinates and does not get recalculated from a
moving posed skin. There is no display smoothing or independent pose generator.
Respiration's display-only displacement cannot be added while retaining exact
agreement with the physical skin contact support.

Other extensive attached surfaces can use this same declared prior. A garment
with its own simulated positions must use those positions instead. At a
mixed-owner rendered triangle, nearest bound vertex selection is an approximate
force-picking rule; sending the selected bone ID preserves the intended segment
and material-point attachment. The frame adds neither forces nor body mass.

This is **not deformable skin FEM**. Hard assignment can expose seams at joints;
canonical anatomical registration is still approximate. The contract fixes the
render/contact disagreement without pretending to solve tissue continuity.

Verification: `PYTHONPATH=. .venv/bin/python scripts/verify_surface_binding.py`.
The retained receipt `data/derived/surface-binding-verification.json` checks all
102,467 full-resolution skin vertex owners against the original registration
ranking, then compares 1,193 contact samples in actual recorded native poses at
20 and 80 ms. Both canonical and world coordinates are bit-for-bit identical
(maximum error 0 m). The fixture `data/derived/surface-binding-fixture.json`
contains frozen metadata and 24 posed vertices for cross-language renderer tests.
These two short poses establish attachment parity, not long-duration stability or
anatomical validation.

## Continuous shared candidate

The hard attachment is retained as a comparison, **not a satisfactory skin
embedding**: an actual native pose at 80 ms stretches a 0.779 mm edge to 325 mm
where adjacent vertices chose hand and pelvis supports. The renderer parity test
correctly exposed this physical prior defect rather than concealing it.

`ContinuousSurfaceBinding` replaces the hard assignment with weights derived
from the retained skin topology. Seeds choose the nearest retained bone mesh
vertex (an explicit vertex-distance approximation to bone surface distance).
A screened weighted graph Laplacian diffuses those seed weights along skin
triangle edges, using an uncalibrated 80 mm length and 1 mm minimum edge length.
Exactly coincident vertices share supports; no distance-based weld joins nearby
body parts. The source has 100 raw connected components and 776 exact duplicate
vertices; the exact weld yields one connected graph with 101,691 vertices.
Negative solver roundoff is removed, tails below 1e-8 are discarded, and weights
are quantized to float32 then normalized once. The frozen sidecar is consumed by
both physics and display; poses never alter its weights.

The position of each skin point is the weighted sum of its segment-specific
posed material stations. A force **F** at that point scatters **wᵢ F** to each
segment at its own station **xᵢ**. This is the transpose of the velocity map. It
preserves total force, total moment about every origin, and virtual power;
applying those distributed forces at the blended point instead would generally
violate the rotational work relationship.

The sidecar `data/derived/canonical/continuous_surface_binding.json.gz` contains
source hashes, algorithm parameters, reference vertices, segment IDs and weights.
`asset_records()` supplies exact immutable JSON bytes keyed by SHA256;
`manifest()` declares `weights_url` and `weights_sha256`. Per-frame transforms
remain the small 22-segment pose set. Other attached surfaces transfer weights
from the nearest canonical rest skin vertex, breaking exact-distance ties by
source vertex index. Independently simulated cloth overrides this attachment.

Reproduce the artifact with:

```
PYTHONPATH=. .venv/bin/python scripts/build_continuous_surface_binding.py \
  --registration <retained-runtime>/mechanics/registration.json
PYTHONPATH=. .venv/bin/python scripts/verify_continuous_surface_binding.py
```

The receipt `data/derived/continuous-surface-binding-verification.json` measures
maximum edge stretch at 80 ms falling from 417.53 to 3.485 and the 99th percentile
from 12.82 to 1.474. Force, moment and virtual-power residuals are below 2e-15;
a shared rigid motion is reproduced within 3e-16 m. These are geometric and
coupling checks. The minimum edge ratio is still 0.0375 in that pose: substantial
local compression remains, and the method neither preserves tissue volume nor
prevents self-intersection. Native body/world regression is required before
promoting the changed contact embedding. This is not skin FEM or a calibrated
material law.

## The 25-body variant: a second blend (18 September 2026)

`data/models/articulated_spine_v1` (`docs/ARTICULATED_SPINE.md`) splits `torso` into `torso`,
`thorax`, `cervical` and `head`. The blend above is solved over the 22 base segments, and
`ContinuousSurfaceBinding.from_root` holds a blend to the exact segment supports it was fitted
on. So a variant `ArticulatedBodyPlant` died here with 'Surface binding registration supports
changed'. That was the second 22-body layer to refuse the variant. The first was the force-frame
registration, fixed in `7ded91b` (`_registration_for_native`).

**There is now a second blend. The 22-segment one is untouched**, and every result above was
measured on it.

    data/derived/continuous-surface-binding-articulated-spine-v1/continuous_surface_binding.json.gz
        sha256 2bfa36afcb27d57cceafdd8b8dcbeaeebffa8808cd7678fa18ecfd48c53de57b
        binding_identity 44ddc6757fcb74664e65171865259114e7244798e43e0c278bf8ee6b3de94b47

The plant picks the blend through `articulated.SURFACE_BINDINGS`, keyed by the anatomy plant
that `_registration_for_native` identified from the engine's body set. This is the same way
`anatomy_pose.PLANTS` picks a segment binding. The base entry is `DEFAULT_ASSET`, which is
exactly what `from_root` used to default to. `continuous_surface_binding.py` is unchanged.

**The same method.** The variant blend is built by the committed `materialize()`, via the same
builder, over the variant registration's 25 groups:

    # a variant plant writes registration.json + canonical_mechanics.json BEFORE it loads the
    # blend, so a construction that refuses for want of this asset leaves the builder's input
    PYTHONPATH=. .venv/bin/python scripts/build_continuous_surface_binding.py \
      --plant articulated_spine_v1 --registration <variant-plant-output>/registration.json

`--plant` sets only the default output and a check that the registration's bodies are that
plant's. The base default is unchanged. The build takes about 3 minutes and peaks at 448 MB RSS.
Two builds produced byte-identical `.gz` files.

**The method reproduces the committed BASE blend, but not bit-identically.** I rebuilt it from a
fresh base plant's `registration.json`. The segments and rest positions are identical, and all
102,467 vertices keep their argmax. The weights differ by at most **2.95e-8**, which is
conjugate-gradient roundoff: the maximum solver residual is 1.5236e-6 against the 1.5179e-6 the
committed artefact records. The committed artefact predates the retained `materialization_source`
field.

**The variant blend is a pure refinement of the base one** (`scripts/verify_variant_runtime.py`, C):

* its 21 non-torso segments have the base's exact supports;
* their weight columns match the base within 3.74e-8;
* `torso + thorax + cervical + head` sums to the base `torso` column within 3.72e-8;
* merged back over that family, the argmax agrees with the base on **all 102,467 vertices**.

This follows from the solve being linear. Each column solves the same screened operator against
its own one-hot seed. Splitting a seed set four ways splits its solution four ways. The
unchanged limb columns show that the new skull anchors claimed only vertices the base had
already seeded to `torso`.

Of the 25,704 vertices whose base argmax was `torso`:

| now argmax | vertices |
|---|---:|
| `head` | **12,218** |
| `thorax` | 11,820 |
| `cervical` | 900 |
| `torso` | 664 |

**Head and neck skin now rides `head`/`cervical`, not `torso`.** Under `head_extension` at its
declared upper limit (0.4189 rad):

* the 60,894 vertices with no head weight stay exactly still (0.0 m);
* each of the 41,573 head-weighted vertices moves by its head weight times the head station's
  motion, to 3.0e-16 m;
* the 12,218 vertices whose largest weight is `head` move 0.07–83.1 mm, following a median 99.4%
  of the rigid head motion.

The screened diffusion leaves every head vertex some cervical or thorax weight: the largest head
weight is 0.99947, and no vertex is wholly on `head`. Head weight above the 1e-8 tail cut
reaches 41% of the skin's vertices.

### What this does not fix

* **The girdle.** The force frame's `torso` is anchored on 9 lumbar vertebrae, because the base
  `torso` group never contained a scapula or clavicle. So shoulder skin blends between the
  humeri and `thorax`, while the humeri hang from `torso` in the model. When the thoracic joint
  moves, that skin is pulled between two parents. This is the same girdle caveat that
  `ANATOMY_SEGMENT_BINDING.md` records for the rigid binding.
* **The poser's skin.** `AnatomyPoser.skin_vertices` (`anatomy_pose.py`, `SKIN_BINDING`) still
  reads the 22-segment blend for every plant. For the variant, the display skin of the head
  therefore still follows `torso`, even though the display skull follows `head`. The fix is one
  line in `anatomy_pose.py` (select by plant, as `PLANTS` does), and `anatomy_pose.py` is outside
  this change. The force-frame skin, `surface_transforms` and the blend a client downloads are
  all the variant's.
* **Segment skin contact.** `mechanical_fidelity={'segment_contact': 'skin'|'skin_carried'}` does
  not refuse the variant, but it is wrong for it. The skin contact meshes
  (`data/derived/segment-contact-meshes/skin`) are the argmax partition of the 22-segment blend,
  and `skin_torso.obj` (the whole trunk, neck and head) is installed on the variant's lumbar-only
  `torso`. It needs meshes rebuilt from this blend. Separately, `joint_stops` stops 22
  coordinates and none of the variant's 15 new ones. Both are in `plant_options.py`, outside
  this change.

### Verification

    PYTHONPATH=. .venv/bin/python scripts/verify_variant_runtime.py
        receipt: data/derived/variant-runtime-verification.json         (3 engine sessions, one at a time)

* **A. Base plant bit-identical.** The current module is compared with `articulated.py` at
  `7ded91b`, which is loaded from git and run in the same process. The construction is upright,
  77.6122029 kg, `whole_body_arm26_v2`, display pose `opensim`. The comparison covers 7 frames
  (5 steps of 10 ms plus a re-projection) and all 6 output files, float by float: **0
  differences**. There are two controls. The pinned source must differ from the current one and
  lack `SURFACE_BINDINGS`. A single 1-ulp change to one centroid must be caught, and it is.
  Before the change, two runs of the committed code also gave 0 differences.
* **B. The variant runs.** It constructs on `articulated_spine_v1` with 25 bodies. The force frame
  registers all 25. It steps, and its frames carry 25 skin transforms and 25 display segments.
  `_project` and `_display_pose` called twice at the same state are bit-identical. It also ran
  with `display_pose='anatomical'` and with `environment='supine'`.
* **C. The blend.** The variant blend is loaded, and loading it twice gives the same identity.
  **The base blend is still refused for the variant registration** (the control that shows the
  selection matters). The refinement result is above.
* **D. The skull.** In the display pose the head's 89 named bones ride `head` and move 26.0–63.8
  mm under `head_extension`, and no other segment's motion changes (0.0). The same assertion on
  the base plant fails, as it must: there the 89 bones ride `torso`.
* **E. The skin in the force frame.** The numbers are above.

`verify_display_pose.py` and `verify_plant_fidelity.py` still ALL PASS. The display payload is
8,111 bytes a frame, the same under the `7ded91b` module. The 8,077 quoted in that commit's
message is not what that code produces now.

The variant runs, and it is not claimed to stand. Why its ankle collapses is a separate
investigation.
