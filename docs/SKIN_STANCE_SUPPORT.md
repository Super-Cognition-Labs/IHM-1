# Why the skin topples the body: separating WHERE the support is

`docs/SEGMENT_CONTACT_SURFACES.md` ("Result: the seat is repaired and it is NOT enough")
left a named hypothesis with its instrument unbuilt: with the per-piece skin replacing the
source feet, the body is 671.8 mm down at 1.00 s, and every skin vertex reaching the floor
is **behind** the centre of mass (−87.5 mm) because the hard partition puts the forefoot
skin on the `toes` body. The identical geometry with all 12 source foot spheres kept
stands (−1.3 mm). That correlate does not separate *where* the support is from *that
there are spheres at all* (their stiffness, their friction, their count).

The engine now keeps a SUBSET of the source spheres beside the meshes
(`segment_contact_replaces_source_feet='posterior'|'anterior'`, native replace modes 2/3,
`scripts/native_mechanical_stream.cpp`): posterior = `contactHeel_*` and
`contactLateralRearfoot_*` (4 spheres, all on `calcn`); anterior = the midfoot and toe
spheres (8, on `calcn` and `toes`).

## PRE-REGISTRATION (written and committed before any arm below runs, 2026-10-29)

**Arms** (`scripts/measure_segment_contact_meshes.py`, stance pose, joint stops, 10 ms):
`spheres`, `skin_per_segment_carried`, `skin_per_segment`, **`skin_per_segment_posterior`**,
**`skin_per_segment_anterior`**. Each is run for **100 steps (1.00 s) and 200 steps
(2.00 s)** — two seconds because twenty-five steps once read a fall as a settle, and one
second is where the known fall was measured, so two is the first horizon that can show
whether a 1 s stand is a stand.

**Outcome per arm:** *stands* = pelvis_ty drift above −20 mm at the horizon with no
integrator failure; *falls* = drift below −100 mm or integrator failure; anything between is
reported as *sagging* and not counted either way.

| `posterior` (support behind CoM only) | `anterior` (support in front) | reading |
|---|---|---|
| falls | stands | **supported**: the skin topples the body because nothing supports it in front of the centre of mass |
| stands | stands | refuted as stated: ANY sphere subset rescues it; the deficit is the skin's own load path (stiffness/friction), not location |
| falls | falls | refuted: neither subset suffices; the full 12 are needed, so count or distribution, not front/back |
| stands | falls | reversed; report it |

**Instrument checks (void if any fails).** (1) `skin_per_segment_carried` reproduces
−1.3 mm at 1.00 s and `skin_per_segment` reproduces the fall (−671.8 mm) to within 1 mm —
the engine edit must not move the existing arms. (2) The number of contact elements each
run reports must be 20 meshes + 0 / 4 / 8 / 12 foot spheres for
replace / posterior / anterior / carried, read out of the run's own snapshot, never the flag.
(3) The arm whose support the hypothesis says is sufficient must be able to fail: if
`anterior` stands, `posterior` is that check, and vice versa.

**Not separated by this run:** a stand with the anterior subset does not say the SKIN could
do it with a continuous carrier; it says a front support is sufficient and its absence is
what the per-piece skin lacks.

## Result (2026-10-30, second host, scored as pre-registered at the 2 s horizon)

`scripts/measure_segment_contact_meshes.py --steps {100,200}` with the five arms →
`out/skin_stance_support_100.json`, `out/skin_stance_support_200.json`. Passive plant: the
script calls `advance(dt)` with no excitations or commands, so the body is held only by its
joint stops and its contacts.

**Instrument checks.** (1) At 1.00 s `skin_per_segment_carried` reproduces **−1.30 mm** and
`skin_per_segment` reproduces **−671.78 mm**: the engine edit did not move the existing arms.
(2) Contact elements read from each run's own initial snapshot: **20 / 24 / 28 / 32** for
replace / posterior / anterior / carried, i.e. 20 meshes + 0 / 4 / 8 / 12 spheres, and 28 for
spheres-only. (3) Both subsets fail at 2 s, so each could fail.

| arm | contact elements | drift at 1.00 s | drift at 2.00 s | outcome (2 s) |
|---|---:|---:|---:|---|
| `spheres` (source feet, no skin) | 28 | 0.00 mm | −49.66 mm | sagging |
| `skin_per_segment_carried` (skin + all 12) | 32 | −1.30 mm | **−743.76 mm** | **falls** |
| `skin_per_segment` (skin only) | 20 | −671.78 mm | −789.81 mm | falls |
| `skin_per_segment_posterior` (skin + 4 heel/rearfoot) | 24 | −578.14 mm | −806.07 mm | falls |
| `skin_per_segment_anterior` (skin + 8 midfoot/toe) | 28 | −25.15 mm | integrator failed at t = 1.0835 s | falls |

**By the pre-registered table: posterior falls, anterior falls → "refuted: neither subset
suffices".** But the result does not stop there, because the reference it was built on did not
survive the longer horizon either: **the carried arm, the "identical geometry with all 12 spheres
that stands (−1.3 mm)", is 744 mm down at 2 s**, and even the source-sphere feet with no skin at all
are sagging 50 mm. The "stand" that motivated this experiment was a 1 s transient. So the question
"where must the support be for the skin to stand?" is not answerable on a passive plant: nothing
stands passively for 2 s, and what the arms differ in is *how fast* they go. In that ordering the
anterior subset is the slowest-falling skin arm at 1 s (−25 mm against −578 mm for posterior), which
is consistent with the front-support hypothesis but is a rate, not the stand/fall outcome that was
pre-registered, and it is not claimed as a confirmation.

What it does settle: (a) the 1 s numbers in `SEGMENT_CONTACT_SURFACES.md` are transients, not
equilibria, for every arm including the reference; (b) a stance comparison of contact models needs
a controller in the loop (as the cortical/LQR stance runs have), or a horizon long enough to see the
passive fall in every arm, and should report time-to-fall rather than a 1 s drift.
