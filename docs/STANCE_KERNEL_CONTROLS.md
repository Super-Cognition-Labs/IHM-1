# Does the learned kernel matter to stance? Eight kernels, twenty-four pushes

`docs/IBM_CURRICULUM16_KERNEL.md` closed on one trial per direction: the curriculum16
kernel recovered a (+5, 0, +5) N push to 1.7469 mm peak COM and a site-permuted
control recovered the same push to 1.7452 mm. That is n = 1 against n = 1 on one
push, and it was read — correctly, as far as it went — as "the dynamics are
load-bearing, the learning is not shown to be". This experiment turns that single
comparison into a measured effect with a stated uncertainty, and adds the control
that was missing: a kernel that was never trained at all.

Engineering sensorimotor simulation. Nothing here is a claim about a biological brain.

## PRE-REGISTRATION (written and committed before any arm below runs, 2026-10-29)

### Policies — one pipeline, eight kernels

Every policy is built by the pipeline in `IBM_CURRICULUM16_KERNEL.md` §Reproduce, line
for line (same teacher `f8740ca0…`, same five trajectories, 1,024 sites, signed-identity
encoder, 100 epochs, lr 0.001, reserve weighting, the same local readout fit and the same
×100 small-signal scaling). Only `dyn.embed` differs.

| id | kernel | learned content |
|---|---|---|
| `fused` | `ibm1_implicit.pt` (existing bundle `ibm_cortical_stance_v1`) | learned (34 sources) |
| `c16` | curriculum16 (existing `ibm_curriculum16_stance_v1`) | learned (16 objectives) |
| `perm-20260908` | curriculum16, rows permuted (existing `ibm_curriculum16_shuffled_stance_v1`) | positions destroyed |
| `perm-1`, `perm-2` | curriculum16, rows permuted, seeds 1 and 2 (new) | positions destroyed |
| `rand-11`, `rand-12`, `rand-13` | i.i.d. Gaussian directions at each row's original norm, seeds 11–13 (new) | everything but per-site norm destroyed |

New kernels come from `scripts/build_control_kernel.py` (derived from the retained
curriculum16 bundle, whose digest is recorded as the parent).

### Conditions — twenty-four pelvis-COM pushes, 100 ms at t = 1 s, 10 s horizon

`scripts/evaluate_native_cortical_stance.py`, registration
`data/models/engineering_stance_v1/registration.json`, 77.6122029 kg. Forces (x, 0, z) N:

| set | forces | n |
|---|---|---:|
| A: axis, 5 N | (±5, 0), (0, ±5) | 4 |
| B: diagonal, as the earlier runs | (±5, ±5) | 4 |
| C: 10 N, eight directions | 10·(cos θ, sin θ), θ = 0°, 45°, …, 315° | 8 |
| D: 20 N, eight directions | 20·(cos θ, sin θ) | 8 |

The policies were trained on axis pushes plus one declared extra, with the negative
diagonal held out entirely, so B–D are mostly outside the training set. The native
plant is deterministic, so the only randomness is the kernel seed; "seeds" in this
experiment are control-kernel seeds, and their spread is reported.

### Measurements, per policy × condition (full arm)

`peak_com_displacement_m`, `final_com_displacement_m`, `completed_horizon` (did not fall
before 10 s). The severed arm is run as the evaluator always runs it and reported, but
it is not part of the kernel comparison (it never depends on the kernel).

### The statistic and the bar — fixed now

Per condition c: **Δ_c = mean over {fused, c16} of ln(peak) − mean over the six controls
of ln(peak)**. Reported: mean Δ over the conditions where every one of the eight full arms
completed the horizon, with a 95% percentile bootstrap CI (10,000 resamples of
conditions, seed 0). A fall in any arm is reported separately as a recovery count per
policy, not folded into Δ.

| outcome | reading |
|---|---|
| mean Δ ≤ −0.05 and CI upper < 0 | the learned content **helps** stance by ≥ 5% peak excursion |
| CI inside (−0.05, +0.05) | **equivalence**: learned and destroyed kernels are indistinguishable at a 5% margin |
| mean Δ ≥ +0.05 and CI lower > 0 | the learned kernels are **worse** |
| anything else | inconclusive at this n; say so |

Recovery counts: learned content "matters for recovery" only if the two learned kernels
recover strictly more conditions than every control kernel does.

**An imported constant, declared as such:** the 5% margin is a judgement of what would be
worth reporting, not a measured noise floor. What measures the floor is the spread among
the six controls themselves (reported per condition as the control range); if that range
exceeds 5% the margin is too tight to resolve and the result will say so.

### Instrument checks (the run is void if any fails)

1. Re-running `c16` and `perm-20260908` on (+5, +5) reproduces 1.7469 mm and 1.7452 mm to
   1e-6 m (the deterministic native plant must give the earlier numbers back).
2. Each policy's `artifact_sha256` is read from its own report, and every new artifact's
   `brain_checkpoint_sha256` matches its kernel's `kernel_sha256`.
3. Each run's force is read back from its own report's `perturbation_force_n`.
4. The random controls preserve per-site norms (checked by the builder) — and are not
   accidentally equal to the parent (cosine of matched rows reported; expected ≈ 0).

### What it cannot say, whatever it returns

It tests one body, one mass, one policy family whose decoder is solved offline onto an
engineered LQR gain. Equivalence here would mean *this pipeline does not use kernel
content*, not that no pipeline could.

## Status: interrupted 2026-10-29, partial and NOT scored

Stopped on request (gb10-direct reassigned) after the build phase and 33 of the 72
existing-policy evaluations. No verdict is drawn from a partial run.

- **Built** (on gb10-direct, IBM-1 pinned at donor commit `267ad26`, equation source
  `9329ce…`, so the controls use the same equations as the existing policies):
  `perm-1`, `perm-2`, `rand-11`, `rand-12`, `rand-13` and the instrument rebuild
  `perm-20260908-rebuilt`, under `data/runtime/motor-learning/stance-kernel-controls/policies/`.
- **Instrument check 1, half met:** this host's native plant reproduces `c16` on
  (+5, +5) at **1.7469 mm** and `fused` at **1.7844 mm**, the earlier numbers to the
  printed digit. The `perm-20260908` half never ran before the stop.
- **Evaluated:** `fused` 24/24 conditions (recovers 15), `c16` 9/24 (recovers 7). Nothing
  else. Rows: `out/stance_kernel_controls_partial.json`.
- **Resume:** `scripts/run_stance_kernel_controls.py --phase eval --jobs 12` (skips done
  conditions), then `--phase collect`. The IBM-1 checkout used for training must stay at
  `267ad26`; evaluation does not read IBM-1.
