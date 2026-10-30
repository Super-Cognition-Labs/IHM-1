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

## Result (2026-10-30, all 216 arm-runs, scored as pre-registered)

Run on gb10-direct (second host), `scripts/run_stance_kernel_controls.py --phase eval --jobs 10`
then `--phase collect` → `out/stance_kernel_controls.json`.

**Instrument checks: all pass.** (1) `c16` reproduces **1.7469 mm** and `perm-20260908`
reproduces **1.7452 mm** on (+5, +5), to the printed digit; the rebuilt `perm-20260908` also gives
1.7452 mm, so the pipeline itself reproduces. (2) Each row carries its own `artifact_sha256`. (3) Every
run's force is read back from its own report (`force_n`). (4) Checked by the builder at build time.

**The pre-registered statistic: inconclusive at this n.** Only **3 of 24 conditions** had all eight
full arms complete the 10 s horizon (A: +5,0; −5,0; 0,+5), because the controls fall often. Over those
three, mean Δ ln(peak) = **−0.032**, 95% bootstrap CI **[−0.070, −0.003]** (10,000 resamples, seed 0).
That is neither ≤ −0.05 with an upper bound below 0, nor inside ±0.05: inconclusive. And the floor the
margin needed is not there: the median spread among the six controls is **0.071** in ln(peak), wider
than the 5% margin, so the margin is too tight to resolve on this pipeline, as the pre-registration
said it would report.

**Recovery (of 24 conditions; completed the horizon without a fall):**

| policy | kernel | A (4) | B (4) | C (8) | D (8) | total |
|---|---|---:|---:|---:|---:|---:|
| `fused` | learned | 4 | 4 | 7 | 0 | **15** |
| `c16` | learned | 4 | 3 | 3 | 0 | **10** |
| `perm-20260908` | permuted | 4 | 4 | 4 | 0 | 12 |
| `perm-1` | permuted | 3 | 3 | 3 | 0 | 9 |
| `perm-2` | permuted | 3 | 1 | 0 | 0 | 4 |
| `rand-11` | random | 4 | 3 | 2 | 0 | 9 |
| `rand-12` | random | 4 | 1 | 0 | 0 | 5 |
| `rand-13` | random | 4 | 3 | 3 | 0 | 10 |
| `perm-20260908-rebuilt` | instrument | 4 | 4 | 4 | 0 | 12 |

The pre-registered recovery rule ("both learned kernels recover strictly more conditions than every
control") is **not met**: `c16` recovers 10, `perm-20260908` 12. Nothing at 20 N (set D) is recovered
by any kernel. Descriptively, the learned pair averages 12.5 recoveries against 8.2 for the six
controls, and `fused` alone is the only policy that recovers most of the 10 N set (7/8) — a pattern
worth a larger, pre-registered n, not a finding.

**What changed from the n = 1 reading.** The earlier single comparison drew a control that happens to
be one of the strong ones. On the very push it used, (+5, +5), three of the six controls fall
(`perm-2` 524 mm, `rand-12` 373 mm, `rand-13` 162 mm) while both learned kernels hold. So "a permuted
kernel does as well" was true of that permutation and is not true of permutations in general; and
"the learned content helps" is not shown either. The defensible statement is: **kernel content
changes whether this pipeline recovers, the control-to-control spread is larger than any
learned-vs-control difference this design can resolve, and the learned kernels are not shown to be
better.**

Runs: `data/runtime/motor-learning/stance-kernel-controls/evals/` on gb10-direct (gitignored).
