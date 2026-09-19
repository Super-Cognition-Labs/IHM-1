# motor recruitment and reflexes past the ankle

Closes the first two items of `WORKBENCH_AUTHENTICITY.md` §3.1 as far as the
evidence reaches, and says exactly where it stops. Code: `ihm/assembly/recruitment.py`,
`ihm/assembly/reflexes.py` (Geyer & Herr section), `ihm/assembly/sensorimotor.py`.
Verification: `scripts/verify_recruitment.py` (53 checks, no native session).

## What changed, in one paragraph

`SensorimotorController` now takes two keyword options. `decoder='recruitment'`
replaces the linear regional-rate gain with a size-principle motor-unit pool.
`reflexes='geyer_herr_2010'` adds the source's knee and hip reflexes (VAS, HAM, GLU,
HFL) to the eight ankle primitives. Both default to `None`, and `None`/`None` runs
the original code path: `verify_recruitment.py` loads the original file from commit
`ab3d2b8` and requires every frame field and the checkpoint to be identical over a
sequence with descending drive, motor blocks and sensory blocks (30/30 frames
identical; the same comparison flags 30/30 frames when the gain is changed by 20%, so
it can fail). Only `model_sha256` differs, because it hashes the file. Nothing
reachable from the workbench changes: `embodied.py` constructs the controller with
defaults and is not touched.

## Sources, and a correction to the audit

| what | source | where in it | pinned |
|---|---|---|---|
| reflex laws, all groups | Geyer & Herr 2010, IEEE TNSRE 18(3):263–273, doi:10.1109/TNSRE.2010.2047592 | Table I p. 269; Appendix I p. 270 | `data/raw/sensorimotor/geyer_herr_2010.pdf`, `SOURCE_SHA256` in `sensorimotor.py` |
| motor-unit pool | Potvin & Fuglevand 2017, PLoS Comput Biol 13(6):e1005581 (CC-BY), restating Fuglevand, Winter & Patla 1993, J Neurophysiol 70:2470 | Methods pp. 19–20, Eqs 1–7; Fig 1B | `data/raw/sensorimotor/potvin_fuglevand_2017.pdf`, `POTVIN_FUGLEVAND_SHA256` in `recruitment.py` |

**The eight ankle constants are Geyer & Herr 2010, not Song & Geyer.** §3.1 and the
task that produced this file both say "Song & Geyer". The retained PDF is fetched
from a URL named `song.pdf`, and the hash-pinned file is the 2010 IEEE TNSRE paper;
`1.1`, `.71`, `1.2`, `.3` are its Table I. The attribution came from the filename.
`docs/research/SENSORIMOTOR_CONNECTIVITY.md` already had it right.

**Fuglevand 1993 itself is not retained**: the publisher returned 403 and a mirror
returned a challenge page. Every pool number is read from the 2017 restatement, whose
Methods give the rested pool completely. Nothing in either model is reconstructed from
memory; every reflex value was read off the page rendered at 300 dpi.

## The recruitment model

A deterministic rested pool of 120 units, identical for every muscle:

- twitch force `P(i) = 100^((i-1)/119)` and recruitment threshold `RTE(i) = 50^((i-1)/119)`,
  so weak units are recruited first (Henneman's size principle, Eqs 1–2);
- rate `R = 1·(E − RTE) + 8 imp/s` once recruited, saturating at a peak falling
  linearly from 35 to 25 imp/s (onion skin, Eq 3); `Emax = 50 + (25−8)/1 = 67`;
- contraction time 90 → 30 ms; force per unit `NF(R·CT)·P`, with `NF = 0.3·NR` up to
  NR = 0.4 and `1 − exp(−2 NR³)` above (Eqs 4–7).

**Output is excitation, nothing more.** The pool's rate-averaged isometric force as a
fraction of its own maximum, `F(E)/F(Emax)`, is what goes to the plant, and native
mechanics owns activation dynamics, fibre dynamics and force. Every frame says so
(`excitation_owner`, `activation_owner: mechanical_plant`). No twitch time course is
simulated, and nothing in the pool is random, so there is no generator to draw.

**How reflex and descending input meet.** Geyer & Herr define a reflex's output as
muscle *stimulation*, which in a pool model is already an output. So the reflex
stimulation `S` is mapped back to the pool input that would produce it alone,
descending input is added *there*, and the excitation is `S` plus what the pool adds on
top. At zero descending input this returns `S` exactly (checked: bit-identical to the
rate decoder with its dial off, both reflex sets), so the pool never reshapes a source
law. Descending drive is recruited by size, and how much it adds depends on how much
the segmental input has already recruited. Descending 0.2 of Emax alone gives 0.132;
on top of a reflex stimulation of 0.5 it adds 0.268.

**What is still engineered.** The descending input's magnitude is
`request × clip(cortical_drive_per_hz × regional rate)`, the same factor the rate
decoder uses. The pool replaces the *gain dial* and the linear decode, not the
cortical scaling. The frame's `decoder_basis` says this.

### One parameter the source gives only as a range

Potvin & Fuglevand say contraction time is "an inverse function of twitch amplitude
(see [19])", 90 ms to 30 ms, and do not print the law. The code uses
`CT = 90 ms · P^(−ln3/ln100)`. It is checked, not assumed: it matches the five
labelled units in their Fig 1B (MU20/40/60/80/100 at ≈76/63/52/43.5/36 ms) to
0.48 ms at worst. The null alternative, CT linear in P, misses by 34 ms, so the
figure check can tell the two laws apart.

### What is not muscle-specific (unmeasured)

The source pool is "generally representative, but not definitive, characterizations of
any specific skeletal muscle" (p. 19), and it is used as such for every muscle.
`max_isometric_force_n` scales each unit's absolute force (the pool at Emax sums to
the muscle's Fmax) and does not change the normalised curve. Ordering units *across*
muscles by Fmax or `optimal_fiber_length_m` is not a size-principle statement any
retained source makes, so it is not done. Per-muscle unit counts and upper
recruitment limits are unmeasured here, and a per-muscle table of both would measure
them. This is a borrowed shape, stated as one.

## The reflexes past the ankle

All from Appendix I of the pinned paper, with gains from Table I:

| group | plant members (declared) | stance | swing |
|---|---|---|---|
| VAS | vasint, vaslat, vasmed | `0.09 + 1.15 F − 2[φk−2.97]·[φk>2.97][φ̇k>0] − 1.2|F_contra|·DSup` | `0.01` |
| HAM | semimem, semiten, bflh | `0.05 + {1.91(θ−0.105) + 0.25θ̇}₊ · 1.2|F_ipsi|` | `0.01 + 0.65 F` |
| GLU | glmax1–3 | `0.05 + {0.68·1.91(θ−0.105) + 0.25θ̇}₊ · 1.2|F_ipsi| − 0.25 DSup` | `0.01 + 0.4 F` |
| HFL | iliacus, psoas | `0.05 + {…}₋ · 1.2|F_ipsi| + 0.25 DSup` (see reading below) | `0.01 + 0.35(l−0.6) − 4(l_HAM−0.85) + {1.15(θ−0.105)}_takeoff` |

Stimulations are limited to 0.01–1. θ is trunk forward lean, from the native torso and
pelvis transforms, and θ̇ is from the torso angular velocity. φk = π − knee_angle, so
180° is straight. |F_leg| is the native foot normal force over body weight. DSup
marks the trailing leg of a double support, taken from touchdown order. The swing HFL
lean term uses θ at that leg's previous take-off.

**Declared readings, recorded rather than silently resolved** (`GEYER_HERR_2010_SOURCE_NOTES`):

1. *Stance HFL sign.* Appendix I prints `+{PD}₋`, and `{}₋` is defined as "only
   negative values". Read literally, the trunk term can only *lower* hip-flexor
   stimulation. That contradicts the same paper on p. 265 (`S_GLU/HFL ∼ ±[kp(θ−θref)+kd θ̇]`,
   HFL taking the minus sign) and Fig. 1(d). The p. 265 reading is used:
   HFL gets `max(0, −PD)`. This is resolved within the source, not by analogy.
2. *k_bw.* Table I prints 1.2 with a tolerance of "1.3 … 5.0", which excludes its own
   value. The value column (1.2) is used.
3. *Take-off lean before the first take-off* is undefined and contributes 0.

**Declared engineering reductions (UNMEASURED).** The source has one muscle–tendon
unit per group and the plant splits it. Group force = Σ tendon force / Σ Fmax (the
rule the retained GAS law already used); group length = the Fmax-weighted mean of
member l/l_opt; every member gets its group's stimulation. Membership follows the
source's articulation. Excluded, so the choice can be argued with: bfsh (HAM is
biarticular), glmed/glmin (abductors), recfem/sart/tfl (not the source's
monoarticular HFL), addmag* (a hip extensor by moment arm, but an adductor the source
never names). Scoring member-resolved laws against lumped ones on EMG would measure
these reductions.

**Delays are not honoured per group.** The source uses 20 ms (SOL, TA, GAS), 10 ms
(VAS) and 5 ms (hip). The controller has one afferent + efferent split (10 + 10 ms)
plus up to one exchange interval, so every group runs at the ankle's 20 ms. Each frame
reports both (`reflex_source.source_delay_s`, `realised_loop_delay_s`).

**Blocks.** A missing *own* input silences a group's law, as the retained GAS rule
did when either head was absent. A missing *cross-group* input (soleus in TA, HAM
length in swing HFL) drops only that term, as the retained TA law did. Posture signals
(trunk, knee, foot load) cannot be blocked by muscle sensory blocks.

## What the verification establishes, and what it does not

It establishes that the pool matches its source, that the default is the original
code bit for bit, that every mode is idempotent and restores exactly, that the
extended set reproduces the eight ankle primitives exactly, that the pool is inert at
zero descending drive, and that every new law has the sign its source equation gives.

It does not establish that any of this improves anything. That needs a measurement
against data, which is the next section.

## Comparison against EMG

### Pre-registration (committed before the first scoring run)

Written after a `--dry-run` of the harness, which reads no EMG and scores nothing. It
exercised the replay plant only: the DGF constants match the on-disk opensim-core
header; mean vertical GRF / model weight = 0.983; harness fibre lengths stay in
0.52–1.26 l_opt; trunk pitch runs 0.077–0.139 rad, centred on the source's
θ_ref = 0.105; and **the rate decoder's dial moves only between 1.005 and 1.103
over the trial.**

**Data.** One measured walking trial on the model the plant is built from
(opensim-core `example3DWalking`: `subject_walk_scaled.osim`, IK `coordinates.sto`,
force-plate `grf_walk.mot`, right-leg EMG envelope `electromyography.sto`, all from
one trial; sha256s are written into the report). There is one subject and about one
stride. Nothing here generalises beyond that, and the report must say so.

**Replay plant, declared.** The native engine has no prescribed-motion interface, and
its static evaluator gives zero-speed fibre length and no tendon force, so the
measurement cannot run on the native plant. Kinematics are prescribed from the
trial. Muscles use the configuration opensim-core's own MocoInverse example uses for
this trial: DeGroote–Fregly 2016 curves, rigid tendon, no passive fibre force,
active width ×1.5, fibre damping 0.01, activation 10/40 ms. The loop closes at the
muscle (excitation → activation → force → afference). Trunk pitch is planar, from
pelvis_tilt + lumbar_extension. This is a kinematically prescribed replay, not
closed-loop locomotion.

**Arms**, each a fresh controller at 5 ms exchanges, with no descending drive (none
exists for this task): B1 rate+ankle (today), B2 rate+Geyer–Herr, B3 pool+ankle,
B4 pool+Geyer–Herr. **Trivial baselines:** A0 zero; A1 train-block mean; A2
*contact phase*: an affine map of the right foot's stance indicator (GRF > 5 N).
A2 is the one that matters, because a stance muscle's EMG is largely "on in stance".
**Control arm:** C1 = B2's prediction circularly shifted by half the window.

**Channels → muscles** (members share one group stimulation): soleus→soleus_r;
gastrocnemius→gasmed_r, gaslat_r; tibialis_anterior→tibant_r;
medial_hamstrings→semimem_r, semiten_r; biceps_femoris→bflh_r;
vastus_lateralis→vaslat_r; vastus_medius→vasmed_r; gluteus_maximus→glmax1–3_r. No
source law: rectus_femoris→recfem_r, gluteus_medius→glmed1–3_r.

**Scoring.** Window 0.55–1.795 s (the first 100 ms are excluded as warm-up). The
prediction at time t is the excitation applied from t. Pearson r and amplitude ratio
are computed with no fit. Skill = 1 − SSE/SSE_baseline after a per-channel affine map,
because EMG is normalised to its own maximum and its scale is arbitrary. The map is
fitted on alternate 100 ms blocks and scored on the others, then the folds swap, and
the script asserts that every scored sample is predicted by a map fitted without it.
Skill is reported against A1 (train mean) and against A0 (zero). Uncertainty comes
from 2,000 bootstrap resamples of the **blocks** (the items, about 13), paired between
arms, with one generator (seed 20260918) drawn in `main()` and asserted by its only
consumer.

**Verdict rules, fixed now:**

- **P1, reflexes past the ankle.** On the five knee/hip channels with a source law,
  the mean skill of B2 against A2 (contact phase) is compared. **PASS** if the paired
  95% CI of B2 − A2 lies above 0. Beating the train mean without beating A2 is
  reported as "carries no more than the contact phase", and P1 is still FAILED.
- **P2, decoder.** On all eight channels with a law, B4 − B2 mean skill is compared.
  "Measurable" if the paired 95% CI excludes 0. **Predicted before scoring:** the two
  differ only through the dial, which stays within 1.005–1.103 here, so the difference
  should be small. A null here is a statement about this task (a reflex replay with no
  descending drive), where the pool is inert by construction. It is not a statement
  about the pool.
- **P3, context, not a gate.** Ankle channels: B1 − A2.
- **Controls; if any fails, the run is void and is reported as void.** C1 must score
  below 0 on the knee/hip channels. The no-law channels must be constant in every
  reflex arm. B2 run twice must be bit-identical.

A gate that fails is recorded as FAILED. The instrument may change after a failure,
but these thresholds may not, and a changed instrument is a new run reported beside
this one.

### Result (run 1, script sha256 `e2f6a179…`, report `out/recruitment_comparison.json`)

250 samples scored in 13 blocks. **All controls pass:** C1 (half-window shift)
scores −0.059 on the knee/hip channels; the no-law channels are constant in every
reflex arm; and B2 run twice is bit-identical. The mean vertical GRF equals 0.983 of
the model's weight.

**P1: FAILED.** The knee/hip source reflexes do not beat the contact-phase baseline.
Mean skill against the train mean over the five knee/hip channels is B2 +0.023
(95% CI −0.094 to +0.084), against A2 contact +0.040; B2 − A2 CI is −0.124 to +0.126.

**P2: not distinguishable, as predicted.** B4 − B2 over the eight law channels has
CI −0.053 to +0.019 (B4 0.204, B2 0.218). With no descending drive, the pool is
inert and the arms differ only by the rate dial (1.005–1.103). Mean |Δexcitation| is
0.013 and the maximum is 0.485.

**P3, context: the eight ankle primitives carry real signal beyond phase.** B1 0.542
against A2 0.091; the B1 − A2 CI is +0.133 to +0.664.

Per channel, raw r (no fit), skill against the train mean, and amplitude ratio
rms(pred)/rms(EMG):

| channel | A2 contact skill | B1 rate+ankle r / skill / amp | B2 rate+GH r / skill / amp | B4 pool+GH r / skill / amp |
|---|---|---|---|---|
| soleus | +0.106 | +0.854 / +0.732 / 0.88 | +0.854 / +0.730 / 0.90 | +0.751 / +0.568 / 0.50 |
| gastrocnemius | +0.179 | +0.836 / +0.594 / 0.60 | +0.835 / +0.594 / 0.61 | +0.797 / +0.545 / 0.37 |
| tibialis anterior | −0.012 | +0.556 / +0.301 / 0.68 | +0.555 / +0.301 / 0.68 | +0.540 / +0.286 / 0.64 |
| medial hamstrings | −0.034 | const / 0 | **−0.253** / +0.065 / 0.24 | **−0.267** / +0.074 / 0.22 |
| biceps femoris | +0.019 | const / 0 | −0.037 / −0.097 / 0.26 | −0.050 / −0.084 / 0.25 |
| vastus lateralis | +0.180 | const / 0 | +0.395 / +0.128 / 1.43 | +0.448 / +0.178 / 1.08 |
| vastus medius | +0.082 | const / 0 | +0.203 / +0.017 / 1.21 | +0.279 / +0.059 / 0.92 |
| gluteus maximus | −0.048 | const / 0 | +0.205 / +0.003 / 0.12 | +0.209 / +0.003 / 0.12 |

**What the table says that the verdicts don't.**

- **Medial hamstrings' positive skill is earned by an inverted predictor.** Both
  folds fit a negative slope (−1.85, −2.19) and raw r is −0.25. The source's HAM law
  (trunk PD in stance, force feedback in swing) predicts this EMG *backwards*, and an
  affine map rescues the number by flipping the sign. **Skill after an affine fit can
  reward a law that is wrong in sign; read the slope and the raw r beside it.** For
  biceps femoris the fitted slope changes sign between folds (+0.93, −1.65), so
  there is no signal.
- **VAS has the right sign and nothing beyond phase.** Its r of 0.40–0.45 and skill of
  0.13–0.18 match what the stance indicator alone gets (0.18). GLU correlates at 0.21,
  but at 0.12× amplitude it has no skill.
- **The ankle result is the strongest evidence here that the existing primitives are
  doing something.** Soleus r = 0.85 and gastrocnemius r = 0.84, against 0.41 for
  contact phase. This holds on one trial, in a replay, with harness muscles.

**Post-hoc, NOT a result, and not tested.** The biggest decoder difference is on the
ankle, and it favours the rate decoder: soleus skill is 0.73 against 0.57, amplitude
0.88 against 0.50. It fits a mechanism that needs no cortex. The dial multiplies the
soleus F+ loop (gain 1.2), and a 0.5–10% change in the gain of a positive-feedback
loop moves its output a lot. If so, the "dial" is acting as a loop-gain adjustment,
not as descending control. This was noticed after the data, and the pre-registered
P2 set, which pools all eight channels, does not separate it. It needs a
pre-registered test on a **different trial**, with this one excluded, before anyone
quotes it. The retained Rajagopal pipeline has running EMG and GRF but no IK
coordinates on this model, so there is no second trial ready.

**What this does and does not show.** The recruitment model is implemented from its
source and cannot be exercised by this task: a reflex replay has no descending drive,
and the pool is exactly inert without one. **A task with a measured descending
command is what would test it**; none is retained here. The knee/hip reflexes are
transcribed faithfully and do not predict their muscles' EMG beyond gait phase on
this trial. One is predicted in the wrong sign. That is a negative result about the
source's hip laws on this subject in this replay, or about the lumping and
membership reductions declared above. This run cannot tell those apart. A
member-resolved variant (semimem, semiten and bflh each on its own length and
force) would separate them.
