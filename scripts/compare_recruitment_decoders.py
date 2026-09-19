"""Pre-registered comparison: recruitment vs rate-gain decoder, ankle vs Geyer-Herr reflexes, against EMG.

The pre-registration is docs/MOTOR_RECRUITMENT.md, "Comparison against EMG", committed
before this script was first run in scoring mode. Nothing below may be changed to
suit a result; a changed instrument is a new, separately reported run.

WHAT IS MEASURED. One measured walking trial on the model the plant is built from
(opensim-core example3DWalking: Rajagopal-derived `subject_walk_scaled.osim`, IK
coordinates, force-plate GRF, and 10 channels of right-leg EMG envelope, all from
the same trial). Each arm is a fresh `SensorimotorController` driven through that
trial at 5 ms exchanges. Its excitations are compared with the EMG.

THE REPLAY PLANT IS NOT THE NATIVE PLANT, and this is why. The native engine has
no prescribed-motion interface, and its static evaluator returns fibre length at
zero speed with no tendon force, so it cannot give the force afference the source
laws need along a measured trajectory. The kinematics are prescribed from
measurement, and muscle states come from the configuration opensim-core's OWN
MocoInverse example uses for this very trial (exampleMocoInverse.cpp):
DeGrooteFregly2016 curves (constants copied from the on-disk opensim-core
source), rigid tendon (ModOpIgnoreTendonCompliance), passive fibre force ignored
(ModOpIgnorePassiveFiberForcesDGF), active force-length width x1.5
(ModOpScaleActiveFiberForceCurveWidthDGF(1.5)), fibre damping 0.01 and Millard
activation constants 10/40 ms carried over by ModOpReplaceMusclesWithDeGroote-
Fregly2016. Path lengths use the model's FunctionBasedPath polynomials, which are
the paths the native plant uses too. The loop closes at the muscle: controller
excitation -> harness activation -> harness tendon force -> afference.
Trunk pitch is planar: torso = Rz(pelvis_tilt + lumbar_extension). List and
rotation are ignored.

SPLIT. The only fitted quantities are a per-channel affine map (EMG is normalised
to its own maximum, so its scale is arbitrary) and the train-mean baseline. They are
fitted on alternate 100 ms blocks and scored on the others, then the folds swap, so
every scored sample is predicted by a map fitted without it (asserted). Pearson r
and amplitude ratio are reported without any fit. Uncertainty comes from bootstrapping
BLOCKS (the items), paired between arms, with a generator drawn once in main().
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.sensorimotor import SensorimotorController
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog
from ihm.native.gait_reference import _read
from ihm.native.moment_arm_control import FittedMomentArms
import xml.etree.ElementTree as ET

TRIAL = ROOT / 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
DGF_SOURCE = ROOT / 'data/raw/mechanics/opensim-core/OpenSim/Actuators/DeGrooteFregly2016Muscle.h'
SEED = 20260918
DT = .005
T0, T_END = .45, 1.795
SCORE_FROM = .55          # first 100 ms excluded: 20 ms reflex delay + controller warm-up
BLOCK_S = .1
BOOTSTRAP = 2000
G = 9.80665
ARMS = {'B1_rate_ankle': (None, None), 'B2_rate_geyer': (None, 'geyer_herr_2010'),
        'B3_pool_ankle': ('recruitment', None), 'B4_pool_geyer': ('recruitment', 'geyer_herr_2010')}
CHANNELS = {'soleus': ('soleus_r',), 'gastrocnemius': ('gasmed_r', 'gaslat_r'), 'tibialis_anterior': ('tibant_r',),
            'medial_hamstrings': ('semimem_r', 'semiten_r'), 'biceps_femoris': ('bflh_r',),
            'vastus_lateralis': ('vaslat_r',), 'vastus_medius': ('vasmed_r',), 'rectus_femoris': ('recfem_r',),
            'gluteus_maximus': ('glmax1_r', 'glmax2_r', 'glmax3_r'), 'gluteus_medius': ('glmed1_r', 'glmed2_r', 'glmed3_r')}
ANKLE_LAW = ('soleus', 'gastrocnemius', 'tibialis_anterior')
KNEE_HIP_LAW = ('medial_hamstrings', 'biceps_femoris', 'vastus_lateralis', 'vastus_medius', 'gluteus_maximus')
NO_LAW = ('rectus_femoris', 'gluteus_medius')

# --- DeGrooteFregly2016Muscle, copied from the on-disk opensim-core header ---------
B = ((0.8150671134243542, 1.055033428970575, 0.162384573599574, 0.063303448465465),
     (0.433004984392647, 0.716775413397760, -0.029947116970696, 0.200356847296188),
     (0.1, 1.0, 0.353553390593274, 0.0))
D1, D2, D3, D4 = -0.3211346127989808, -8.149, -0.374, 0.8825327733249912
WIDTH_SCALE = 1.5            # exampleMocoInverse: ModOpScaleActiveFiberForceCurveWidthDGF(1.5)
TAU_ACT, TAU_DEACT = .010, .040   # Millard2012EquilibriumMuscle defaults, carried over
SMOOTHING = .1               # DGF activation_dynamics_smoothing default
V_MAX = 10.                  # Muscle max_contraction_velocity default, l_opt/s
FIBER_DAMPING = .01          # this model's Millard fiber_damping, carried over


def dgf_constants_match_source():
    text = DGF_SOURCE.read_text()
    for value in [v for row in B for v in row] + [D1, D2, D3, D4]:
        if repr(value) not in text and f'{value:.15f}'.rstrip('0') not in text:
            return False
    return True


def f_active(l):
    x = (l - 1.) / WIDTH_SCALE + 1.
    return sum(b1 * np.exp(-.5 * (x - b2) ** 2 / (b3 + b4 * x) ** 2) for b1, b2, b3, b4 in B)


def f_velocity(v):
    t = D2 * v + D3
    return D1 * np.log(t + np.sqrt(t * t + 1.)) + D4


def activation_step(a, e, h):
    """DGF activation ODE (DeGrooteFregly2016Muscle.cpp), explicit Euler at 1 ms substeps."""
    for _ in range(int(round(h / .001))):
        z = .5 + 1.5 * a
        f = .5 * np.tanh(SMOOTHING * (e - a))
        a = a + .001 * ((f + .5) / (TAU_ACT * z) + (-f + .5) * z / TAU_DEACT) * (e - a)
    return np.clip(a, 0., 1.)


def unique_times(table):
    """Drop repeated time rows, refusing any whose data disagree (as gait_reference._read does)."""
    keep = np.r_[True, np.diff(table[:, 0]) > 0]
    for i in np.flatnonzero(~keep):
        if not np.array_equal(table[i], table[i - 1]):
            raise ValueError('Conflicting duplicate source time')
    table = table[keep]
    if not (np.diff(table[:, 0]) > 0).all() or not np.isfinite(table).all():
        raise ValueError('Source time must increase strictly')
    return table


class Trial:
    def __init__(self):
        names, q, qinfo = _read(TRIAL / 'coordinates.sto')
        self.coordinate = {path.split('/')[-2]: i for i, path in enumerate(names) if path.endswith('/value')}
        self.q = q
        self.qd = np.gradient(q, q[:, 0], axis=0)
        self.qd[:, 0] = q[:, 0]   # gradient of the time column is 1; keep it as the time axis
        grf_lines = (TRIAL / 'grf_walk.mot').read_text().splitlines()
        h = grf_lines.index('endheader')
        self.grf_names = grf_lines[h + 1].split()
        self.grf = unique_times(np.array([list(map(float, l.split())) for l in grf_lines[h + 2:] if l.strip()]))
        emg_lines = (TRIAL / 'electromyography.sto').read_text().splitlines()
        h = emg_lines.index('endheader')
        self.emg_names = emg_lines[h + 1].split()
        self._emg = unique_times(np.array([list(map(float, l.split())) for l in emg_lines[h + 2:] if l.strip()]))
        model = ET.parse(TRIAL / 'subject_walk_scaled.osim').getroot()
        self.mass_kg = sum(float(b.findtext('mass')) for b in model.iter('Body'))
        self.params = {m.attrib['name']: dict(fmax=float(m.findtext('max_isometric_force')), lopt=float(m.findtext('optimal_fiber_length')),
                                              lslack=float(m.findtext('tendon_slack_length')), alpha=float(m.findtext('pennation_angle_at_optimal')))
                       for m in model.iter('Millard2012EquilibriumMuscle')}
        self.paths = FittedMomentArms(TRIAL / 'subject_walk_scaled_FunctionBasedPathSet.xml')
        # Musculotendon length at every coordinate sample, velocity by central difference.
        self.muscles = sorted(self.params)
        L = np.array([[self.paths.length_and_moment_arms(m, self._qdict(row))[0] for m in self.muscles] for row in q])
        self.L, self.Ld = L, np.gradient(L, q[:, 0], axis=0)
        self.sources = {p: hashlib.sha256((TRIAL / p).read_bytes()).hexdigest() for p in
                        ('coordinates.sto', 'grf_walk.mot', 'electromyography.sto', 'subject_walk_scaled.osim',
                         'subject_walk_scaled_FunctionBasedPathSet.xml', 'grf_walk.xml')}

    def _qdict(self, row):
        return {name: row[i] for name, i in self.coordinate.items()}

    def at(self, table, t):
        i = min(max(int(np.searchsorted(table[:, 0], t, side='right')) - 1, 0), len(table) - 2)
        w = (t - table[i, 0]) / (table[i + 1, 0] - table[i, 0])
        return table[i] + w * (table[i + 1] - table[i])

    def value(self, name, t, speed=False):
        return float(self.at(self.qd if speed else self.q, t)[self.coordinate[name]])

    def vertical_grf(self, t):
        row = self.at(self.grf, t)
        return {s: max(0., float(row[self.grf_names.index(f'ground_force_{s}_vy')])) for s in ('r', 'l')}

    def emg(self, t):
        return dict(zip(self.emg_names[1:], self.at(self._emg, t)[1:].tolist()))

    def muscle_states(self, t, activation):
        """Rigid-tendon DGF states at the measured pose; activation is the harness's own."""
        L = self.at(np.column_stack([self.q[:, 0], self.L]), t)[1:]
        Ld = self.at(np.column_stack([self.q[:, 0], self.Ld]), t)[1:]
        out = {}
        for k, m in enumerate(self.muscles):
            p = self.params[m]
            along = L[k] - p['lslack']
            height = p['lopt'] * math.sin(p['alpha'])
            lf = math.sqrt(along * along + height * height)
            cos = along / lf
            v = float(np.clip(Ld[k] * cos / (V_MAX * p['lopt']), -1., 1.))
            ln = lf / p['lopt']
            force = p['fmax'] * (activation[m] * f_active(ln) * f_velocity(v) + FIBER_DAMPING * v) * cos
            out[m] = {'fiber_length_m': lf, 'optimal_fiber_length_m': p['lopt'], 'tendon_force_n': max(0., float(force)),
                      'max_isometric_force_n': p['fmax'], 'normalized_fiber_velocity': v,
                      'sensor_basis': 'replay harness: measured kinematics, FunctionBasedPath length, rigid-tendon DGF (MocoInverse example3DWalking configuration); NOT native'}
        return out

    def observation(self, t, time_s, activation):
        tilt = self.value('pelvis_tilt', t) + self.value('lumbar_extension', t)
        rate = self.value('pelvis_tilt', t, True) + self.value('lumbar_extension', t, True)
        rz = lambda phi: [[math.cos(phi), -math.sin(phi), 0., 0.], [math.sin(phi), math.cos(phi), 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]]
        return {'time_s': time_s, 'muscles': self.muscle_states(t, activation), 'foot_contact_force_n': self.vertical_grf(t),
                'effective_native_body_mass_kg': self.mass_kg,
                'native_bodies': {'torso': {'transform_ground': rz(tilt), 'angular_velocity_rad_s': [0., 0., rate]},
                                  'pelvis': {'transform_ground': rz(self.value('pelvis_tilt', t)), 'angular_velocity_rad_s': [0., 0., 0.]}},
                'joints': {f'knee_angle_{s}': {'value': self.value(f'knee_angle_{s}', t), 'speed': self.value(f'knee_angle_{s}', t, True)} for s in ('r', 'l')}}


def run_arm(trial, catalog, decoder, reflexes):
    controller = SensorimotorController.from_root(ROOT, muscle_catalog=catalog, decoder=decoder, reflexes=reflexes)
    muscles = controller.muscles
    missing = set(muscles) - set(trial.params)
    if missing:
        raise ValueError(f'Catalog muscles absent from the trial model: {sorted(missing)}')
    activation = {m: .01 for m in muscles}
    applied = {m: 0. for m in muscles}
    steps = int(round((T_END - T0) / DT))
    times, excitation, gains, thetas, afferent = [], [], [], [], []
    for k in range(steps):
        t = T0 + k * DT
        obs = trial.observation(t, controller.time_s, activation)
        frame = controller.step(DT, obs)
        a = activation_step(np.array([activation[m] for m in muscles]), np.array([applied[m] for m in muscles]), DT)
        activation = dict(zip(muscles, a.tolist()))
        applied = dict(frame['motor_excitations'])
        times.append(t + DT)
        excitation.append([applied[m] for m in muscles])
        gains.append([frame['descending_gain'][m] for m in muscles])
        thetas.append(float('nan') if 'reflex_state' not in frame else (frame['reflex_state']['posture'] or {}).get('theta', float('nan')))
        afferent.append({m: (obs['muscles'][m]['fiber_length_m'] / obs['muscles'][m]['optimal_fiber_length_m'],
                             obs['muscles'][m]['tendon_force_n'] / obs['muscles'][m]['max_isometric_force_n'])
                         for m in ('soleus_r', 'tibant_r', 'vaslat_r', 'semimem_r', 'glmax1_r', 'iliacus_r')})
    return {'times': np.array(times), 'muscles': list(muscles), 'excitation': np.array(excitation), 'gain': np.array(gains),
            'theta': np.array(thetas, float), 'afferent': afferent, 'model_sha256': controller.model_sha256,
            'decoder': frame.get('decoder', 'rate_gain'), 'reflex_set': frame.get('reflex_set', 'ankle')}


def channel_prediction(run, channel):
    idx = [run['muscles'].index(m) for m in CHANNELS[channel]]
    return run['excitation'][:, idx].mean(axis=1)


def crossfit(pred, target, blocks):
    """Per-channel affine map fitted on alternate blocks, scored on the others; folds swap."""
    out = np.full_like(target, np.nan)
    fits = []
    for fold in (0, 1):
        train = blocks % 2 == fold
        test = ~train
        assert not np.any(train & test) and train.any() and test.any()
        if np.std(pred[train]) < 1e-12:
            coef = np.array([0., target[train].mean()])
        else:
            X = np.column_stack([pred[train], np.ones(train.sum())])
            coef = np.linalg.lstsq(X, target[train], rcond=None)[0]
        out[test] = coef[0] * pred[test] + coef[1]
        fits.append(coef.tolist())
    assert np.isfinite(out).all(), 'every scored sample must be predicted by a map fitted without it'
    return out, fits


def score(yhat, baseline_mean, target, index):
    """Skill against the cross-fitted train mean and against zero, on the given samples."""
    sse = np.sum((yhat[index] - target[index]) ** 2)
    return {'skill_vs_train_mean': float(1 - sse / np.sum((baseline_mean[index] - target[index]) ** 2)),
            'skill_vs_zero': float(1 - sse / np.sum(target[index] ** 2))}


def pearson(a, b):
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def bootstrap(rng, blocks, per_arm_channels, target, means, arms, channels):
    """Resample BLOCKS (items), not samples; returns per-resample mean skill per arm."""
    assert isinstance(rng, np.random.Generator), 'the generator must be drawn once in main() and passed in'
    ids = np.unique(blocks)
    members = {b: np.flatnonzero(blocks == b) for b in ids}
    draws = rng.integers(0, len(ids), size=(BOOTSTRAP, len(ids)))
    out = {arm: np.empty(BOOTSTRAP) for arm in arms}
    for r in range(BOOTSTRAP):
        index = np.concatenate([members[ids[j]] for j in draws[r]])
        for arm in arms:
            out[arm][r] = np.mean([score(per_arm_channels[arm][c], means[c], target[c], index)['skill_vs_train_mean'] for c in channels])
    return out, draws


def ci(x):
    return [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dry-run', action='store_true', help='exercise the replay plant only; reads no EMG, scores nothing')
    ap.add_argument('--out', default='out/recruitment_comparison.json')
    a = ap.parse_args()
    rng = np.random.default_rng(SEED)   # drawn ONCE; passed to the only consumer
    catalog = native_muscle_catalog(ROOT)
    trial = Trial()
    checks = {'dgf_constants_match_opensim_source': dgf_constants_match_source()}
    grf_total = np.array([sum(trial.vertical_grf(t).values()) for t in np.arange(T0, T_END, DT)])
    checks['mean_vertical_grf_over_model_weight'] = float(grf_total.mean() / (trial.mass_kg * G))
    print(json.dumps(checks), flush=True)
    if not checks['dgf_constants_match_opensim_source']:
        raise SystemExit('DGF constants differ from the on-disk opensim-core source')
    if a.dry_run:
        run = run_arm(trial, catalog, None, 'geyer_herr_2010')
        aff = run['afferent']
        for m in aff[0]:
            ln = [x[m][0] for x in aff]; fn = [x[m][1] for x in aff]
            print(f'{m:10s} l/lopt {min(ln):.3f}..{max(ln):.3f}   F/Fmax {min(fn):.3f}..{max(fn):.3f}', flush=True)
        th = run['theta'][np.isfinite(run['theta'])]
        print(f'trunk pitch {th.min():.3f}..{th.max():.3f} rad  (source theta_ref 0.105)', flush=True)
        g = run['gain']
        print(f'rate-decoder dial range {g.min():.4f}..{g.max():.4f}', flush=True)
        return

    runs = {name: run_arm(trial, catalog, *cfg) for name, cfg in ARMS.items()}
    again = run_arm(trial, catalog, *ARMS['B2_rate_geyer'])
    checks['idempotent_B2_rerun'] = bool(np.array_equal(again['excitation'], runs['B2_rate_geyer']['excitation']))
    times = runs['B1_rate_ankle']['times']
    assert all(np.array_equal(times, r['times']) for r in runs.values())
    window = (times >= SCORE_FROM - 1e-12) & (times <= T_END + 1e-12)
    t = times[window]
    blocks = np.floor((t - SCORE_FROM) / BLOCK_S + 1e-9).astype(int)
    stance = np.array([1. if trial.vertical_grf(x)['r'] > 5. else 0. for x in t])
    target = {c: np.array([trial.emg(x)[c] for x in t]) for c in CHANNELS}
    means = {c: crossfit(np.zeros_like(t), target[c], blocks)[0] for c in CHANNELS}
    everything = np.arange(len(t))

    predictions = {name: {c: channel_prediction(r, c)[window] for c in CHANNELS} for name, r in runs.items()}
    half = len(t) // 2
    predictions['C1_B2_shifted_half_window'] = {c: np.roll(predictions['B2_rate_geyer'][c], half) for c in CHANNELS}
    predictions['A2_contact'] = {c: stance for c in CHANNELS}
    crossfitted = {name: {} for name in predictions}
    table = {}
    for name, per in predictions.items():
        table[name] = {}
        for c in CHANNELS:
            yhat, fits = crossfit(per[c], target[c], blocks)
            crossfitted[name][c] = yhat
            row = score(yhat, means[c], target[c], everything)
            raw = per[c]
            row.update(r=pearson(raw, target[c]), amplitude_ratio=float(np.sqrt(np.mean(raw ** 2)) / np.sqrt(np.mean(target[c] ** 2))),
                       affine_fits=fits, prediction_constant=bool(np.std(raw) < 1e-12))
            table[name][c] = row
    crossfitted['A1_train_mean'] = means
    table['A1_train_mean'] = {c: score(means[c], means[c], target[c], everything) for c in CHANNELS}
    table['A0_zero'] = {c: score(np.zeros_like(t), means[c], target[c], everything) for c in CHANNELS}

    arms = list(ARMS) + ['A2_contact', 'C1_B2_shifted_half_window']
    boot_kh, draws = bootstrap(rng, blocks, crossfitted, target, means, arms, KNEE_HIP_LAW)
    boot_law, _ = bootstrap(rng, blocks, crossfitted, target, means, arms, ANKLE_LAW + KNEE_HIP_LAW)
    boot_ankle, _ = bootstrap(rng, blocks, crossfitted, target, means, arms, ANKLE_LAW)
    mean_skill = lambda name, chans: float(np.mean([table[name][c]['skill_vs_train_mean'] for c in chans]))

    P1 = {'B2_minus_A2_ci': ci(boot_kh['B2_rate_geyer'] - boot_kh['A2_contact']),
          'B2_ci': ci(boot_kh['B2_rate_geyer']), 'B2': mean_skill('B2_rate_geyer', KNEE_HIP_LAW),
          'A2': mean_skill('A2_contact', KNEE_HIP_LAW), 'B1': mean_skill('B1_rate_ankle', KNEE_HIP_LAW),
          'B2_minus_B1_ci': ci(boot_kh['B2_rate_geyer'] - boot_kh['B1_rate_ankle'])}
    P1['verdict'] = ('PASS: knee/hip source reflexes beat the contact-phase baseline' if P1['B2_minus_A2_ci'][0] > 0 else
                     'FAIL: knee/hip source reflexes do not beat the contact-phase baseline')
    diff = runs['B4_pool_geyer']['excitation'] - runs['B2_rate_geyer']['excitation']
    P2 = {'B4_minus_B2_ci': ci(boot_law['B4_pool_geyer'] - boot_law['B2_rate_geyer']),
          'B4': mean_skill('B4_pool_geyer', ANKLE_LAW + KNEE_HIP_LAW), 'B2': mean_skill('B2_rate_geyer', ANKLE_LAW + KNEE_HIP_LAW),
          'max_abs_excitation_difference': float(np.abs(diff).max()), 'mean_abs_excitation_difference': float(np.abs(diff).mean()),
          'rate_dial_range': [float(runs['B2_rate_geyer']['gain'].min()), float(runs['B2_rate_geyer']['gain'].max())]}
    lo, hi = P2['B4_minus_B2_ci']
    P2['verdict'] = ('decoder difference MEASURABLE on this task' if lo > 0 or hi < 0 else
                     'decoder difference NOT distinguishable on this task')
    P3 = {'B1_minus_A2_ci': ci(boot_ankle['B1_rate_ankle'] - boot_ankle['A2_contact']),
          'B1': mean_skill('B1_rate_ankle', ANKLE_LAW), 'A2': mean_skill('A2_contact', ANKLE_LAW)}
    controls = {'C1_shifted_skill_knee_hip': mean_skill('C1_B2_shifted_half_window', KNEE_HIP_LAW),
                'C1_pass': mean_skill('C1_B2_shifted_half_window', KNEE_HIP_LAW) < 0,
                'C2_no_law_channels_constant_in_every_reflex_arm': all(table[n][c]['prediction_constant'] for n in ARMS for c in NO_LAW),
                'C3_idempotent': checks['idempotent_B2_rerun'],
                'bootstrap_items': int(len(np.unique(blocks))), 'bootstrap_draws_shape': list(draws.shape)}
    report = {'schema': 'ihm.recruitment-comparison.v1', 'preregistration': 'docs/MOTOR_RECRUITMENT.md#comparison-against-emg',
              'seed': SEED, 'dt_s': DT, 'window_s': [SCORE_FROM, T_END], 'block_s': BLOCK_S, 'bootstrap': BOOTSTRAP,
              'samples_scored': int(len(t)), 'checks': checks, 'controls': controls, 'P1_reflexes_past_ankle': P1,
              'P2_decoder': P2, 'P3_ankle_context': P3, 'table': table, 'sources_sha256': trial.sources,
              'arm_model_sha256': {n: r['model_sha256'] for n, r in runs.items()},
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'replay_plant': 'measured kinematics + GRF; rigid-tendon DGF muscles (MocoInverse example3DWalking configuration); NOT the native plant'}
    out = ROOT / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x)) + '\n')
    print(json.dumps({k: report[k] for k in ('controls', 'P1_reflexes_past_ankle', 'P2_decoder', 'P3_ankle_context')}, indent=1), flush=True)
    for name in ['A0_zero', 'A1_train_mean', 'A2_contact', *ARMS, 'C1_B2_shifted_half_window']:
        row = '  '.join(f"{c[:8]}:{table[name][c]['skill_vs_train_mean']:+.3f}" for c in CHANNELS)
        print(f'{name:26s} {row}', flush=True)
    for name in ARMS:
        row = '  '.join(f"{c[:8]}:{'  n/a ' if table[name][c]['r'] is None else format(table[name][c]['r'], '+.3f')}" for c in CHANNELS)
        print(f'r {name:24s} {row}', flush=True)


if __name__ == '__main__':
    main()
