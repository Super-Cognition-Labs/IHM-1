"""Verify the motor-unit pool, the Geyer & Herr 2010 reflex set, and that the default is untouched.

No native session, no optimiser. Every check here can fail, and several are built
to fail for the obvious wrong implementation (see `null` cases).

  1. Pool known answers from Potvin & Fuglevand 2017: Emax = 67, last unit
     recruited at 50/67 = 74.6% of Emax, Eq 5/6 continuity, contraction times
     against the five labelled points of Fig 1B -- and the NULL: a CT law linear
     in P must MISS Fig 1B, or the figure check cannot tell laws apart.
  2. Default = the original implementation, EXACTLY, on a sequence with
     descending drive, sensory and motor blocks (loads ORIGINAL_COMMIT's file).
  3. Idempotence: checkpoint -> step(x) -> restore -> step(x) is identical, per mode.
  4. Identities between modes: the extended reflex set reproduces the eight ankle
     primitives exactly; recruitment at zero descending drive equals the rate
     decoder with its dial disabled (the pool does not re-shape the source law).
  5. Sign tests for every new source law, from the equations in Appendix I.
"""
import importlib.util
import math
import subprocess
import sys
import tempfile
from copy import deepcopy
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly.recruitment import MotorUnitPool, normalized_force, POTVIN_FUGLEVAND_PATH, POTVIN_FUGLEVAND_SHA256
from ihm.assembly.reflexes import geyer_herr_2010_stimulation, GEYER_HERR_2010_TABLE_I as T
from ihm.assembly.sensorimotor import SensorimotorController, SensorimotorParameters
from ihm.assembly.sensorimotor_catalog import native_muscle_catalog
import hashlib

CATALOG = native_muscle_catalog(ROOT)
FAILED = []
# The last commit holding the ORIGINAL decoder. Not HEAD: once this change is
# committed, HEAD is the new file and the equivalence test would compare the code
# with itself -- a control that cannot fail.
ORIGINAL_COMMIT = 'ab3d2b81eb816284adc00445e3f4ecde0290ede0'



def check(name, condition, detail=''):
    print(('PASS ' if condition else 'FAIL ') + name + (f'  [{detail}]' if detail else ''), flush=True)
    if not condition:
        FAILED.append(name)


def rz(phi):
    c, s = math.cos(phi), math.sin(phi)
    return [[c, -s, 0., 0.], [s, c, 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]]


def observation(time=0., extension=1.0, load=100., *, loads=None, theta=0., theta_rate=0., knee=(.3, .3), knee_speed=(0., 0.),
                lengths=None, forces=None):
    """verify_sensorimotor's fixture plus the posture fields the hip/knee laws read."""
    muscles = {}
    for row in CATALOG:
        name = row['id']
        length = .1 * (extension if name == 'tibant_r' else 1.)
        force = load if name == 'soleus_r' else 0.
        if lengths and name in lengths: length = .1 * lengths[name]
        if forces and name in forces: force = forces[name]
        muscles[name] = {'fiber_length_m': length, 'optimal_fiber_length_m': .1, 'tendon_force_n': force,
                         'max_isometric_force_n': 1000., 'sensor_basis': 'synthetic unit fixture CE'}
    loads = loads or {'r': 100., 'l': 0.}
    return {'time_s': time, 'muscles': muscles, 'foot_contact_force_n': dict(loads),
            'effective_native_body_mass_kg': 75.,
            'native_bodies': {'torso': {'transform_ground': rz(-theta), 'angular_velocity_rad_s': [0., 0., -theta_rate]},
                              'pelvis': {'transform_ground': rz(0.), 'angular_velocity_rad_s': [0., 0., 0.]}},
            'joints': {'knee_angle_r': {'value': knee[0], 'speed': knee_speed[0]},
                       'knee_angle_l': {'value': knee[1], 'speed': knee_speed[1]}}}


def pool_known_answers():
    pool = MotorUnitPool()
    check('pool source pinned', hashlib.sha256((ROOT / POTVIN_FUGLEVAND_PATH).read_bytes()).hexdigest() == POTVIN_FUGLEVAND_SHA256)
    check('Emax = 67 (Potvin & Fuglevand p.20)', abs(pool.max_excitation - 67.) < 1e-12, pool.max_excitation)
    check('last unit recruited at 74.6% of Emax', abs(pool.threshold[-1] / pool.max_excitation - 50 / 67) < 1e-12)
    check('Eq 5/6 continuous at NR=0.4', abs(normalized_force(.4) - normalized_force(.4 + 1e-12)) < 2e-3)
    # Fig 1B, read off the rendered figure (every 20th unit is labelled).
    figure = {20: 76., 40: 63., 60: 52., 80: 43.5, 100: 36.}
    err = max(abs(1000 * pool.contraction_time_s[i - 1] - ct) for i, ct in figure.items())
    check('CT power law reproduces Fig 1B', err < 1.5, f'max |err| {err:.2f} ms')
    linear = {i: 90 - 60 * (pool.twitch_force[i - 1] - 1) / 99 for i in figure}
    err_linear = max(abs(linear[i] - ct) for i, ct in figure.items())
    check('NULL: CT linear in P misses Fig 1B (figure check can discriminate)', err_linear > 20, f'{err_linear:.1f} ms')
    check('CT endpoints 90 / 30 ms', abs(pool.contraction_time_s[0] - .09) < 1e-15 and abs(pool.contraction_time_s[-1] - .03) < 1e-12)
    x = np.linspace(0, 1, 2001)
    y = pool.output(x)
    check('pool output 0 at 0, 1 at Emax, nondecreasing', y[0] == 0 and abs(y[-1] - 1) < 1e-12 and np.all(np.diff(y) >= 0))
    check('units at Emax sum to Fmax', abs(pool.unit_forces_n(1., 1234.5).sum() - 1234.5) < 1e-9)
    check('size principle: thresholds and twitch forces rise together', np.all(np.diff(pool.threshold) > 0) and np.all(np.diff(pool.twitch_force) > 0))
    check('onion skin: peak rates fall with threshold', np.all(np.diff(pool.peak_rate_hz) < 0))
    s = np.linspace(0, 1, 101)
    e0, _, _ = pool.combine(s, np.zeros_like(s))
    check('combine at zero descending returns the source stimulation exactly', np.array_equal(e0, s))
    a = pool.combine(.3, .2)[0]; b = pool.combine(.3, .2)[0]
    check('combine idempotent', float(a) == float(b))
    d = np.linspace(0, 1, 51)
    check('combine monotone in descending input', np.all(np.diff(pool.combine(np.full_like(d, .3), d)[0]) >= 0))
    inv = pool.output(pool.input_for(s))
    check('output(input_for(y)) >= y, within one recruitment step', np.all(inv >= s - 1e-15) and (inv - s).max() < .005, f'{(inv-s).max():.4f}')
    return pool


def load_head_module():
    """The original sensorimotor.py (ORIGINAL_COMMIT), imported as a sibling module of ihm.assembly."""
    raw = subprocess.run(['git', 'show', f'{ORIGINAL_COMMIT}:ihm/assembly/sensorimotor.py'], cwd=ROOT, capture_output=True, check=True).stdout
    current = (ROOT / 'ihm/assembly/sensorimotor.py').read_bytes()
    check('equivalence reference is NOT the file under test', raw != current)
    tmp = Path(tempfile.mkdtemp(prefix='sensorimotor-head-'))
    (tmp / 'sensorimotor.py').write_bytes(raw)
    (tmp / 'brain.py').write_bytes(subprocess.run(['git', 'show', f'{ORIGINAL_COMMIT}:ihm/assembly/brain.py'], cwd=ROOT,
                                                  capture_output=True, check=True).stdout)
    spec = importlib.util.spec_from_file_location('ihm.assembly._sensorimotor_head', tmp / 'sensorimotor.py',
                                                  submodule_search_locations=None)
    module = importlib.util.module_from_spec(spec)
    module.__package__ = 'ihm.assembly'
    spec.loader.exec_module(module)
    return module


def legacy_sequence(controller):
    frames = []
    for i in range(30):
        kwargs = {}
        if i % 7 == 3: kwargs['motor_blocks'] = ['tibant_r']
        if 10 <= i < 14: kwargs['sensory_blocks'] = ['soleus_r', 'gasmed_l']
        if i % 5 == 0: kwargs['descending'] = {'vaslat_r': .4, 'soleus_l': .2, 'glmax1_r': 1.}
        obs = observation(controller.time_s, 1 + .02 * math.sin(i / 3), 50 + 40 * math.cos(i / 4),
                          loads={'r': 100. if i % 9 < 6 else 0., 'l': 0. if i % 9 < 3 else 100.})
        frames.append(controller.step(.01, obs, **kwargs))
    return frames


def default_equals_head():
    head = load_head_module()
    old = head.SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG)
    new = SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG)
    a, b = legacy_sequence(old), legacy_sequence(new)
    mismatches = 0
    for fa, fb in zip(a, b):
        fa = dict(fa); fb = dict(fb)
        fa.pop('model_sha256'); fb.pop('model_sha256')
        if fa != fb: mismatches += 1
    check('default decoder/reflexes == original implementation, every frame field', mismatches == 0, f'{mismatches} of {len(a)} frames differ')
    perturbed = legacy_sequence(head.SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG,
                                parameters=head.SensorimotorParameters(cortical_gain_per_hz=.006)))
    moved = sum(dict(fa, model_sha256=None) != dict(fb, model_sha256=None) for fa, fb in zip(perturbed, b))
    check('NULL: the equivalence comparison detects a 20% gain change', moved > 0, f'{moved} of {len(b)} frames differ')
    check('default frame declares the original scope', b[-1]['scope'].startswith('all catalog effectors') and 'eight ankle' in b[-1]['scope'])
    ck_old, ck_new = old.checkpoint(), new.checkpoint()
    ck_old.pop('model_sha256'); ck_new.pop('model_sha256')
    check('default checkpoint == original implementation', ck_old == ck_new)


def stepping(controller, n, start=0):
    frames = []
    for i in range(start, start + n):
        t = controller.time_s
        stance_r = (i // 8) % 2 == 0
        obs = observation(t, 1 + .05 * math.sin(i / 2), 300 + 200 * math.sin(i / 5),
                          loads={'r': 600. if stance_r else 0., 'l': 600. if (not stance_r or i % 8 < 2) else 0.},
                          theta=.105 + .1 * math.sin(i / 6), theta_rate=.3 * math.cos(i / 6),
                          knee=(.05 + .2 * (1 + math.sin(i / 4)), .3), knee_speed=(-1. if i % 4 < 2 else 1., 0.),
                          lengths={'iliacus_r': 1.1, 'psoas_r': 1.0, 'semimem_r': .9 + .1 * math.sin(i), 'semiten_r': 1., 'bflh_r': 1.},
                          forces={'vaslat_r': 2000. + 500 * math.sin(i), 'semimem_r': 400., 'glmax1_r': 300., 'gasmed_r': 800.})
        frames.append(controller.step(.01, obs, descending={'vaslat_r': .3, 'recfem_r': .5} if i % 3 == 0 else None))
    return frames


def strip(frame):
    frame = deepcopy(frame)
    frame.pop('brain', None)
    return frame


def idempotence_and_restore():
    for decoder in (None, 'recruitment'):
        for reflexes in (None, 'geyer_herr_2010'):
            c = SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, decoder=decoder, reflexes=reflexes)
            stepping(c, 25)
            saved = c.checkpoint()
            first = stepping(c, 12, 25)
            c.restore(saved)
            second = stepping(c, 12, 25)
            check(f'restore -> same steps identical ({decoder},{reflexes})', first == second)
            c.restore(saved)
            x = observation(c.time_s, 1.03, 250., loads={'r': 500., 'l': 500.}, theta=.2)
            one = c.step(.01, x, descending={'vaslat_r': .5}); c.restore(saved)
            two = c.step(.01, x, descending={'vaslat_r': .5})
            check(f'idempotent: step(x) twice from one state ({decoder},{reflexes})', one == two)


def mode_identities():
    ankle = [f'{n}_{s}' for s in 'rl' for n in ('tibant', 'soleus', 'gasmed', 'gaslat')]
    a = stepping(SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG), 40)
    b = stepping(SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, reflexes='geyer_herr_2010'), 40)
    same = all(fa['motor_excitations'][m] == fb['motor_excitations'][m] for fa, fb in zip(a, b) for m in ankle)
    check('extended reflex set reproduces the eight ankle primitives exactly', same)
    moved = sum(fa['motor_excitations']['vaslat_r'] != fb['motor_excitations']['vaslat_r'] for fa, fb in zip(a, b))
    check('extended reflex set actually drives VAS (control can fail)', moved > 10, f'{moved} of 40 frames differ')
    flat = SensorimotorParameters(cortical_gain_per_hz=0.)
    for reflexes in (None, 'geyer_herr_2010'):
        r = SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, parameters=flat, reflexes=reflexes)
        q = SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, parameters=flat, reflexes=reflexes, decoder='recruitment')
        fr, fq = [], []
        for i in range(30):
            obs = observation(r.time_s, 1 + .04 * math.sin(i), 400., loads={'r': 500., 'l': 0.}, theta=.2)
            fr.append(r.step(.01, obs)); fq.append(q.step(.01, obs))
        diff = max(abs(x['motor_excitations'][m] - y['motor_excitations'][m]) for x, y in zip(fr, fq) for m in r.muscles)
        check(f'recruitment at zero descending == rate decoder with dial off ({reflexes})', diff == 0., f'max |diff| {diff:g}')
    # The difference the pool makes is ONLY on descending drive, and it is size-ordered.
    q = SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, parameters=flat, decoder='recruitment')
    r = SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, parameters=flat)
    for _ in range(5):
        obs = observation(q.time_s)
        fq = q.step(.01, obs, descending={'recfem_r': 1.}); fr = r.step(.01, obs, descending={'recfem_r': 1.})
    check('descending drive decodes differently under the pool', fq['requested_excitations']['recfem_r'] != fr['requested_excitations']['recfem_r'],
          f"rate {fr['requested_excitations']['recfem_r']:.4f} vs pool {fq['requested_excitations']['recfem_r']:.4f}")
    check('recruitment frame says it produces excitation only', 'excitation only' in fq['excitation_owner'] and fq['activation_owner'] == 'mechanical_plant')


def law_signs():
    base = dict(load_bw=1., contra_load_bw=0., dsup=False, knee_phi=2.5, knee_phi_rate=0., theta=T['theta_ref'], theta_rate=0., theta_takeoff=None)
    f = {'SOL': .3, 'GAS': .3, 'VAS': .5, 'HAM': .2, 'GLU': .2}
    l = {'TA': 1., 'HFL': .9, 'HAM': .9}
    S = lambda g, stance=True, **kw: geyer_herr_2010_stimulation(g, stance=stance, force=f, length=l, side_state=dict(base, **kw))[0]
    fwd, back = T['theta_ref'] + .1, T['theta_ref'] - .1
    check('forward lean in stance: HAM, GLU above S0,BAL; HFL at S0,BAL', S('HAM', theta=fwd) > T['S0_BAL'] and S('GLU', theta=fwd) > T['S0_BAL'] and S('HFL', theta=fwd) == T['S0_BAL'])
    check('backward lean in stance: HFL above S0,BAL; HAM, GLU at S0,BAL', S('HFL', theta=back) > T['S0_BAL'] and S('HAM', theta=back) == T['S0_BAL'] and S('GLU', theta=back) == T['S0_BAL'])
    check('trunk terms scale with leg load (unloaded stance leg: no trunk balance)', S('HAM', theta=fwd, load_bw=0.) == T['S0_BAL'])
    check('GLU uses 0.68 k_p', abs((S('GLU', theta=fwd) - T['S0_BAL']) - .68 * T['k_p'] * .1 * T['k_bw']) < 1e-12)
    check('VAS F+ in stance', abs(S('VAS') - (T['S0_VAS'] + T['G_VAS'] * .5)) < 1e-12)
    check('VAS knee overextension inhibition only when extending past 170 deg',
          S('VAS', knee_phi=3.05, knee_phi_rate=1.) < S('VAS') and S('VAS', knee_phi=3.05, knee_phi_rate=-1.) == S('VAS'))
    check('trailing leg: VAS inhibited by contra load, GLU -dS, HFL +dS',
          S('VAS', dsup=True, contra_load_bw=.5) < S('VAS') and S('GLU', dsup=True) == S('GLU') - T['delta_S'] and S('HFL', dsup=True) == S('HFL') + T['delta_S'])
    check('swing: SOL/GAS/VAS at S0 = 0.01', all(S(g, stance=False) == .01 for g in ('SOL', 'GAS', 'VAS')))
    check('swing: HAM and GLU force feedback', abs(S('HAM', stance=False) - (.01 + .65 * .2)) < 1e-12 and abs(S('GLU', stance=False) - (.01 + .4 * .2)) < 1e-12)
    check('swing HFL: HAM stretch past 0.85 suppresses HFL', S('HFL', stance=False) > geyer_herr_2010_stimulation('HFL', stance=False, force=f, length=dict(l, HAM=1.1), side_state=base)[0])
    check('swing HFL: take-off lean bias', abs(S('HFL', stance=False, theta_takeoff=.3) - S('HFL', stance=False) - T['k_lean'] * (.3 - T['theta_ref'])) < 1e-12)
    check("own input missing -> no law; cross-group missing -> term dropped",
          geyer_herr_2010_stimulation('VAS', stance=True, force={}, length=l, side_state=base) is None
          and geyer_herr_2010_stimulation('TA', stance=True, force={}, length=l, side_state=base)[0] == .01 + 1.1 * (1. - .71))


def trailing_memory():
    c = SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, reflexes='geyer_herr_2010')
    def go(loads, theta=.1):
        return c.step(.01, observation(c.time_s, loads=loads, theta=theta))
    for _ in range(4): f = go({'r': 500., 'l': 500.})
    check('double support from t0: order unknown -> neither leg trailing', f['reflex_state']['trailing'] == {'r': False, 'l': False})
    for _ in range(4): f = go({'r': 0., 'l': 900.}, theta=.25)
    check('take-off lean remembered at the stance->swing transition', abs(f['reflex_state']['memory']['r']['theta_takeoff'] - .25) < 1e-9,
          f['reflex_state']['memory']['r']['theta_takeoff'])
    for _ in range(4): f = go({'r': 500., 'l': 500.})
    check('right lands second -> left is the trailing leg', f['reflex_state']['trailing'] == {'r': False, 'l': True})
    try:
        SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, reflexes='geyer_herr_2010').step(.01, {k: v for k, v in observation(0.).items() if k != 'native_bodies'})
        check('extended reflexes refuse a frame without bodies', False)
    except ValueError:
        check('extended reflexes refuse a frame without bodies', True)
    for bad in ({'decoder': 'linear'}, {'reflexes': 'song'}, {'pool_parameters': object()}):
        try:
            SensorimotorController.from_root(ROOT, muscle_catalog=CATALOG, **bad)
            check(f'rejects {bad}', False)
        except (ValueError, TypeError):
            check(f'rejects {list(bad)[0]}', True)


def main():
    pool_known_answers()
    default_equals_head()
    idempotence_and_restore()
    mode_identities()
    law_signs()
    trailing_memory()
    print(f'\n{len(FAILED)} FAILED' if FAILED else '\nALL PASSED', flush=True)
    for name in FAILED: print('  FAILED:', name)
    sys.exit(1 if FAILED else 0)


if __name__ == '__main__':
    main()
