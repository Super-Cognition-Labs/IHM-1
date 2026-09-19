"""Does the upper body MOVE under muscle excitation alone, with every torque port at zero?

A moment arm says a path crosses a joint. It does not say the joint moves. This
drives the plant for a bounded interval with the stance bundle's own equilibrium
excitations held everywhere except one named agonist group, which is driven to
full excitation, and reports the coordinate excursion.

`coordinate_actuation` is NEVER passed, so all 13 declared `CoordinateActuator`
torque ports stay at zero for the whole run (`NativeMechanicalStream.advance`:
"they start at zero"). Anything that moves in the arms or trunk therefore moved
because a muscle pulled it.

Controls:

  PAIRED NULL   the same interval is run first with the equilibrium excitations
                unchanged. The driven excursion is reported as a DIFFERENCE from
                that null, because the bundle's equilibrium is not exactly static
                and a raw excursion would credit the drive with the drift.
  CONTRALATERAL the undriven side is reported beside the driven one. A drive that
                moves both sides equally is not a drive, it is the whole body
                settling.
  STEP BUDGET   a hard wall-clock and step cap. `advance()` cost is unbounded by
                dt (IHM-1 CLAUDE.md); a run that blows the budget is recorded as
                a budget failure, never silently truncated into a result.

    OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
      .venv/bin/python scripts/probe_upper_body_muscle_drive.py \
      --group elbow_r --output data/derived/upper-body-drive-elbow-r
"""
import argparse
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STANCE = 'data/models/engineering_stance_v1/registration.json'

GROUPS = {
    'elbow_r': (('arm26_BIClong_r', 'arm26_BICshort_r', 'arm26_BRA_r'), 'elbow_flex_r', 'elbow_flex_l'),
    'elbow_l': (('arm26_BIClong_l', 'arm26_BICshort_l', 'arm26_BRA_l'), 'elbow_flex_l', 'elbow_flex_r'),
    'triceps_r': (('arm26_TRIlong_r', 'arm26_TRIlat_r', 'arm26_TRImed_r'), 'elbow_flex_r', 'elbow_flex_l'),
    'shoulder_r': (('arm26_BICshort_r', 'arm26_TRIlong_r'), 'arm_add_r', 'arm_add_l'),
    'lumbar_ext': (('gait2392_ercspn_r', 'gait2392_ercspn_l'), 'lumbar_extension', 'pelvis_tilt'),
}


def run(stream, excitations, steps, dt, watched, budget_s):
    """Advance `steps` and return the watched coordinate trace, or a budget failure."""
    started = time.monotonic()
    trace = []
    for index in range(steps):
        stream.advance(dt, actuation=excitations)
        trace.append({name: stream.state['coordinates'][name]['value'] for name in watched})
        if time.monotonic() - started > budget_s:
            return trace, {'budget_exceeded': True, 'steps_completed': index + 1,
                           'wall_clock_s': time.monotonic() - started}
    return trace, {'budget_exceeded': False, 'steps_completed': steps,
                   'wall_clock_s': time.monotonic() - started}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--group', default='elbow_r', choices=sorted(GROUPS))
    parser.add_argument('--steps', type=int, default=40)
    parser.add_argument('--dt', type=float, default=0.01)
    parser.add_argument('--budget-s', type=float, default=240.)
    parser.add_argument('--output', required=True, help='root-relative fresh output directory')
    args = parser.parse_args()

    from ihm.native.mechanical_stream import NativeMechanicalStream

    manifest = json.loads((ROOT / STANCE).read_text())
    equilibrium = dict(manifest['equilibrium_excitations'])
    driven, watched_coordinate, mirror_coordinate = GROUPS[args.group]
    watched = (watched_coordinate, mirror_coordinate)
    output = ROOT / args.output
    report = {'schema': 'ihm.upper-body-muscle-drive.v1', 'group': args.group,
              'driven_muscles': list(driven), 'registration': STANCE,
              'dt_s': args.dt, 'steps': args.steps,
              'coordinate_actuation_commanded': False,
              'watched_coordinate': watched_coordinate, 'mirror_coordinate': mirror_coordinate}

    # Control arm first: the same interval with nothing changed.
    stream = NativeMechanicalStream(ROOT, output / 'null', environment='upright',
                                    target_mass_kg=float(manifest['target_mass_kg']),
                                    augmented_registration=STANCE,
                                    initial_pose=manifest['initial_pose'])
    try:
        start = {name: stream.state['coordinates'][name]['value'] for name in watched}
        null_trace, null_budget = run(stream, equilibrium, args.steps, args.dt, watched, args.budget_s)
    finally:
        stream.close()
    report['start_rad'] = start
    report['null'] = {'budget': null_budget,
                      'end_rad': null_trace[-1] if null_trace else None,
                      'excursion_rad': {n: null_trace[-1][n] - start[n] for n in watched} if null_trace else None}

    if null_budget['budget_exceeded']:
        report['verdict'] = 'BUDGET FAILURE in the null arm; no drive result is reported'
        output.mkdir(parents=True, exist_ok=True)
        (output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
        print(json.dumps(report['null'], indent=2))
        return

    # Drive arm: identical except the named agonists at full excitation.
    excited = dict(equilibrium)
    for name in driven:
        if name not in excited:
            raise ValueError('Driven muscle is not in the bundle equilibrium: ' + name)
        excited[name] = 1.0
    stream = NativeMechanicalStream(ROOT, output / 'driven', environment='upright',
                                    target_mass_kg=float(manifest['target_mass_kg']),
                                    augmented_registration=STANCE,
                                    initial_pose=manifest['initial_pose'])
    try:
        start_driven = {name: stream.state['coordinates'][name]['value'] for name in watched}
        if start_driven != start:
            raise ValueError('The two arms did not start at the same state')
        drive_trace, drive_budget = run(stream, excited, args.steps, args.dt, watched, args.budget_s)
    finally:
        stream.close()

    report['driven'] = {'budget': drive_budget,
                        'end_rad': drive_trace[-1] if drive_trace else None,
                        'excursion_rad': {n: drive_trace[-1][n] - start[n] for n in watched} if drive_trace else None}
    if drive_budget['budget_exceeded'] or len(drive_trace) != len(null_trace):
        report['verdict'] = 'BUDGET FAILURE in the driven arm; the paired difference is not comparable'
    else:
        difference = {n: drive_trace[-1][n] - null_trace[-1][n] for n in watched}
        report['paired_difference_rad'] = difference
        report['paired_difference_deg'] = {n: v * 57.29577951308232 for n, v in difference.items()}
        moved = difference[watched_coordinate]
        mirrored = difference[mirror_coordinate]
        report['driven_minus_mirror_rad'] = moved - mirrored
        report['verdict'] = ('muscle drive moved the watched coordinate'
                             if abs(moved) > 1e-3 and abs(moved) > 3 * abs(mirrored)
                             else 'no separable muscle-driven excursion at this budget')

    output.mkdir(parents=True, exist_ok=True)
    (output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: report[k] for k in report
                      if k not in ('null', 'driven')}, indent=2))
    print('null   ', report['null'])
    print('driven ', report['driven'])
    print('report:', (output / 'report.json').relative_to(ROOT))


if __name__ == '__main__':
    main()
