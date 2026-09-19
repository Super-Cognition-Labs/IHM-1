"""Which coordinates of the live plant have MUSCLE authority, and which have none.

The register (`docs/WORKBENCH_AUTHENTICITY.md` 0.1) says the upper body is driven
by 13 `CoordinateActuator` torque motors. That is a statement about the SOURCE
model `subject_walk_scaled.osim`. It is not a statement about the plant that runs:
`EmbodiedRuntime.from_workspace` defaults `augmented_registration` to
`whole_body_arm26_v2` (92 muscles) and the stance controllers force
`data/models/engineering_stance_v1` (98), and `NativeMechanicalStream.advance`
leaves every torque port at zero unless a caller commands it, which nothing on the
live path does.

So the question that actually decides 0.1 is not "are there torque motors" but
"which rotational coordinates does a MUSCLE cross, and with what moment arm".
This asks the native engine directly, through `moment_arms`, which reads the
actual source paths (the 80 fitted `FunctionBasedPath` polynomials for the legs
and the registered `GeometryPath` + wrap objects for the additions).

Reported per coordinate:

  n_muscles          muscles with |moment arm| >= MIN_ARM_M
  agonist/antagonist peak isometric torque BOUND, sum(Fmax * |r|) by sign
  present_torque_n_m sum(tendon_force * r) in the plant's CURRENT state
  actuator_optimal   the declared CoordinateActuator optimal force, for contrast

The peak torque is an UPPER BOUND and nothing else. It multiplies max isometric
force by a moment arm at ONE pose: it ignores force-length, force-velocity, tendon
state, activation dynamics and the fact that a tension-feasible torque cone is a
3-D object, not a per-axis sum. `docs/research/LUMBAR_SHOULDER_MUSCLE_COVERAGE.md`
warns about exactly this multiplication; it is printed because a coordinate with a
bound of ZERO is decisive (no muscle crosses it at all) while a large bound is
merely an upper limit. Read the zeros, not the magnitudes.

Controls, both of which must pass or the numbers mean nothing:

  IDEMPOTENCE  `moment_arms` is called twice at the identical state and the two
               dicts must be bit-equal. A native query that is not a function of
               its arguments produces confident wrong conclusions (IHM-1
               CLAUDE.md, "Call it twice at the same input").
  KNOWN ANSWER the ankle must come back muscle-dominated with soleus plantarflexing
               at roughly 4-5 cm. If the legs do not read as legs, the instrument
               is wrong and no upper-body zero can be believed.

Run ONE native session at a time.

    OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
      .venv/bin/python scripts/measure_upper_body_actuation.py \
      --registration data/derived/mechanics/whole_body_arm26_v2/registration.json \
      --output data/derived/upper-body-actuation-arm26v2
"""
import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_MODEL = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'
MIN_ARM_M = 1e-3  # 1 mm: below this a crossing carries no usable torque
UPPER_BODY_BODIES = ('torso', 'humerus', 'ulna', 'radius', 'hand')


def torque_ports(model_path):
    """Declared CoordinateActuator ports of the model actually being loaded."""
    ports = {}
    for element in ET.parse(model_path).getroot().find('.//ForceSet/objects'):
        if element.tag != 'CoordinateActuator':
            continue
        coordinate = element.findtext('coordinate')
        if coordinate is None:
            raise ValueError('CoordinateActuator without a coordinate')
        ports[coordinate] = {'actuator': element.get('name'),
                             'optimal_force': float(element.findtext('optimal_force'))}
    return ports


def muscle_bodies(model_path):
    """Attachment bodies per muscle, read off the model that is being loaded."""
    bodies = {}
    for element in ET.parse(model_path).getroot().find('.//ForceSet/objects'):
        if 'Muscle' not in element.tag:
            continue
        found = sorted({p.text.rsplit('/', 1)[-1]
                        for p in element.findall('.//PathPointSet/objects/*/socket_parent_frame')
                        if p.text})
        bodies[element.get('name')] = found
    return bodies


def summarize(arms, muscles, ports, bodies):
    rows = {}
    for coordinate in sorted({c for row in arms.values() for c in row}):
        crossing = []
        agonist = antagonist = present = 0.
        for name, row in arms.items():
            arm = row[coordinate]
            if abs(arm) < MIN_ARM_M:
                continue
            state = muscles[name]
            fmax = float(state['max_isometric_force_n'])
            bound = fmax * abs(arm)
            if arm > 0:
                agonist += bound
            else:
                antagonist += bound
            present += float(state['tendon_force_n']) * arm
            crossing.append({'muscle': name, 'moment_arm_m': arm,
                             'max_isometric_force_n': fmax,
                             'peak_torque_bound_n_m': bound,
                             'present_torque_n_m': float(state['tendon_force_n']) * arm,
                             'attachment_bodies': bodies.get(name, [])})
        crossing.sort(key=lambda r: -abs(r['moment_arm_m']))
        port = ports.get(coordinate)
        rows[coordinate] = {
            'n_muscles_crossing': len(crossing),
            'positive_peak_torque_bound_n_m': agonist,
            'negative_peak_torque_bound_n_m': antagonist,
            'present_net_muscle_torque_n_m': present,
            'declared_torque_port': None if port is None else port['actuator'],
            'declared_torque_port_optimal_force': None if port is None else port['optimal_force'],
            'muscle_authority': 'none' if not crossing else ('unidirectional' if not (agonist and antagonist) else 'bidirectional'),
            'crossing': crossing,
        }
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--registration', default='data/derived/mechanics/whole_body_arm26_v2/registration.json',
                        help='root-relative augmented registration, or "none" for the bare source plant')
    parser.add_argument('--environment', default='upright', choices=('free', 'supine', 'upright'))
    parser.add_argument('--mass-kg', type=float, default=77.6122029)
    parser.add_argument('--output', required=True, help='root-relative fresh output directory')
    args = parser.parse_args()

    from ihm.native.mechanical_stream import NativeMechanicalStream

    registration = None if args.registration == 'none' else args.registration
    pose = None
    if registration is not None:
        manifest = json.loads((ROOT / registration).read_text())
        model_path = ROOT / manifest['model_path']
        pose = manifest.get('initial_pose')
        mass = float(manifest.get('target_mass_kg', args.mass_kg))
    else:
        model_path = ROOT / SOURCE_MODEL
        mass = args.mass_kg

    output = ROOT / args.output
    stream = NativeMechanicalStream(ROOT, output / 'plant', environment=args.environment,
                                    target_mass_kg=mass, augmented_registration=registration,
                                    **({'initial_pose': pose} if pose else {}))
    try:
        snapshot = stream.snapshot()
        muscles = snapshot['muscles']
        coordinates = {name: meta for name, meta in stream.state['coordinates'].items()
                       if meta['unit'] == 'rad'}
        names = sorted(muscles)
        first = stream.moment_arms(muscles=names, coordinates=sorted(coordinates))
        second = stream.moment_arms(muscles=names, coordinates=sorted(coordinates))
        idempotent = first['moment_arms_m'] == second['moment_arms_m']
        arms = first['moment_arms_m']

        ports = torque_ports(model_path)
        bodies = muscle_bodies(model_path)
        rows = summarize(arms, muscles, ports, bodies)

        # Known answer: the ankle must read as a muscle-driven ankle.
        control = {}
        for coordinate in ('ankle_angle_r', 'ankle_angle_l'):
            if coordinate in rows:
                control[coordinate] = rows[coordinate]['n_muscles_crossing']
        soleus = arms.get('soleus_r', {}).get('ankle_angle_r')
        control['soleus_r_ankle_moment_arm_m'] = soleus
        control['soleus_r_arm_in_4_to_6_cm'] = (
            soleus is not None and 0.04 <= abs(soleus) <= 0.06)
        control['ankles_muscle_driven'] = all(
            rows.get(c, {}).get('n_muscles_crossing', 0) >= 6 for c in ('ankle_angle_r', 'ankle_angle_l'))

        upper = {c: r for c, r in rows.items()
                 if any(b.startswith(UPPER_BODY_BODIES)
                        for m in r['crossing'] for b in m['attachment_bodies'])
                 or c in ports}
        unmuscled_ports = sorted(c for c in ports if rows.get(c, {}).get('n_muscles_crossing', 0) == 0)

        report = {
            'schema': 'ihm.upper-body-actuation.v1',
            'registration': registration,
            'model_path': str(model_path.relative_to(ROOT)),
            'environment': args.environment,
            'initial_pose_applied': bool(pose),
            'n_muscles': len(names),
            'n_rotational_coordinates': len(coordinates),
            'controls': {
                'moment_arms_idempotent': idempotent,
                'known_answer_ankle': control,
            },
            'declared_torque_ports': ports,
            'torque_ports_with_no_muscle_crossing': unmuscled_ports,
            'per_coordinate': rows,
            'upper_body_coordinates': sorted(upper),
            'bound_caveat': 'peak_torque_bound_n_m is Fmax*|r| at one pose: an upper bound only. '
                            'It ignores force-length, force-velocity, tendon and activation state, '
                            'and a per-axis sum is not a tension-feasible 3-D torque cone. '
                            'A ZERO is decisive; a magnitude is not.',
        }
        if not idempotent:
            report['controls']['verdict'] = 'VOID: native moment_arms returned different answers at the identical state'
        (output / 'report.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')

        print('muscles', len(names), 'registration', registration)
        print('idempotent', idempotent, '| ankle control', control['ankles_muscle_driven'],
              '| soleus_r arm', None if soleus is None else round(soleus, 4),
              'in 4-6cm', control['soleus_r_arm_in_4_to_6_cm'])
        print()
        header = f"{'coordinate':22s} {'nmus':>4s} {'+bound':>9s} {'-bound':>9s} {'now':>8s} {'authority':>14s} {'port':>16s}"
        print(header)
        for coordinate in sorted(ports):
            row = rows.get(coordinate)
            if row is None:
                continue
            print(f"{coordinate:22s} {row['n_muscles_crossing']:4d} "
                  f"{row['positive_peak_torque_bound_n_m']:9.1f} {row['negative_peak_torque_bound_n_m']:9.1f} "
                  f"{row['present_net_muscle_torque_n_m']:8.2f} {row['muscle_authority']:>14s} "
                  f"{str(row['declared_torque_port']):>16s}")
        print()
        print('torque ports with NO muscle crossing:', unmuscled_ports or 'none')
        print('report:', (output / 'report.json').relative_to(ROOT))
    finally:
        stream.close()


if __name__ == '__main__':
    main()
