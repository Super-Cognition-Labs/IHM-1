#!/usr/bin/env python3
"""How far apart are the two poses of the same 4,000 entities?

The live plant carries `CanonicalRegistration.project()`: ONE shared global rigid
fit of 22 approximate COM/bone-envelope-centre correspondences, with each entity
assigned to a segment by bounding-box distance. `ihm/assembly/anatomy_pose.py`
carries `AnatomyPoser`: the per-entity segment binding, the declared 0.963
similarity, and forward kinematics verified against Simbody to 7.8e-16
(`scripts/verify_anatomy_pose.py`).

Both pose the SAME 4,000 ids in the SAME frame (`bodyparts3d-display-m`), so the
disagreement is a number rather than an argument. This measures it over real
native frames, because at the reference pose the two agree by construction and a
registration error is invisible until something rotates -- which is the whole
failure mode this repo keeps recording.

    PYTHONPATH=. .venv/bin/python scripts/measure_display_pose_disagreement.py
"""
import argparse, json, time, uuid
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REGISTRATION = 'data/derived/mechanics/whole_body_arm26_v2/registration.json'


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--steps', type=int, default=25)
    ap.add_argument('--dt', type=float, default=0.02)
    ap.add_argument('--excitation', type=float, default=0.12)
    ap.add_argument('--out', default='data/derived/display-pose-disagreement/report.json')
    a = ap.parse_args()

    from ihm.assembly.articulated import ArticulatedBodyPlant
    from ihm.assembly.anatomy_pose import AnatomyPoser

    poser = AnatomyPoser.from_workspace(ROOT)
    out_dir = ROOT / f'data/derived/scratch-posecmp-{uuid.uuid4().hex[:8]}'
    plant = ArticulatedBodyPlant(ROOT, out_dir, environment='upright',
                                 target_mass_kg=77.6122029,
                                 augmented_registration=REGISTRATION)
    rows, costs = [], []
    try:
        native = plant.native.snapshot()
        commands = {k: a.excitation for k in native['muscles']}
        for step in range(a.steps + 1):
            if step:
                native = plant.native.advance(a.dt, actuation=commands)
            legacy = plant.registration.project(native)
            t0 = time.perf_counter()
            pose = poser.pose_from_native(native)
            costs.append(time.perf_counter() - t0)
            ids = [i for i in pose.entity_ids if i in legacy]
            index = {i: n for n, i in enumerate(pose.entity_ids)}
            # the SAME rest point for both, so the only thing that differs is the transform
            rest = np.array([plant.registration.specs[i]['centroid_m'] for i in ids])
            verified = np.array([pose.rotation[index[i]] @ rest[n] + pose.translation[index[i]]
                                 for n, i in enumerate(ids)])
            shipped = np.array([legacy[i]['centroid_m'] for i in ids])
            d = np.linalg.norm(verified - shipped, axis=1)
            rows.append({'step': step, 'time_s': native['time_s'], 'entities': len(ids),
                         'median_mm': float(np.median(d) * 1e3),
                         'mean_mm': float(d.mean() * 1e3),
                         'p95_mm': float(np.percentile(d, 95) * 1e3),
                         'max_mm': float(d.max() * 1e3),
                         'worst_entity': ids[int(d.argmax())]})
    finally:
        plant.native.close()

    last = rows[-1]
    report = {'schema': 'ihm.display-pose-disagreement.v1',
              'compared': 'CanonicalRegistration.project (shipped, one global rigid fit, '
                          'bounding-box segment assignment) against AnatomyPoser.pose_from_native '
                          '(per-entity binding.json, FK verified to 7.8e-16)',
              'frame': 'bodyparts3d-display-m', 'entities': last['entities'],
              'steps': a.steps, 'dt_s': a.dt, 'excitation': a.excitation,
              'poser_cost_ms_median': float(np.median(costs) * 1e3),
              'poser_cost_ms_max': float(max(costs) * 1e3),
              'per_step': rows,
              'reading': 'At the reference pose the two agree by construction; the disagreement '
                         'is what the shared rigid fit costs once segments rotate. Neither is a '
                         'measurement of the real body -- but only one of them has a known answer '
                         'against the integrator that produced the motion.'}
    target = ROOT / a.out
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2, default=lambda o: getattr(o, 'tolist', lambda: str(o))()) + '\n')
    print(f'entities compared: {last["entities"]}')
    print(f'poser cost: {report["poser_cost_ms_median"]:.2f} ms median, '
          f'{report["poser_cost_ms_max"]:.2f} ms max')
    print(f'{"step":>5} {"t (s)":>7} {"median":>9} {"mean":>9} {"p95":>9} {"max":>10}  worst')
    for r in rows[:: max(1, len(rows) // 8)] + [last]:
        print(f'{r["step"]:5d} {r["time_s"]:7.3f} {r["median_mm"]:8.2f}m {r["mean_mm"]:8.2f}m '
              f'{r["p95_mm"]:8.2f}m {r["max_mm"]:9.2f}m  {r["worst_entity"]}')
    print(f'written: {a.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
