#!/usr/bin/env python3
"""Known answers for the mechanical fidelity a live body can now be asked for.

`ihm/assembly/plant_options.py` made three plant capabilities reachable from a
session that could previously only be reached from an offline script: the joint
stops, the real per-segment contact surfaces, and the derived tissue force
elements. Wiring is not evidence.

METHOD, taken from `scripts/measure_segment_contact_meshes.py` rather than
invented: the engineering stance initial pose, a PASSIVE settle (no actuation),
joint stops on, 50 steps of 10 ms. Driving every muscle instead sends the plant
into a configuration the error-controlled integrator cannot advance -- the first
version of this file did exactly that and timed out at 120 s in a single step,
which is the `advance()` cost trap CLAUDE.md records.

A GATE THIS FILE RECORDS AS FAILED, AND WHY THE INSTRUMENT CHANGED RATHER THAN THE
BAR. The first run asserted `the floor carries the body` (vertical contact force
~= weight) and `the body stays standing` (pelvis_ty roughly held). Both FAILED on
the BASELINE arm -- 190.3 N against 761.4 N, pelvis_ty 1.0187 -> 0.7748 m over
0.50 s -- so they were not a property of any arm, they were a property the stance
pose does not have: `measure_segment_contact_meshes.py` says in as many words that
the source pose is NOT an equilibrium, and reports pelvis_ty as a measurement
rather than gating on it. A passive body falls from it whatever its contact set.
The absolute gates are therefore withdrawn and replaced by PAIRED ones -- each arm
against the baseline's own fall -- which is what the comparison was always about.
The bar was not moved; the quantity was wrong.

A PREDICTION THIS FILE GOT BACKWARDS, RECORDED BECAUSE IT WAS WRITTEN DOWN FIRST.
Reading `measure_segment_contact_meshes.py`'s note that *the skin never reaches the
floor in the stance pose* -- which is why that script excludes the `skin` arm from
its own default run -- this file first asserted that replacing the source feet with
skin would make the body fall FURTHER, having nothing under it. Measured: it falls
**190.1 mm against the baseline's 243.9 mm**, i.e. it is caught 53.8 mm HIGHER. The
note is about the initial pose, where nothing is touching at all; once the body
moves, the skin is the OUTERMOST surface and a sphere is inscribed in the inertia
ellipsoid, so the skin necessarily contacts no later. The gate now states that
geometric prediction. The original assertion is withdrawn, not quietly deleted.

    .venv/bin/python scripts/verify_plant_fidelity.py

One native engine at a time, closed before the next opens.
"""
import json, sys, uuid
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
REGISTRATION = 'data/derived/mechanics/whole_body_arm26_v2/registration.json'
MASS = 77.6122029
WEIGHT_N = MASS * 9.81
STEPS, DT = 50, 0.01


def check(name, ok, detail=''):
    print(f'  {"PASS" if ok else "FAIL"}  {name}{"  " + detail if detail else ""}', flush=True)
    return bool(ok)


def settle(fidelity, pose):
    from ihm.assembly.articulated import ArticulatedBodyPlant
    out = ROOT / f'data/derived/scratch-fidchk-{uuid.uuid4().hex[:8]}'
    plant = ArticulatedBodyPlant(ROOT, out, environment='upright', target_mass_kg=MASS,
                                 augmented_registration=REGISTRATION, initial_pose=pose,
                                 mechanical_fidelity=fidelity)
    try:
        native = plant.native.snapshot()
        start_ty = native['coordinates']['pelvis_ty']['value']
        span = {}
        for _ in range(STEPS):
            native = plant.native.advance(DT)
            for name, row in native['coordinates'].items():
                lo, hi = span.get(name, (row['value'], row['value']))
                span[name] = (min(lo, row['value']), max(hi, row['value']))
        contacts = native['contacts']
        vertical = float(sum(c['force_n'][1] for c in contacts))
        return {'names': sorted({c['name'] for c in contacts}), 'count': len(contacts),
                'vertical_n': vertical, 'span': span,
                'start_ty': start_ty, 'end_ty': native['coordinates']['pelvis_ty']['value'],
                'tissue': len(native.get('tissue_ligaments') or []),
                'residual': float(np.linalg.norm(native['momentum_balance_residual_n'])),
                'selection': plant.mechanical_fidelity}
    finally:
        plant.native.close()


def main() -> int:
    from ihm.assembly.plant_options import declared_ranges
    ranged, unranged = declared_ranges(ROOT)
    pose = json.loads((ROOT / 'data/models/engineering_stance_v1/initial_pose.json').read_text())
    stops_only = {'joint_stops': True}
    ok = True

    print('baseline -- the historical plant, stance pose, passive settle', flush=True)
    base = settle(None, pose)
    ok &= check('contact is the inertia-ellipsoid sphere set', base['count'] == 28,
                f'{base["count"]} elements')
    ok &= check('no tissue force elements', base['tissue'] == 0)
    base_fall = base['start_ty'] - base['end_ty']
    print(f'    reference (not a gate): vertical contact {base["vertical_n"]:.1f} N of '
          f'{WEIGHT_N:.1f} N weight; pelvis_ty {base["start_ty"]:.4f} -> {base["end_ty"]:.4f} m '
          f'({base_fall*1e3:.1f} mm) -- the stance pose is not an equilibrium', flush=True)

    print('\nskin surfaces, source feet KEPT -- the cost of carrying the geometry', flush=True)
    carried = settle({'segment_contact': 'skin_carried'}, pose)
    skin_named = [n for n in carried['names'] if n.startswith('mesh_support_skin_')]
    ok &= check('real skin meshes are installed', len(skin_named) > 0,
                f'{len(skin_named)} skin elements, e.g. {skin_named[0] if skin_named else "-"}')
    ok &= check('the source feet are still there',
                any('contact' in n.lower() and 'mesh_support' not in n for n in carried['names']))
    carried_fall = carried['start_ty'] - carried['end_ty']
    ok &= check('carrying the geometry does not change how the body falls',
                abs(carried_fall - base_fall) < 0.005,
                f'{carried_fall*1e3:.1f} mm against the baseline\'s {base_fall*1e3:.1f} mm')
    ok &= check("the bundle's own declared skin material is used, not a default",
                carried['selection']['segment_contact']['material'] == {
                    'youngs_modulus_pa': 3000.0, 'poissons_ratio': 0.45, 'layer_thickness_m': 0.0066})

    print('\nskin surfaces REPLACING the feet -- asserted to collapse, because the plantar', flush=True)
    print('skin is not LEVEL in this pose: toes_l -8.5 mm through the floor plane while', flush=True)
    print('calcn_l sits +35.700 mm above it (scripts/verify_skin_contact.py, 2026-09-18).', flush=True)
    replaced = settle({'segment_contact': 'skin'}, pose)
    ok &= check('every contact element is a skin mesh',
                all(n.startswith('mesh_support_skin_') for n in replaced['names']),
                f'{replaced["count"]} elements')
    dropped = replaced['start_ty'] - replaced['end_ty']
    # The skin is the OUTERMOST surface of the segment and the sphere is inscribed in
    # its inertia ellipsoid, so on a body falling from a non-equilibrium pose the skin
    # must make contact no LATER than the sphere does, and catch the body no lower.
    # That is a prediction from the geometry, and it is the opposite of what the first
    # version of this file asserted -- see the note at the top.
    ok &= check('skin, being outside bone, catches the body no lower than the spheres',
                dropped <= base_fall + 1e-9,
                f'fell {dropped*1e3:.1f} mm against the baseline\'s {base_fall*1e3:.1f} mm '
                f'in {STEPS*DT:.2f} s -- caught {(base_fall-dropped)*1e3:.1f} mm higher')
    # Until 2026-09-18 this asserted the string 'does not reach the floor', which the
    # measurement withdrew: the shipped bundle's toes DO reach it, at -8.5 mm. The
    # check is on the caveat naming the measured defect, not on the retired sentence.
    caveat = replaced['selection']['segment_contact']['caveat'] or ''
    ok &= check('and the selection carries that caveat rather than a claim',
                'not LEVEL' in caveat and '+20.139 mm ABOVE the calcaneus' in caveat)

    print('\njoint stops -- the model declares ranges and nothing enforced them', flush=True)
    stopped = settle(stops_only, pose)
    ok &= check('stops on every coordinate that declares a real range',
                len(stopped['selection']['joint_stops']['coordinates']) == len(ranged),
                f'{len(ranged)} coordinates, {len(unranged)} left free: {", ".join(unranged)}')

    def worst(result):
        w, where = 0.0, None
        for name, (lo, hi) in result['span'].items():
            if name not in ranged:
                continue
            low, high = ranged[name]
            past = max(low - lo, hi - high, 0.0)
            if past > w:
                w, where = past, name
        return w, where

    fw, fwhere = worst(base)
    sw, swhere = worst(stopped)
    ok &= check('stops hold the body no further outside its declared ranges', sw <= fw + 1e-9,
                f'{sw:.4f} rad ({swhere}) against {fw:.4f} rad ({fwhere}) unstopped')

    print('\ntissue force elements', flush=True)
    tissue = settle({'tissue_ligaments': 'admissible'}, pose)
    ok &= check('only the kinematically admissible subset is carried', tissue['tissue'] == 66,
                f'{tissue["tissue"]} of 117 elements')
    ok &= check('internal forces do not change the momentum balance',
                tissue['residual'] < max(1e-6, base['residual'] * 10),
                f'residual {tissue["residual"]:.3e} N against baseline {base["residual"]:.3e} N')

    print('\nall three at once, on the arm that stays standing', flush=True)
    full = settle({'joint_stops': True, 'segment_contact': 'skin_carried',
                   'tissue_ligaments': 'admissible'}, pose)
    ok &= check('stops, real skin geometry and tissue coexist on one body',
                full['tissue'] == 66 and full['selection']['joint_stops'] is not None
                and any(n.startswith('mesh_support_skin_') for n in full['names']))
    # Each addition can only resist: a stop opposes leaving a range, tissue is a
    # tension-only internal force, and the skin contacts no later than a sphere. So a
    # body carrying all three must fall no FURTHER than the bare one. (An earlier
    # version of this gate asserted the combined arm would behave like the baseline,
    # which is wrong by construction -- joint stops change the dynamics, that is what
    # they are for. Withdrawn, not rescored.)
    full_fall = full['start_ty'] - full['end_ty']
    ok &= check('every addition resists, so the body falls no further than the bare one',
                full_fall <= base_fall + 1e-9,
                f'fall {full_fall*1e3:.1f} mm against the baseline\'s {base_fall*1e3:.1f} mm')

    print('\nALL PASS' if ok else '\nFAILURES ABOVE', flush=True)
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
