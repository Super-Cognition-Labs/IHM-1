"""What the reassigned muscles can actually do about the thoracic joint.

    OPENBLAS_NUM_THREADS=1 prlimit --as=4294967296 -- nice -n 10 \
        .venv/bin/python -m scripts.measure_thoracic_drive_capacity

Reads every muscle's moment arm about the three thoracic coordinates out of the
engine, on `thoracic_drive_v1` and on its base `shoulder_girdle_v1`, and sums
`Fmax x |r|` in each direction.

READ THE ZEROS AND THE SIGNS, NOT THE MAGNITUDES.  `Sum Fmax x |r|` at one pose
is full activation at optimal fibre length with no force-length, no
force-velocity, no tendon state and no activation dynamics: a CEILING. A
per-axis sum is not a tension-feasible torque cone. docs/UPPER_BODY_ACTUATION.md
section 2 attaches that caution to this exact quantity and it applies unchanged
here. Nothing in this script integrates, so nothing here is behaviour.

TWO POPULATIONS, REPORTED SEPARATELY AND NEVER SUMMED INTO ONE COUNT.
  * ATTACHMENT-DRIVEN: the ten muscles whose torso attachment moved onto the
    thorax. Their thoracic moment arm comes from an endpoint that the
    pre-registered measurement placed on a rib, the sternum or an intercostal.
  * WRAP-ONLY: muscles that keep both attachments where they were but wrap on
    the rib-cage ellipsoid, which moved with the rib cage. Their thoracic
    moment arm comes from a wrap surface alone. That is a real mechanism and a
    much weaker claim, and it is reported as its own row.

OUTPUT: data/models/thoracic_drive_v1/thoracic_capacity.json
"""
from __future__ import annotations

import json
import pathlib
import shutil
import sys
import uuid
import xml.etree.ElementTree as ET

from ihm.native.mechanical_stream import NativeMechanicalStream

ROOT = pathlib.Path(__file__).resolve().parents[1]
VARIANT_DIR = ROOT / 'data/models/thoracic_drive_v1'
VARIANT = 'data/models/thoracic_drive_v1/registration.json'
BASE = 'data/models/shoulder_girdle_v1/registration.json'
TARGET_MASS_KG = 77.6122029
THORACIC = ('thoracic_extension', 'thoracic_bending', 'thoracic_rotation')
#: Below this a "moment arm" is not a lever, it is the wrap solver's residual.
#: docs/UPPER_BODY_ACTUATION.md section 2 uses the same 1 mm floor.
ARM_FLOOR_M = 1e-3


def classify():
    root = ET.parse(VARIANT_DIR / 'model.osim').getroot()
    attachment, wrap_only = [], []
    for muscle in root.iter():
        if not (muscle.tag.endswith('Muscle') and muscle.get('name')):
            continue
        bodies = {p.findtext('socket_parent_frame') for p in muscle.iter('PathPoint')}
        wraps = {w.findtext('wrap_object') for w in muscle.iter('PathWrap')}
        if '/bodyset/thorax' in bodies:
            attachment.append(muscle.get('name'))
        elif wraps & {'Thorax_r', 'Thorax_l'}:
            wrap_only.append(muscle.get('name'))
    return sorted(attachment), sorted(wrap_only)


def read(registration, muscles):
    out = 'data/derived/thoracic-capacity-' + uuid.uuid4().hex[:10]
    stream = NativeMechanicalStream(ROOT, ROOT / out, environment='free',
                                    target_mass_kg=TARGET_MASS_KG,
                                    augmented_registration=registration)
    try:
        state = stream.snapshot()
        first = stream.moment_arms(muscles=muscles, coordinates=list(THORACIC))
        second = stream.moment_arms(muscles=muscles, coordinates=list(THORACIC))
        if first['moment_arms_m'] != second['moment_arms_m']:
            raise SystemExit('moment_arms is not a function of its input')
        fmax = {m: state['muscles'][m]['max_isometric_force_n'] for m in muscles}
        return first['moment_arms_m'], fmax
    finally:
        stream.close()
        shutil.rmtree(ROOT / out, ignore_errors=True)


def bounds(arms, fmax, names):
    out = {}
    for coordinate in THORACIC:
        plus = sum(fmax[m] * arms[m][coordinate]
                   for m in names if arms[m][coordinate] > ARM_FLOOR_M)
        minus = sum(fmax[m] * -arms[m][coordinate]
                    for m in names if arms[m][coordinate] < -ARM_FLOOR_M)
        crossing = [m for m in names if abs(arms[m][coordinate]) > ARM_FLOOR_M]
        out[coordinate] = {'muscles_crossing': len(crossing),
                           'plus_bound_nm': plus, 'minus_bound_nm': minus,
                           'crossing': sorted(crossing)}
    return out


def main() -> int:
    attachment, wrap_only = classify()
    everything = sorted(set(attachment) | set(wrap_only))

    arms, fmax = read(VARIANT, everything)
    base_arms, _ = read(BASE, everything)

    base_worst = max((abs(base_arms[m][c]), m, c)
                     for m in everything for c in THORACIC)
    if base_worst[0] >= 1e-9:
        raise SystemExit('the base plant already had a thoracic arm of %.3e m at '
                         '%s/%s; the control that makes this measurement mean '
                         'anything has failed' % base_worst)

    report = {
        'schema': 'ihm.thoracic-capacity.v1',
        'model': 'data/models/thoracic_drive_v1/model.osim',
        'pose': "each model's own rest pose; thoracic_* = 0",
        'quantity': 'Sum Fmax x |moment arm| per direction, a CEILING at one '
                    'pose: full activation at optimal fibre length, no '
                    'force-length, no force-velocity, no tendon state, no '
                    'activation dynamics, and not a tension-feasible torque cone',
        'arm_floor_m': ARM_FLOOR_M,
        'control_base_plant_worst_thoracic_arm_m': base_worst[0],
        'control_basis': 'every muscle probed must read ZERO about every thoracic '
                         'coordinate on the base plant, which is what '
                         'docs/ARTICULATED_SPINE.md measured; a nonzero reading '
                         'there would void this whole measurement',
        'attachment_driven': {
            'muscles': attachment,
            'basis': 'torso attachment moved onto the thorax under the rule in '
                     'scripts/measure_thoracic_reassignability.py',
            'bounds': bounds(arms, fmax, attachment)},
        'wrap_only': {
            'muscles': wrap_only,
            'basis': 'both attachments unchanged; the thoracic dependence comes '
                     'from the rib-cage wrap ellipsoid alone, which moved with '
                     'the rib cage. A WEAKER mechanism than an attachment and '
                     'deliberately not added to the count above.',
            'bounds': bounds(arms, fmax, wrap_only)},
        'moment_arms_m': {m: arms[m] for m in everything},
        'max_isometric_force_n': fmax,
    }
    (VARIANT_DIR / 'thoracic_capacity.json').write_text(
        json.dumps(report, indent=2) + '\n')

    print('control: worst thoracic moment arm on the BASE plant = %.2e m '
          '(must be zero)\n' % base_worst[0])
    for label, group in (('ATTACHMENT-DRIVEN', 'attachment_driven'),
                         ('WRAP-ONLY', 'wrap_only')):
        print('%s  (%d muscles)' % (label, len(report[group]['muscles'])))
        print('  %-20s %9s %12s %12s' % ('coordinate', 'crossing', '+bound N.m',
                                         '-bound N.m'))
        for coordinate in THORACIC:
            row = report[group]['bounds'][coordinate]
            print('  %-20s %9d %12.1f %12.1f'
                  % (coordinate, row['muscles_crossing'],
                     row['plus_bound_nm'], row['minus_bound_nm']))
        print()
    print('per-muscle thoracic arms, attachment-driven (mm):')
    for muscle in attachment:
        print('  %-34s %8.1f %8.1f %8.1f'
              % (muscle, 1000 * arms[muscle]['thoracic_extension'],
                 1000 * arms[muscle]['thoracic_bending'],
                 1000 * arms[muscle]['thoracic_rotation']))
    print('\nreport: %s'
          % (VARIANT_DIR / 'thoracic_capacity.json').relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
