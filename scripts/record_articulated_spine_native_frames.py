#!/usr/bin/env python3
"""Native frames of `articulated_spine_v1`, for checking the anatomy poser's FK against Simbody.

    OPENBLAS_NUM_THREADS=1 nice -n 10 .venv/bin/python scripts/record_articulated_spine_native_frames.py

`scripts/verify_anatomy_pose.py` checks the pure-python forward kinematics in
`ihm/assembly/anatomy_pose.py` against `transform_ground` in frames the engine itself wrote.
Every stored frame on disk is of `engineering_stance_v1` (22 bodies); none exercises the
variant's 25 bodies, its two UniversalJoints (wrists) or its un-welded subtalars.  This records
some, and nothing else.

Two sessions, opened ONE AT A TIME (NativeMechanicalStream prlimits each engine to 4 GB), in the
`free` environment.  Each starts from an explicit initial pose in which EVERY one of the fifteen
new coordinates is off zero -- one session on the upper side of its declared range, one on the
lower -- together with a handful of existing coordinates, so that a wrong axis order, a wrong
sign or a wrong joint frame on any new joint shows up as a transform mismatch rather than
cancelling at zero.  The fractions (0.5 of each declared bound) are chosen here, not measured;
they only need to be inside the range and far from zero.  Six frames per session: the initial
state and five 10 ms steps with no excitation.

Writes data/derived/anatomy-segment-binding-articulated-spine-v1/native_frames.json (gitignored).
Only `bodies[*].transform_ground` and `coordinates` are kept, which is what the FK check reads.
"""
from __future__ import annotations

import json
import shutil
import sys
import time
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm.body_constants import MECHANICAL_TARGET_MASS_KG  # noqa: E402
from ihm.native.mechanical_stream import NativeMechanicalStream  # noqa: E402

VARIANT = 'data/models/articulated_spine_v1'
REGISTRATION = VARIANT + '/registration.json'
BASE_MODEL = ROOT / 'data/models/engineering_stance_v1/model.osim'
OUT = ROOT / 'data/derived/anatomy-segment-binding-articulated-spine-v1/native_frames.json'
FRACTION = 0.5
EXISTING = {'knee_angle_r': 0.6, 'hip_flexion_l': 0.3, 'elbow_flex_r': 0.8, 'lumbar_bending': 0.1,
            'ankle_angle_l': 0.2, 'pro_sup_l': 0.5, 'mtp_angle_r': -0.3, 'arm_add_l': -0.2}


def ranges(path):
    model = ET.parse(path).getroot().find('Model')
    out = {}
    for c in model.iter('Coordinate'):
        r = c.findtext('range')
        out[c.get('name')] = tuple(float(v) for v in r.split()) if r else None
    return out


def main() -> int:
    variant, base = ranges(ROOT / VARIANT / 'model.osim'), ranges(BASE_MODEL)
    new = sorted(set(variant) - set(base))
    if len(new) != 15:
        raise SystemExit(f'expected the 15 coordinates the variant adds, found {len(new)}: {new}')
    sessions = {}
    for side, pick in (('upper', 1), ('lower', 0)):
        pose = dict(EXISTING)
        pose.update({c: FRACTION * variant[c][pick] for c in new})
        output = 'data/derived/record-spine-native-' + uuid.uuid4().hex[:10]
        stream = None
        try:
            stream = NativeMechanicalStream(ROOT, ROOT / output, environment='free',
                                            target_mass_kg=MECHANICAL_TARGET_MASS_KG,
                                            augmented_registration=REGISTRATION, initial_pose=pose)
            frames = [stream.snapshot()]
            for _ in range(5):
                frames.append(stream.advance(0.01))
        finally:
            if stream is not None:
                stream.close()
            shutil.rmtree(ROOT / output, ignore_errors=True)
        if len(frames[0]['bodies']) != 25:
            raise SystemExit(f'{side}: engine reports {len(frames[0]["bodies"])} bodies, not 25')
        read_back = {c: frames[0]['coordinates'][c]['value'] for c in new}
        sessions[side] = dict(
            initial_pose_requested=pose, new_coordinates_read_back_at_t0=read_back,
            frames=[dict(kind=f.get('kind'), time_s=f['time_s'],
                         bodies={b: {'transform_ground': r['transform_ground']} for b, r in f['bodies'].items()},
                         coordinates={k: {'value': v['value']} for k, v in f['coordinates'].items()})
                    for f in frames])
        print(f'{side}: {len(frames)} frames, t = {frames[-1]["time_s"]:.3f} s; new coordinates read back '
              f'at t0: ' + ', '.join(f'{c}={read_back[c]:+.4f}' for c in new))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(dict(
        schema='ihm.articulated-spine-native-frames.v1', model=VARIANT + '/model.osim',
        registration=REGISTRATION, environment='free', recorded_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        purpose='FK check of ihm/assembly/anatomy_pose.py against Simbody on the 25-body variant',
        sessions=sessions)) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
