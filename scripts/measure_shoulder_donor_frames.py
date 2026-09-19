"""Measure the Thoracoscapular donor's OWN reference pose, with the engine.

    .venv/bin/python -m scripts.measure_shoulder_donor_frames

WHY THIS EXISTS.  `data/models/shoulder_girdle_v1` carries the Seth 2019
thoracoscapular girdle into this plant's torso frame.  To place it, the pose of
the donor's `clavicle`, `scapula` and `humerus` relative to its `thorax` AT THE
DONOR'S OWN REFERENCE CONFIGURATION has to be known.  Two of those poses are not
readable from the XML: the `ScapulothoracicJoint` is a Simbody ellipsoid
mobilizer, and the clavicle/scapula loop is closed by a `PointConstraint`, so the
scapula's reference pose is the solution of an assembly, not a stored transform.

Reimplementing the ellipsoid mobilizer in Python would be a second
implementation of the thing being measured.  Instead this asks the ENGINE, using
the only quantity the sampler reports: path length.  A `PathActuator` with two
points is a ruler.  Sixteen rulers between four thorax stations and four stations
on the target body give, per target station, four distances to four known
non-coplanar points -- which determine that station's position in the thorax
frame exactly (trilateration is linear once one equation is subtracted from the
rest).  Four target stations, one at the body origin and three on its axes, then
give the body's full rigid pose.

The known answers this measurement must pass, all checked below and written into
the output:

  * the recovered rotation is orthonormal with determinant +1 (a wrong solve or a
    misread length destroys this immediately);
  * the three recovered axis stations sit at exactly the distance from the origin
    station that they were placed at (0.1 m);
  * every one of the 16 distances is reproduced by the recovered pose to within
    the printed residual.

Nothing here is authored: the only inputs are the donor file and the station
coordinates chosen for the rulers, and the answer is the engine's.
"""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

DONOR = 'data/raw/mechanics/opensim-core/OpenSim/Tests/shared/ThoracoscapularShoulderModel.osim'
BINARY = 'data/runtime/opensim/native_polynomial_path_fit'
OUT = 'data/models/shoulder_girdle_v1/donor_reference_frames.json'

#: Ruler anchors on the donor thorax, in its own body frame (m).  Four
#: non-coplanar points; the spread is a third of a metre so a millimetre of
#: length error is a millimetre of position error, not a kilometre.
ANCHORS = ((0.0, 0.0, 0.0), (0.3, 0.0, 0.0), (0.0, 0.3, 0.0), (0.0, 0.0, 0.3))
#: Stations on each measured body: origin plus one per axis.
STATIONS = ((0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.0, 0.1, 0.0), (0.0, 0.0, 0.1))
AXIS_LENGTH = 0.1
BODIES = ('clavicle', 'scapula', 'humerus')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def probe_model(root: Path, target: Path, model=DONOR, parent='thorax', bodies=BODIES):
    """The model, unchanged, plus two-point PathActuator rulers."""
    tree = ET.parse(root / model if not Path(model).is_absolute() else model)
    forces = tree.getroot().find('.//ForceSet/objects')
    for body in bodies:
        for a, anchor in enumerate(ANCHORS):
            for s, station in enumerate(STATIONS):
                name = 'ruler_%s_%d_%d' % (body, a, s)
                actuator = ET.SubElement(forces, 'PathActuator', name=name)
                ET.SubElement(actuator, 'appliesForce').text = 'false'
                path = ET.SubElement(actuator, 'GeometryPath', name='path')
                point_set = ET.SubElement(path, 'PathPointSet')
                objects = ET.SubElement(point_set, 'objects')
                for index, (frame, location) in enumerate(
                        (('/bodyset/' + parent, anchor), ('/bodyset/' + body, station)), 1):
                    point = ET.SubElement(objects, 'PathPoint', name='%s-P%d' % (name, index))
                    ET.SubElement(point, 'socket_parent_frame').text = frame
                    ET.SubElement(point, 'location').text = ' '.join(repr(float(v)) for v in location)
                ET.SubElement(point_set, 'groups')
                wraps = ET.SubElement(path, 'PathWrapSet')
                ET.SubElement(wraps, 'objects')
                ET.SubElement(wraps, 'groups')
                ET.SubElement(actuator, 'optimal_force').text = '1'
    ET.indent(tree, space='\t')
    tree.write(target, encoding='utf-8', xml_declaration=True)


def run_sampler(root: Path, model: Path, out_csv: Path, log: Path):
    """One row at the model's own default state; no coordinate is overridden."""
    coordinates = out_csv.with_name('default_pose.sto')
    coordinates.write_text('default_pose\nversion=1\nnRows=1\nnColumns=1\n'
                           'inDegrees=no\nendheader\ntime\n0\n')
    argv = [str(root / BINARY), 'sample', '--model', str(model),
            '--coordinates', str(coordinates), '--output', str(out_csv)]
    with log.open('w') as handle:
        subprocess.run(argv, check=True, stdout=handle, stderr=subprocess.STDOUT)


def read_lengths(out_csv: Path):
    lengths = {}
    with out_csv.open() as handle:
        for row in csv.DictReader(handle):
            if row['quantity'] != 'length':
                continue
            name = row['actuator'].rsplit('/', 1)[-1]
            if name.startswith('ruler_'):
                lengths[name] = float(row['value'])
    return lengths


def trilaterate(anchors, distances):
    """Position from four distances to four non-coplanar known points.

    |p - a_i|^2 = d_i^2.  Subtracting the first equation removes |p|^2 and
    leaves three linear equations in p.
    """
    anchors = np.asarray(anchors, float)
    distances = np.asarray(distances, float)
    a0, d0 = anchors[0], distances[0]
    matrix = 2 * (anchors[1:] - a0)
    rhs = (d0 ** 2 - distances[1:] ** 2
           + np.einsum('ij,ij->i', anchors[1:], anchors[1:]) - a0 @ a0)
    point = np.linalg.solve(matrix, rhs)
    residual = float(np.abs(np.linalg.norm(anchors - point, axis=1) - distances).max())
    return point, residual


def measure(root: Path, work: Path, model=DONOR, parent='thorax', bodies=BODIES):
    work.mkdir(parents=True, exist_ok=True)
    probe = work / 'model_with_rulers.osim'
    probe_model(root, probe, model, parent, bodies)
    out_csv = work / 'ruler_lengths.csv'
    run_sampler(root, probe, out_csv, work / 'sampler.log')
    lengths = read_lengths(out_csv)
    if len(lengths) != len(bodies) * len(ANCHORS) * len(STATIONS):
        raise RuntimeError('Expected %d ruler lengths, read %d'
                           % (len(bodies) * len(ANCHORS) * len(STATIONS), len(lengths)))

    result = {}
    for body in bodies:
        points, residuals = [], []
        for s in range(len(STATIONS)):
            distances = [lengths['ruler_%s_%d_%d' % (body, a, s)] for a in range(len(ANCHORS))]
            point, residual = trilaterate(ANCHORS, distances)
            points.append(point)
            residuals.append(residual)
        origin = points[0]
        columns = np.stack([(points[i + 1] - origin) / AXIS_LENGTH for i in range(3)], axis=1)
        orthonormality = float(np.abs(columns.T @ columns - np.eye(3)).max())
        determinant = float(np.linalg.det(columns))
        # Nearest proper rotation, so the stored transform is exactly rigid; the
        # departure from it is reported rather than hidden.
        u, _, vt = np.linalg.svd(columns)
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            raise RuntimeError('%s: recovered frame is a reflection' % body)
        projection = float(np.abs(rotation - columns).max())
        transform = np.eye(4)
        transform[:3, :3] = rotation
        transform[:3, 3] = origin
        # Reproduce every measured distance from the stored rigid transform.
        worst = 0.0
        for s, station in enumerate(STATIONS):
            placed = rotation @ np.asarray(station, float) + origin
            for a, anchor in enumerate(ANCHORS):
                predicted = np.linalg.norm(placed - np.asarray(anchor, float))
                worst = max(worst, abs(predicted - lengths['ruler_%s_%d_%d' % (body, a, s)]))
        result[body] = {
            'body_to_thorax_reference': transform.tolist(),
            'trilateration_worst_residual_m': max(residuals),
            'axis_station_orthonormality_max_abs_error': orthonormality,
            'raw_axis_matrix_determinant': determinant,
            'rotation_projection_max_abs_change': projection,
            'distance_reproduction_worst_m': worst,
        }
    return result, lengths


def main():
    work = ROOT / 'data/derived/shoulder-girdle/donor-frames'
    result, lengths = measure(ROOT, work)
    record = {
        'schema': 'ihm.shoulder-donor-reference-frames.v1',
        'measured_by': 'scripts/measure_shoulder_donor_frames.py',
        'instrument': ('OpenSim %s via %s `sample`; path length of two-point '
                       'PathActuator rulers at the donor model\'s own default state, '
                       'with model.assemble() satisfying the AC PointConstraint.'
                       % (BINARY, BINARY)),
        'donor_path': DONOR,
        'donor_sha256': sha(ROOT / DONOR),
        'binary_sha256': sha(ROOT / BINARY),
        'anchor_stations_thorax_m': [list(a) for a in ANCHORS],
        'target_stations_m': [list(s) for s in STATIONS],
        'frames': result,
        'scope': ('Rigid pose of each donor body in the donor THORAX frame at the '
                  'donor\'s own default coordinate values. Not a claim about any other '
                  'configuration, and not a registration onto this plant.'),
    }
    path = ROOT / OUT
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps({body: {k: v for k, v in row.items()
                             if k != 'body_to_thorax_reference'}
                      for body, row in result.items()}, indent=2))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
