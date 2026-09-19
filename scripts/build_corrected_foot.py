"""Build `data/models/corrected_foot_v1`: engineering_stance_v1 with the upstream
mtp geometry defect repaired, and nothing else changed.

THE DEFECT, and where every number below comes from.  `subject_walk_scaled.osim`
is `Rajagopal2016.osim` put through OpenSim's ScaleTool.  ScaleTool writes the
per-body factors it used into the model it emits, as `<Mesh><scale_factors>`, and
it multiplies every `PhysicalOffsetFrame` translation by its OWN PARENT BODY's
factors.  Checked here over all 44 joint offset frames the two models share: 20
are the zero vector, 1 is on `ground` and has no factors, 19 obey the rule to
1e-12, the two `walker_knee_*/femur_*_offset` frames were re-placed by the knee's
own spline rebuild, and **`mtp_{l,r}/calcn_{l,r}_offset` is byte-identical to the
unscaled generic value** -- 0.1788 -0.002 +-0.00108 m, the only joint offset in
the file that the scaling step missed.  The same two joints are the only ones
whose `orientation` was not carried over: Rajagopal declares
`-3.14159 +-0.619901 0` on BOTH mtp frames (2016 and LaiUhlrich2023 agree to the
byte), and the scaled model carries `0 0 0`, which turns the oblique metatarsal
break into a pure calcn-z hinge.

So the toe muscles' calcaneal via points were scaled by calcn_r's
(1.155608663, 1.1170848239, 1.3052527858) and the axis they are meant to straddle
was not moved with them.  In Rajagopal the last calcaneal point of each toe muscle
sits 5-17 mm PROXIMAL of the mtp axis; in the scaled model it sits 8-22 mm DISTAL
of it (`docs/FOOT_JOINTS.md`, `docs/WORKBENCH_AUTHENTICITY.md` 0.2a).  Rajagopal
ships `mtp_angle` LOCKED, so the defect could not show there; this repo's engine
runs it free.

THE FIX, derived and not chosen:

    translation <- Rajagopal2016's own translation * this model's own calcn
                   scale_factors, componentwise -- the rule the other 19 non-zero frames obey
    orientation <- Rajagopal2016's own orientation, on both frames of the joint
                   (angles are dimensionless; ScaleTool carried every other
                   joint's orientation over unchanged)

    mtp_r calcn offset  0.1788    -> 0.20662282894439998  m   (+27.8228 mm)
                       -0.002     -> -0.0022341696478     m   ( -0.2342 mm)
                       +0.00108   -> +0.001409673008664   m   ( +0.3297 mm)
    mtp_r orientation   0 0 0     -> -3.14159 +0.619901 0
    mtp_l orientation   0 0 0     -> -3.14159 -0.619901 0

Nothing in this file is typed: both halves are read out of
`Rajagopal2016.osim` and out of `engineering_stance_v1/model.osim` at build time,
and `--check` re-derives them.

CONSEQUENCE THAT IS NOT A MOMENT ARM, and which the runner must control for: the
mtp axis carries the `toes_{l,r}` BODY with it, so the toes, their two source
contact spheres (`lateralToe`, `medialToe`) and the toe muscles' toes-frame points
all move 27.82 mm distally.  The corrected foot is a longer foot.  That is why
`registration.json` here (corrected geometry, shipped fitted paths) exists beside
`registration_toe_paths.json`: it is the arm that carries the geometry change with
NO new mtp moment arm, and the 2x2 in `docs/FOOT_GEOMETRY.md` separates the two.

`engineering_stance_v1` is read and never written.  It is the identified plant and
`docs/FOOT_JOINTS.md`'s recommendation is explicitly to leave it alone.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_corrected_foot --check
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

BASE = 'data/models/engineering_stance_v1'
OUT = 'data/models/corrected_foot_v1'
RAJAGOPAL = 'data/raw/anatomy/opensim-models/source/Models/Rajagopal/Rajagopal2016.osim'
LAI = 'data/raw/anatomy/opensim-models/source/Models/Rajagopal/RajagopalLaiUhlrich2023.osim'
UPSTREAM = ('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/'
            'subject_walk_scaled.osim')
PATHS = ('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/'
         'subject_walk_scaled_FunctionBasedPathSet.xml')

#: The four muscles per side whose own GeometryPath crosses mtp.  Same tuple as
#: scripts/build_foot_joint_arms.TOE_MUSCLES, restated so this builder does not
#: import a module another worker is editing.
TOE_MUSCLES = tuple('%s_%s' % (m, s) for s in ('r', 'l') for m in ('edl', 'ehl', 'fdl', 'fhl'))


def sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(element) -> str:
    return ET.canonicalize(ET.tostring(element, encoding='unicode'), strip_text=True)


def triple(text):
    return tuple(float(x) for x in (text or '0 0 0').split())


def joint_frames(root, name):
    """{frame name: element} for one joint's PhysicalOffsetFrames."""
    for joint in root.iter():
        if joint.get('name') == name and joint.tag.endswith('Joint'):
            return {f.get('name'): f for f in joint.iter('PhysicalOffsetFrame')}
    raise KeyError(name)


def mesh_scale(root, body):
    """The ScaleTool factors this model itself records for one body."""
    for b in root.iter('Body'):
        if b.get('name') == body:
            for mesh in b.iter('Mesh'):
                return triple(mesh.findtext('scale_factors'))
    raise KeyError(body)


def derive(base_root, rajagopal_root, side):
    """The corrected calcn-side offset for mtp_<side>, from the two source files."""
    generic = joint_frames(rajagopal_root, 'mtp_%s' % side)['calcn_%s_offset' % side]
    scale = mesh_scale(base_root, 'calcn_%s' % side)
    generic_t = triple(generic.findtext('translation'))
    return {'generic_translation_m': list(generic_t),
            'calcn_scale_factors': list(scale),
            'translation_m': [generic_t[i] * scale[i] for i in range(3)],
            'orientation_rad': list(triple(generic.findtext('orientation')))}


def fmt(values):
    return ' '.join(repr(float(v)) for v in values)


def write_tree(tree, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(tree, space='\t')
    tree.write(path, encoding='utf-8', xml_declaration=True)


def corrected_model(root: Path):
    tree = ET.parse(root / BASE / 'model.osim')
    model = tree.getroot().find('Model')
    generic = ET.parse(root / RAJAGOPAL).getroot()
    record = {}
    for side in ('r', 'l'):
        fix = derive(model, generic, side)
        frames = joint_frames(model, 'mtp_%s' % side)
        calcn = frames['calcn_%s_offset' % side]
        fix['was_translation_m'] = list(triple(calcn.findtext('translation')))
        fix['was_orientation_rad'] = list(triple(calcn.findtext('orientation')))
        fix['delta_translation_mm'] = [1000.0 * (fix['translation_m'][i] - fix['was_translation_m'][i])
                                       for i in range(3)]
        calcn.find('translation').text = fmt(fix['translation_m'])
        # BOTH frames take the orientation, so the calcn->toes transform at q = 0 is
        # unchanged apart from the translation, and only the hinge AXIS moves.
        for frame in frames.values():
            frame.find('orientation').text = fmt(fix['orientation_rad'])
        record['mtp_%s' % side] = fix
    return tree, record


def filtered_paths(root: Path, destination: Path, drop):
    tree = ET.parse(root / PATHS)
    containers = [c for c in tree.getroot().iter('objects')
                  if c.find('FunctionBasedPath') is not None]
    assert len(containers) == 1
    dropped = []
    for element in list(containers[0]):
        if element.get('name', '').rsplit('/', 1)[-1] in drop:
            dropped.append(element.get('name').rsplit('/', 1)[-1])
            containers[0].remove(element)
    write_tree(tree, destination)
    return sorted(dropped)


def registration(root: Path, out: Path, name, model_path: Path, record, overrides=None):
    catalog_rel = BASE + '/catalog.json'
    body = {
        'schema': 'ihm.corrected-foot.v1',
        'arm': name,
        'model_path': str(model_path.relative_to(root)),
        'model_sha256': sha(model_path),
        'catalog_path': catalog_rel,
        'catalog_sha256': sha(root / catalog_rel),
        'sources': {p: sha(root / p) for p in
                    (BASE + '/model.osim', RAJAGOPAL, LAI, UPSTREAM, PATHS,
                     'scripts/build_corrected_foot.py')},
        'default_enabled': False,
        'native_acceptance_complete': False,
        'scope': ('An experimental variant for docs/FOOT_GEOMETRY.md: the upstream mtp '
                  'offset/axis defect repaired. No linearization, no stance acceptance, '
                  'not the identified plant.'),
        **record}
    if overrides:
        body['source_overrides'] = {k: {'path': str(v.relative_to(root)), 'sha256': sha(v)}
                                    for k, v in overrides.items()}
    path = out / {'corrected_fitted_paths': 'registration.json',
                  'base_toe_paths': 'registration_base_toe_paths.json',
                  'corrected_toe_paths': 'registration_toe_paths.json'}[name]
    path.write_text(json.dumps(body, indent=2, sort_keys=True) + '\n')
    return path


def build(root: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    tree, fixes = corrected_model(root)
    model_path = out / 'model.osim'
    write_tree(tree, model_path)
    written = {'model': model_path}

    shared = {'derived_from': 'data/models/engineering_stance_v1/model.osim, mtp_{r,l} only',
              'mtp_correction': fixes,
              'derivation': ('Rajagopal2016 mtp calcn-offset translation x this model\'s own '
                             'calcn <Mesh><scale_factors>, componentwise; orientation copied '
                             'from Rajagopal2016 onto both frames of the joint. The rule is '
                             'the one the other 42 joint offset frames of '
                             'subject_walk_scaled.osim already obey.')}
    written['corrected_fitted_paths'] = registration(
        root, out, 'corrected_fitted_paths', model_path,
        dict(shared, role=('corrected geometry, shipped fitted path set: the geometry and '
                           'contact change with NO mtp moment arm')))

    paths = out / 'muscle_paths.xml'
    dropped = filtered_paths(root, paths, set(TOE_MUSCLES))

    # The shipped-geometry control arm of the 2x2, so the whole experiment can run
    # without rebuilding data/derived/foot-joint-arms, which another worker is using.
    # It selects the BASE model unchanged with the same filtered path set, i.e. it is
    # docs/FOOT_JOINTS.md's `base_toe_geometry_paths` reconstructed; the runner's I3
    # check requires it to reproduce that arm's recorded moment arms to 1e-6 m.
    written['base_toe_paths'] = registration(
        root, out, 'base_toe_paths', root / BASE / 'model.osim',
        {'derived_from': BASE + '/model.osim, unchanged',
         'role': ('shipped geometry, toe muscles on their own GeometryPaths: '
                  "docs/FOOT_JOINTS.md's base_toe_geometry_paths arm, rebuilt here"),
         'geometry_path_muscles': dropped},
        {'subject_walk_scaled_FunctionBasedPathSet.xml': paths})

    written['corrected_toe_paths'] = registration(
        root, out, 'corrected_toe_paths', model_path,
        dict(shared, role=('corrected geometry, and edl/ehl/fdl/fhl on their own '
                           'GeometryPaths so the toe muscles can see mtp'),
             geometry_path_muscles=dropped),
        {'subject_walk_scaled_FunctionBasedPathSet.xml': paths})
    return written


# --------------------------------------------------------------------- checks
def check(root: Path):
    """Static known answers.  Each can fail for the reason it exists."""
    out = root / OUT
    build(root, out)

    def digest(where: Path):
        prefix = str(where.relative_to(root))
        found = {}
        for p in sorted(where.rglob('*')):
            if p.is_file():
                data = p.read_bytes()
                if p.name.startswith('registration'):
                    data = data.replace(prefix.encode(), b'<OUT>')
                found[str(p.relative_to(where))] = hashlib.sha256(data).hexdigest()
        return found

    first = digest(out)
    with tempfile.TemporaryDirectory(dir=root / 'data/derived') as scratch:
        build(root, Path(scratch) / 'again')
        again = digest(Path(scratch) / 'again')
    assert first == again, 'C1 the builder is not a function of its inputs'

    base = ET.parse(root / BASE / 'model.osim').getroot().find('Model')
    fixed = ET.parse(out / 'model.osim').getroot().find('Model')
    generic = ET.parse(root / RAJAGOPAL).getroot()
    lai = ET.parse(root / LAI).getroot()
    upstream = ET.parse(root / UPSTREAM).getroot()

    # C2  the two Rajagopal releases declare the SAME mtp frames, so "Rajagopal's
    #     value" is not a choice between two sources
    for side in ('r', 'l'):
        a = joint_frames(generic, 'mtp_%s' % side)
        b = joint_frames(lai, 'mtp_%s' % side)
        for k in a:
            assert triple(a[k].findtext('translation')) == triple(b[k].findtext('translation')), (side, k)
            assert triple(a[k].findtext('orientation')) == triple(b[k].findtext('orientation')), (side, k)

    # C3  THE DERIVATION.  Every joint offset frame of the upstream scaled model is
    #     its generic translation times its own parent body's recorded scale factors
    #     -- except the two mtp calcn frames, which are the unscaled generic value.
    scaled, unscaled, other, zero = [], [], [], []
    for joint in generic.iter():
        if not (joint.get('name') and joint.tag.endswith('Joint')):
            continue
        try:
            up = joint_frames(upstream, joint.get('name'))
        except KeyError:
            continue
        for name, frame in joint_frames(generic, joint.get('name')).items():
            if name not in up:
                continue
            body = (frame.findtext('socket_parent') or '').strip().split('/')[-1]
            try:
                s = mesh_scale(upstream, body)
            except KeyError:
                continue
            g, u = triple(frame.findtext('translation')), triple(up[name].findtext('translation'))
            key = (joint.get('name'), name)
            if g == (0.0, 0.0, 0.0) and u == (0.0, 0.0, 0.0):
                zero.append(key)
            elif max(abs(g[i] * s[i] - u[i]) for i in range(3)) < 1e-12:
                scaled.append(key)
            elif g == u:
                unscaled.append(key)
            else:
                other.append(key)
    assert sorted(unscaled) == [('mtp_l', 'calcn_l_offset'), ('mtp_r', 'calcn_r_offset')], unscaled
    assert sorted(other) == [('walker_knee_l', 'femur_l_offset'),
                             ('walker_knee_r', 'femur_r_offset')], other
    # measured on this pair of files, not a target: 19 scaled + 20 zero + 2 mtp
    # + 2 walker_knee + 1 ground offset = the 44 shared frames
    assert (len(scaled), len(zero)) == (19, 20), (len(scaled), len(zero))

    # C4  only mtp_{r,l} differs from the base model, element for element
    def elements(model, tag):
        return {e.get('name'): canonical(e) for e in model.find(tag)}
    for tag in ('BodySet/objects', 'ForceSet/objects', 'ConstraintSet/objects'):
        if base.find(tag) is not None:
            assert elements(base, tag) == elements(fixed, tag), tag
    bj, fj = elements(base, 'JointSet/objects'), elements(fixed, 'JointSet/objects')
    assert set(bj) == set(fj)
    assert {k for k in bj if bj[k] != fj[k]} == {'mtp_r', 'mtp_l'}, \
        sorted(k for k in bj if bj[k] != fj[k])

    # C5  the corrected values ARE the derived ones, and the rest pose moves by a
    #     pure translation: both frames of each joint carry the same orientation,
    #     so calcn->toes at q = 0 is unchanged apart from that translation
    for side in ('r', 'l'):
        want = derive(base, generic, side)
        frames = joint_frames(fixed, 'mtp_%s' % side)
        got = frames['calcn_%s_offset' % side]
        assert triple(got.findtext('translation')) == tuple(want['translation_m']), side
        for frame in frames.values():
            assert triple(frame.findtext('orientation')) == tuple(want['orientation_rad']), side
        was = triple(joint_frames(base, 'mtp_%s' % side)['calcn_%s_offset' % side]
                     .findtext('translation'))
        assert abs(1000 * (want['translation_m'][0] - was[0]) - 27.8228289444) < 1e-6, side

    # C6  the toe path set: exactly the eight removed, the other 72 byte-equal
    shipped = {e.get('name'): canonical(e) for e in ET.parse(root / PATHS).getroot().iter('FunctionBasedPath')}
    kept = {e.get('name'): canonical(e) for e in ET.parse(out / 'muscle_paths.xml').getroot().iter('FunctionBasedPath')}
    assert len(shipped) == 80 and len(kept) == 72, (len(shipped), len(kept))
    assert all(kept[k] == shipped[k] for k in kept)
    assert sorted(json.loads((out / 'registration_toe_paths.json').read_text())
                  ['geometry_path_muscles']) == sorted(TOE_MUSCLES)

    # C7  the base plant is untouched
    assert sha(root / BASE / 'model.osim') == \
        json.loads((out / 'registration.json').read_text())['sources'][BASE + '/model.osim']

    print('static checks C1-C7 pass (%d files, builder is a function)' % len(first))
    print('  C3: 19 non-zero upstream joint offsets scaled exactly, 20 zero, '
          '2 walker_knee re-placed, 2 mtp calcn offsets UNSCALED -- the defect')
    for side in ('r', 'l'):
        w = derive(base, generic, side)
        print('  mtp_%s  %s -> %s   orient -> %s' % (
            side, w['generic_translation_m'], [round(v, 12) for v in w['translation_m']],
            w['orientation_rad']))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check:
        check(ROOT)
    else:
        for name, path in build(ROOT, ROOT / OUT).items():
            print(name, Path(path).relative_to(ROOT))


if __name__ == '__main__':
    main()
