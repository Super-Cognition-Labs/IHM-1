"""Build the plants that separate WHY the ankle collapses on articulated_spine_v1,
and the plants that price the base plant's two undriven toe hinges.

Pre-registration, arms and reading rules: docs/FOOT_JOINTS.md.  This script
only writes model files; it integrates nothing.

Q1 -- the ankle.  articulated_spine_v1 differs from engineering_stance_v1 in
three separable ways, and each becomes a factor:

  T  the torso repartitioned into torso/thorax/cervical/head, with the thoracic,
     neck and atlanto-occipital joints (nine coordinates, stops, damping)
  W  radius_hand_{l,r} WeldJoint -> UniversalJoint (+ stops, damping)
  S  subtalar_{l,r}    WeldJoint -> PinJoint       (+ stop, damping)

Arm `t{T}w{W}s{S}` carries exactly the factors set to 1.  T=0 arms are built
FROM THE BASE with the variant builder's own `unweld`/`limit_force`/
`damping_force`, so a W or S joint is element-for-element the variant's.  T=1
arms are built FROM THE COMMITTED VARIANT by re-welding W and/or S with the
joint's own frames and removing that coordinate's stop and damping.  t1w1s1 is
not rebuilt: it is the committed registration.json.  t0w0s0 is rebuilt through
the round trip base -> unweld(S, W) -> reweld(S, W) and must equal the base
(checked below, and again in the engine by the runner).

`tweld` is T with its three new joints WELDED as well: the repartitioned masses
and the three extra bodies (and therefore their supine contact balls, and the
torso's larger ball) with no new articulation.  It splits T into "mass and
contact redistribution" against "spine articulation".

Q2 -- the toes (base plant only; engineering_stance_v1 is read, never written):

  base_mtp_welded           mtp_{l,r} PinJoint -> WeldJoint at q=0, as upstream
                            did before fitting the paths
                            (exampleMocoInverse.cpp:51).  The source passive
                            set's two `mtp_angle` terms (-25*q - 2*qdot) are
                            removed through a source override, because the
                            coordinate no longer exists.
  base_toe_geometry_paths   the model unchanged; the shipped FunctionBasedPathSet
                            with edl/ehl/fdl/fhl (both sides) REMOVED, so those
                            eight run on the model's own GeometryPaths -- the
                            route 1e29ee3 took for the subtalar.

Everything is written under data/derived/foot-joint-arms/ (gitignored): the
builder is deterministic, and `--check` builds twice and compares bytes.

    OPENBLAS_NUM_THREADS=1 .venv/bin/python -m scripts.build_foot_joint_arms --check
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import shutil
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scripts.build_articulated_spine import (  # noqa: E402
    DAMPING_LOWER, DAMPING_UPPER, NEW_JOINTS, SUBTALAR_RANGE, WRIST_RANGE,
    damping_force, find_joint, limit_force, sub, unweld,
)

BASE = 'data/models/engineering_stance_v1'
VARIANT = 'data/models/articulated_spine_v1'
OUT = 'data/derived/foot-joint-arms'
SOURCE_DIR = 'data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking'
PASSIVE = SOURCE_DIR + '/subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml'
PATHS = SOURCE_DIR + '/subject_walk_scaled_FunctionBasedPathSet.xml'

#: The four muscles per side whose own GeometryPath crosses mtp
#: (verify_articulated_spine.MTP_CROSSING, measured there by M1).
TOE_MUSCLES = tuple('%s_%s' % (m, s) for s in ('r', 'l') for m in ('edl', 'ehl', 'fdl', 'fhl'))

FACTORIAL = tuple('t%dw%ds%d' % (t, w, s) for t in (0, 1) for w in (0, 1) for s in (0, 1))


def sha(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(element) -> str:
    """Whitespace-insensitive XML text, for element equality."""
    return ET.canonicalize(ET.tostring(element, encoding='unicode'), strip_text=True)


def reweld(jointset, forceset, name):
    """Real joint -> WeldJoint with the joint's OWN frames, so the rest pose is
    the joint at q = 0.  Removes every force element that names one of the
    joint's coordinates and returns (coordinates, removed force names)."""
    index, joint = find_joint(jointset, name)
    names = [c.get('name') for c in joint.find('coordinates')]
    weld = ET.Element('WeldJoint', name=name)
    sub(weld, 'socket_parent_frame', joint.findtext('socket_parent_frame'))
    sub(weld, 'socket_child_frame', joint.findtext('socket_child_frame'))
    weld.append(joint.find('frames'))
    jointset[index] = weld
    removed = []
    for force in list(forceset):
        if force.findtext('coordinate') in names:
            removed.append(force.get('name'))
            forceset.remove(force)
    return names, removed


def parts(tree):
    model = tree.getroot().find('Model')
    return model, model.find('JointSet/objects'), model.find('ForceSet/objects')


def unweld_from_base(jointset, forceset, *, wrists, subtalar):
    """The variant builder's own unweld, in the builder's own per-side order."""
    for side in ('r', 'l'):
        if wrists:
            unweld(jointset, 'radius_hand_%s' % side, 'UniversalJoint',
                   [('wrist_flex_%s' % side, WRIST_RANGE['wrist_flex']),
                    ('wrist_dev_%s' % side, WRIST_RANGE['wrist_dev'])])
        if subtalar:
            unweld(jointset, 'subtalar_%s' % side, 'PinJoint',
                   [('subtalar_angle_%s' % side, SUBTALAR_RANGE)])
    stops = []
    for side in ('r', 'l'):
        if wrists:
            stops.append(('wrist_flex_%s' % side, WRIST_RANGE['wrist_flex'], DAMPING_UPPER))
            stops.append(('wrist_dev_%s' % side, WRIST_RANGE['wrist_dev'], DAMPING_UPPER))
        if subtalar:
            stops.append(('subtalar_angle_%s' % side, SUBTALAR_RANGE, DAMPING_LOWER))
    for name, bounds, damping in stops:
        forceset.append(limit_force(name, bounds))
        forceset.append(damping_force(name, damping))


def write_tree(tree, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(tree, space='\t')
    tree.write(path, encoding='utf-8', xml_declaration=True)


def factorial_model(root: Path, arm: str):
    t, w, s = (int(arm[i]) for i in (1, 3, 5))
    record = {'factors': {'T': t, 'W': w, 'S': s}}
    if t == 0:
        tree = ET.parse(root / BASE / 'model.osim')
        model, jointset, forceset = parts(tree)
        if arm == 't0w0s0':
            # the round trip: the variant builder's unweld, then this reweld
            unweld_from_base(jointset, forceset, wrists=True, subtalar=True)
            for side in ('r', 'l'):
                reweld(jointset, forceset, 'radius_hand_%s' % side)
                reweld(jointset, forceset, 'subtalar_%s' % side)
            record['derived_from'] = 'engineering_stance_v1 via unweld(S,W) then reweld(S,W)'
        else:
            unweld_from_base(jointset, forceset, wrists=bool(w), subtalar=bool(s))
            record['derived_from'] = 'engineering_stance_v1 + the variant builder\'s unweld'
    else:
        tree = ET.parse(root / VARIANT / 'model.osim')
        model, jointset, forceset = parts(tree)
        removed = []
        for side in ('r', 'l'):
            if not w:
                removed += reweld(jointset, forceset, 'radius_hand_%s' % side)[1]
            if not s:
                removed += reweld(jointset, forceset, 'subtalar_%s' % side)[1]
        record['derived_from'] = 'articulated_spine_v1/model.osim with re-welds'
        record['removed_forces'] = removed
    return tree, record


def tweld_model(root: Path):
    tree = ET.parse(root / VARIANT / 'model.osim')
    model, jointset, forceset = parts(tree)
    removed = []
    for name in NEW_JOINTS:
        removed += reweld(jointset, forceset, name)[1]
    for side in ('r', 'l'):
        removed += reweld(jointset, forceset, 'radius_hand_%s' % side)[1]
        removed += reweld(jointset, forceset, 'subtalar_%s' % side)[1]
    return tree, {'factors': {'T': 'masses/bodies only, joints welded', 'W': 0, 'S': 0},
                  'derived_from': 'articulated_spine_v1/model.osim with all 7 new joints re-welded',
                  'removed_forces': removed}


def mtp_welded_model(root: Path):
    tree = ET.parse(root / BASE / 'model.osim')
    model, jointset, forceset = parts(tree)
    removed = []
    for side in ('r', 'l'):
        removed += reweld(jointset, forceset, 'mtp_%s' % side)[1]
    return tree, {'derived_from': 'engineering_stance_v1 with mtp_{r,l} welded at q=0',
                  'removed_forces_in_model': removed}


def filtered_passive(root: Path, destination: Path, drop_coordinates):
    tree = ET.parse(root / PASSIVE)
    containers = [c for c in tree.getroot().iter('objects')
                  if c.find('ExpressionBasedCoordinateForce') is not None]
    assert len(containers) == 1
    dropped = []
    for element in list(containers[0]):
        if element.findtext('coordinate') in drop_coordinates:
            dropped.append(element.get('name') or element.findtext('coordinate'))
            containers[0].remove(element)
    write_tree(tree, destination)
    return dropped


def filtered_paths(root: Path, destination: Path, drop_muscles):
    tree = ET.parse(root / PATHS)
    containers = [c for c in tree.getroot().iter('objects')
                  if c.find('FunctionBasedPath') is not None]
    assert len(containers) == 1
    dropped = []
    for element in list(containers[0]):
        if element.get('name', '').rsplit('/', 1)[-1] in drop_muscles:
            dropped.append(element.get('name').rsplit('/', 1)[-1])
            containers[0].remove(element)
    write_tree(tree, destination)
    return sorted(dropped)


def registration(root: Path, out: Path, arm: str, model_path: Path, record, overrides=None):
    model_rel = str(model_path.relative_to(root))
    catalog_rel = BASE + '/catalog.json'
    sources = {BASE + '/model.osim': sha(root / BASE / 'model.osim'),
               VARIANT + '/model.osim': sha(root / VARIANT / 'model.osim'),
               'scripts/build_foot_joint_arms.py': sha(root / 'scripts/build_foot_joint_arms.py'),
               'scripts/build_articulated_spine.py': sha(root / 'scripts/build_articulated_spine.py'),
               model_rel: sha(model_path), catalog_rel: sha(root / catalog_rel)}
    body = {'schema': 'ihm.foot-joint-arm.v1', 'arm': arm,
            'model_path': model_rel, 'model_sha256': sha(model_path),
            'catalog_path': catalog_rel, 'catalog_sha256': sha(root / catalog_rel),
            'sources': sources, 'default_enabled': False,
            'native_acceptance_complete': False,
            'scope': ('An experimental arm for docs/FOOT_JOINTS.md. Not a plant for '
                      'any other use.'),
            **record}
    if overrides:
        body['source_overrides'] = {name: {'path': str(p.relative_to(root)), 'sha256': sha(p)}
                                    for name, p in overrides.items()}
    path = out / arm / 'registration.json'
    path.write_text(json.dumps(body, indent=2, sort_keys=True) + '\n')
    return path


def build(root: Path, out: Path):
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    written = {}
    for arm in FACTORIAL:
        if arm == 't1w1s1':
            continue
        tree, record = factorial_model(root, arm)
        model_path = out / arm / 'model.osim'
        write_tree(tree, model_path)
        written[arm] = registration(root, out, arm, model_path, record)
    tree, record = tweld_model(root)
    write_tree(tree, out / 'tweld' / 'model.osim')
    written['tweld'] = registration(root, out, 'tweld', out / 'tweld' / 'model.osim', record)

    tree, record = mtp_welded_model(root)
    model_path = out / 'base_mtp_welded' / 'model.osim'
    write_tree(tree, model_path)
    passive = out / 'base_mtp_welded' / 'passive.xml'
    record['removed_passive_terms'] = filtered_passive(root, passive, {'mtp_angle_r', 'mtp_angle_l'})
    written['base_mtp_welded'] = registration(
        root, out, 'base_mtp_welded', model_path, record,
        {'subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml': passive})

    model_path = out / 'base_toe_geometry_paths' / 'model.osim'
    model_path.parent.mkdir(parents=True)
    shutil.copyfile(root / BASE / 'model.osim', model_path)
    paths = out / 'base_toe_geometry_paths' / 'muscle_paths.xml'
    dropped = filtered_paths(root, paths, set(TOE_MUSCLES))
    written['base_toe_geometry_paths'] = registration(
        root, out, 'base_toe_geometry_paths', model_path,
        {'derived_from': 'engineering_stance_v1 model unchanged; shipped path set minus the '
                         '8 mtp-crossing muscles, which run on their own GeometryPaths',
         'geometry_path_muscles': dropped},
        {'subject_walk_scaled_FunctionBasedPathSet.xml': paths})
    return written


# --------------------------------------------------------------------- checks
def check(root: Path):
    """Static known answers.  Each can fail for the reason it exists."""
    out = root / OUT
    build(root, out)
    def digest(where: Path):
        # registration.json names its files by path, and the second build lives
        # elsewhere, so it is compared with that one prefix substituted
        prefix = str(where.relative_to(root))
        found = {}
        for p in sorted(where.rglob('*')):
            if p.is_file():
                data = p.read_bytes()
                if p.name == 'registration.json':
                    data = data.replace(prefix.encode(), b'<OUT>')
                found[str(p.relative_to(where))] = hashlib.sha256(data).hexdigest()
        return found

    first = digest(out)
    with tempfile.TemporaryDirectory(dir=root / 'data/derived') as scratch:
        build(root, Path(scratch) / 'again')
        again = digest(Path(scratch) / 'again')
    assert first == again, 'the arm builder is not a function of its inputs'

    def model_of(path):
        return ET.parse(path).getroot().find('Model')

    base = model_of(root / BASE / 'model.osim')
    variant = model_of(root / VARIANT / 'model.osim')
    # 1. the round trip reproduces the base, element for element
    assert canonical(model_of(out / 't0w0s0/model.osim')) == canonical(base), \
        'unweld -> reweld does not return the base model'

    def joints(model):
        return {j.get('name'): canonical(j) for j in model.find('JointSet/objects')}

    def forces(model):
        return {f.get('name'): canonical(f) for f in model.find('ForceSet/objects')
                if f.findtext('coordinate') is not None}

    def bodies(model):
        return {b.get('name'): canonical(b) for b in model.find('BodySet/objects')}

    vj, vf = joints(variant), forces(variant)
    bj, bb = joints(base), bodies(base)
    wrist = {'radius_hand_r', 'radius_hand_l'}
    subtalar = {'subtalar_r', 'subtalar_l'}
    wrist_c = {'wrist_flex_r', 'wrist_dev_r', 'wrist_flex_l', 'wrist_dev_l'}
    subtalar_c = {'subtalar_angle_r', 'subtalar_angle_l'}
    for arm in FACTORIAL:
        if arm == 't1w1s1':
            continue
        t, w, s = (int(arm[i]) for i in (1, 3, 5))
        m = model_of(out / arm / 'model.osim')
        aj, af = joints(m), forces(m)
        free = (wrist if w else set()) | (subtalar if s else set())
        # 2. every freed joint is the variant's joint, element for element; every
        #    welded one is the base's
        for name in wrist | subtalar:
            expected = vj[name] if name in free else bj[name]
            assert aj[name] == expected, (arm, name)
        # 3. the stops and damping present are exactly the freed coordinates'
        coords = (wrist_c if w else set()) | (subtalar_c if s else set())
        foot_hand = {k for k, v in af.items() if any(c in k for c in wrist_c | subtalar_c)}
        assert {k.split('_', 2)[2] for k in foot_hand} == coords, (arm, sorted(foot_hand))
        for k in foot_hand:
            assert af[k] == vf[k], (arm, k)
        # 4. T=0 arms keep the base bodies; T=1 arms the variant's
        mb = bodies(m)
        assert mb == (bodies(variant) if t else bb), arm
        # 5. no reference to a welded coordinate survives anywhere in the file
        text = (out / arm / 'model.osim').read_text()
        for c in (wrist_c | subtalar_c) - coords:
            assert c not in text, (arm, c)
    tw = model_of(out / 'tweld/model.osim')
    assert bodies(tw) == bodies(variant)
    for name in NEW_JOINTS:
        assert tw.find("JointSet/objects/WeldJoint[@name='%s']" % name) is not None
    text = (out / 'tweld/model.osim').read_text()
    for c in ('thoracic_extension', 'neck_extension', 'head_rotation') + tuple(wrist_c | subtalar_c):
        assert c not in text, c
    # 6. mtp welded: no mtp coordinate in the model or in the passive override;
    #    everything else in the passive set is kept
    for c in ('mtp_angle_r', 'mtp_angle_l'):
        assert c not in (out / 'base_mtp_welded/model.osim').read_text()
        assert c not in (out / 'base_mtp_welded/passive.xml').read_text()
    shipped = ET.parse(root / PASSIVE).getroot()
    kept = ET.parse(out / 'base_mtp_welded/passive.xml').getroot()
    n_shipped = len(list(shipped.iter('ExpressionBasedCoordinateForce')))
    n_kept = len(list(kept.iter('ExpressionBasedCoordinateForce')))
    assert n_shipped - n_kept == 2, (n_shipped, n_kept)
    # 7. toe paths: exactly the eight removed, the other 72 byte-equal in content
    reg = json.loads((out / 'base_toe_geometry_paths/registration.json').read_text())
    assert sorted(reg['geometry_path_muscles']) == sorted(TOE_MUSCLES)
    shipped_paths = {e.get('name'): canonical(e) for e in ET.parse(root / PATHS).getroot().iter('FunctionBasedPath')}
    kept_paths = {e.get('name'): canonical(e) for e in ET.parse(out / 'base_toe_geometry_paths/muscle_paths.xml').getroot().iter('FunctionBasedPath')}
    assert len(shipped_paths) == 80 and len(kept_paths) == 72
    assert all(kept_paths[k] == shipped_paths[k] for k in kept_paths)
    print('static checks: 7 pass; builder is a function (%d files byte-identical)' % len(first))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check:
        check(ROOT)
    else:
        for arm, path in build(ROOT, ROOT / OUT).items():
            print(arm, path.relative_to(ROOT))


if __name__ == '__main__':
    main()
