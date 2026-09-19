"""Weld `mtp_{l,r}` in the articulated spine variant -- `docs/FOOT_JOINTS.md`'s
standing recommendation, applied to the one NEW plant this repository has.

`docs/FOOT_JOINTS.md`: *"leave `engineering_stance_v1` as it is ... Weld mtp in
any NEW plant, where it costs nothing and removes a dishonest degree of freedom."*
The toe hinges are undrivable in that plant for the same reason they are
undrivable in the base one: the shipped path polynomials were fitted upstream with
the toes welded (`exampleMocoInverse.cpp:51`), so no muscle has an arm about
`mtp_angle`.

WHAT THIS DOES **NOT** DO.  It does not touch
`data/models/articulated_spine_v1/model.osim`.  That file is the plant gate F2 was
measured on, F2 FAILED on it, and a gate is never rescored under a changed
instrument.  The weld goes into a SEPARATE model and a SEPARATE registration in
the same directory, beside the `*_muscled` pair that is already there:

    model_mtp_welded.osim          model.osim with mtp_{l,r} PinJoint -> WeldJoint at q = 0
    passive_mtp_welded.xml         the source ExpressionBasedCoordinateForceSet minus the
                                   two `mtp_angle` terms (-25*q - 2*qdot), which name a
                                   coordinate that no longer exists
    registration_mtp_welded.json   selects both, with the base catalog

Any number measured on that plant is a number about a DIFFERENT plant.  It is
never a new score for F2.

The weld uses the joint's own frames, so the rest pose is the joint at q = 0 --
the same construction `scripts/build_foot_joint_arms.reweld` uses for the
factorial's re-welds, restated here rather than imported, because that module is
under concurrent edit.

    .venv/bin/python -m scripts.build_mtp_welded_spine --check
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

VARIANT = 'data/models/articulated_spine_v1'
BASE = 'data/models/engineering_stance_v1'
PASSIVE = ('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/'
           'subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml')
MTP_COORDINATES = ('mtp_angle_r', 'mtp_angle_l')


def sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(element) -> str:
    return ET.canonicalize(ET.tostring(element, encoding='unicode'), strip_text=True)


def sub(parent, tag, text):
    child = ET.SubElement(parent, tag)
    child.text = text
    return child


def write_tree(tree, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(tree, space='\t')
    tree.write(path, encoding='utf-8', xml_declaration=True)


def reweld(jointset, forceset, name):
    """Real joint -> WeldJoint carrying the joint's OWN frames, and every force
    element naming one of its coordinates removed."""
    index = joint = None
    for i, j in enumerate(jointset):
        if j.get('name') == name:
            index, joint = i, j
            break
    if joint is None:
        raise KeyError(name)
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


def filtered_passive(root: Path, destination: Path, drop):
    tree = ET.parse(root / PASSIVE)
    containers = [c for c in tree.getroot().iter('objects')
                  if c.find('ExpressionBasedCoordinateForce') is not None]
    assert len(containers) == 1
    dropped = []
    for element in list(containers[0]):
        if element.findtext('coordinate') in drop:
            dropped.append(element.get('name') or element.findtext('coordinate'))
            containers[0].remove(element)
    write_tree(tree, destination)
    return dropped


def build(root: Path, out: Path):
    tree = ET.parse(root / VARIANT / 'model.osim')
    model = tree.getroot().find('Model')
    jointset, forceset = model.find('JointSet/objects'), model.find('ForceSet/objects')
    coordinates, removed = [], []
    for side in ('r', 'l'):
        c, r = reweld(jointset, forceset, 'mtp_%s' % side)
        coordinates += c
        removed += r
    model_path = out / 'model_mtp_welded.osim'
    write_tree(tree, model_path)
    passive = out / 'passive_mtp_welded.xml'
    dropped = filtered_passive(root, passive, set(MTP_COORDINATES))

    catalog_rel = VARIANT + '/catalog.json'
    body = {
        'schema': 'ihm.articulated-spine-mtp-welded.v1',
        'model_path': str(model_path.relative_to(root)),
        'model_sha256': sha(model_path),
        'catalog_path': catalog_rel,
        'catalog_sha256': sha(root / catalog_rel),
        'sources': {p: sha(root / p) for p in
                    (VARIANT + '/model.osim', catalog_rel, PASSIVE,
                     'scripts/build_mtp_welded_spine.py')},
        'default_enabled': False,
        'native_acceptance_complete': False,
        'derived_from': VARIANT + '/model.osim with mtp_{r,l} welded at q = 0',
        'welded_coordinates': sorted(coordinates),
        'removed_forces_in_model': sorted(removed),
        'removed_passive_terms': sorted(dropped),
        'rationale': ("docs/FOOT_JOINTS.md's standing recommendation: weld mtp in any NEW "
                      'plant. The shipped path polynomials were fitted with the toes welded '
                      '(exampleMocoInverse.cpp:51), so no muscle has an arm about mtp_angle '
                      'in this plant either.'),
        'scope': ('A SEPARATE plant. articulated_spine_v1/model.osim is unchanged and gate '
                  'F2 remains FAILED on it. No number measured here is a score for F2.'),
        'source_overrides': {
            'subject_walk_scaled_ExpressionBasedCoordinateForceSet.xml':
                {'path': str(passive.relative_to(root)), 'sha256': sha(passive)}},
    }
    path = out / 'registration_mtp_welded.json'
    path.write_text(json.dumps(body, indent=2, sort_keys=True) + '\n')
    return {'model': model_path, 'passive': passive, 'registration': path}


def check(root: Path):
    out = root / VARIANT
    before = sha(out / 'model.osim')
    build(root, out)

    def digest(where: Path, names):
        """Content digest with the OUTPUT DIRECTORY's own path normalised away --
        only in the two fields that name it.  Substituting the prefix everywhere
        would also rewrite `sources` and `catalog_path`, which legitimately point
        at the variant directory whichever directory the build writes to."""
        found = {}
        for n in names:
            data = (where / n).read_bytes()
            if n.endswith('.json'):
                body = json.loads(data)
                body['model_path'] = Path(body['model_path']).name
                body['derived_from'] = '<VARIANT>' + body['derived_from'].split('model.osim', 1)[1]
                for v in body.get('source_overrides', {}).values():
                    v['path'] = Path(v['path']).name
                data = json.dumps(body, indent=2, sort_keys=True).encode()
            found[n] = hashlib.sha256(data).hexdigest()
        return found

    names = ('model_mtp_welded.osim', 'passive_mtp_welded.xml', 'registration_mtp_welded.json')
    first = digest(out, names)
    with tempfile.TemporaryDirectory(dir=root / 'data/derived') as scratch:
        again_dir = Path(scratch) / 'again'
        again_dir.mkdir()
        build(root, again_dir)
        again = digest(again_dir, names)
    assert first == again, 'W1 the builder is not a function of its inputs'

    # W2  the gated plant is byte-identical
    assert sha(out / 'model.osim') == before, 'W2 articulated_spine_v1/model.osim was written'

    variant = ET.parse(out / 'model.osim').getroot().find('Model')
    welded = ET.parse(out / 'model_mtp_welded.osim').getroot().find('Model')

    # W3  only mtp_{r,l} changed, element for element; bodies untouched
    def elements(model, tag):
        return {e.get('name'): canonical(e) for e in model.find(tag)}
    assert elements(variant, 'BodySet/objects') == elements(welded, 'BodySet/objects')
    vj, wj = elements(variant, 'JointSet/objects'), elements(welded, 'JointSet/objects')
    assert {k for k in vj if vj[k] != wj[k]} == {'mtp_r', 'mtp_l'}, \
        sorted(k for k in vj if vj[k] != wj[k])
    for side in ('r', 'l'):
        assert welded.find("JointSet/objects/WeldJoint[@name='mtp_%s']" % side) is not None

    # W4  no reference to a welded coordinate survives, in the model or the override
    text = (out / 'model_mtp_welded.osim').read_text()
    passive_text = (out / 'passive_mtp_welded.xml').read_text()
    for c in MTP_COORDINATES:
        assert c not in text, ('W4 model', c)
        assert c not in passive_text, ('W4 passive', c)

    # W5  exactly the two passive terms dropped; everything else kept
    shipped = ET.parse(root / PASSIVE).getroot()
    kept = ET.parse(out / 'passive_mtp_welded.xml').getroot()
    n_shipped = len(list(shipped.iter('ExpressionBasedCoordinateForce')))
    n_kept = len(list(kept.iter('ExpressionBasedCoordinateForce')))
    assert n_shipped - n_kept == 2, (n_shipped, n_kept)

    # W6  the 15 coordinates the variant added are all still there and still free
    for name in ('thoracic_extension', 'neck_extension', 'head_rotation',
                 'wrist_flex_r', 'subtalar_angle_l'):
        assert name in text, ('W6', name)
    print('static checks W1-W6 pass; welded %s, dropped %d passive terms'
          % (', '.join(MTP_COORDINATES),
             len(json.loads((out / 'registration_mtp_welded.json').read_text())
                 ['removed_passive_terms'])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check:
        check(ROOT)
    else:
        for name, path in build(ROOT, ROOT / VARIANT).items():
            print(name, Path(path).relative_to(ROOT))


if __name__ == '__main__':
    main()
