"""Exact retained 80-muscle source catalog; cortical assignments are priors."""
import hashlib
from pathlib import Path
import xml.etree.ElementTree as ET
SOURCE_MODEL='data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled.osim'
# Schemas that describe a separately materialized native body. Both are real
# registered plants: `upperbody-registration.v1` is the base model plus the 12
# bilateral Arm26 actuators (92 muscles, the EmbodiedRuntime default), and
# `lumbar-muscle-variant.v1` is that body plus a 6-muscle Gait2392 trunk insert
# (98 muscles, what every stance controller forces via
# `data/models/engineering_stance_v1`). `ihm/assembly/embodied.py` has accepted
# both since the variant existed; this reader accepted only the first, so the
# fail-closed catalog and `peripheral_coverage` could not see the 98-muscle body
# the live stance path actually integrates. Two readers disagreeing about what a
# registered body is, is the disagreement itself -- keep this tuple and
# `_prepare_mechanical_registration` in step.
REGISTERED_EFFECTOR_SCHEMAS=('ihm.upperbody-registration.v1','ihm.lumbar-muscle-variant.v1')

def native_muscle_catalog(root):
    path=Path(root)/SOURCE_MODEL
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    rows=[]
    for muscle in ET.parse(path).getroot().iter('Millard2012EquilibriumMuscle'):
        name=muscle.attrib['name'];side=name[-1]
        if side not in ('r','l') or name[-2]!='_':raise ValueError('Unresolved source laterality')
        bodies=sorted({p.text.rsplit('/',1)[-1] for p in muscle.findall('.//PathPointSet/objects/*/socket_parent_frame')})
        if not bodies:raise ValueError('Missing source attachment bodies')
        distal=any(b.startswith(('calcn','talus','toes')) for b in bodies)
        knee=any(b.startswith(('tibia','patella')) for b in bodies)
        group='ankle_foot' if distal else 'knee' if knee else 'hip_pelvis'
        hemi='lh' if side=='r' else 'rh'
        rows.append({'id':name,'side':side,'body_group':group,'attachment_bodies':bodies,
            'sensory_region':f'brain-{hemi}-postcentral','motor_region':f'brain-{hemi}-precentral',
            'max_isometric_force_n':float(muscle.findtext('max_isometric_force')),
            'optimal_fiber_length_m':float(muscle.findtext('optimal_fiber_length')),
            'source_path':SOURCE_MODEL,'source_sha256':digest,
            'assignment_basis':'contralateral regional cortical assignment engineering prior; body group from source attachment names; no measured somatotopic recruitment'})
    if not rows or len({m['id'] for m in rows})!=len(rows):raise ValueError('Invalid source muscle catalog')
    return rows


def whole_body_effector_catalog(root,manifest):
    """Fail-closed catalog for the exact separately materialized native body.

    This does not mutate the default 80-muscle plant. Caller must instantiate
    the matching manifest model; XML construction alone is not native validation.
    """
    import json
    root=Path(root).resolve()
    if not isinstance(manifest,dict) or manifest.get('schema') not in REGISTERED_EFFECTOR_SCHEMAS:
        raise ValueError('Invalid registered effector manifest')
    def verified(relative,digest):
        if not isinstance(relative,str) or Path(relative).is_absolute():raise ValueError('Source path must be root-relative')
        path=(root/relative).resolve()
        if not path.is_relative_to(root) or hashlib.sha256(path.read_bytes()).hexdigest()!=digest:
            raise ValueError('Stale/outside registered effector source')
        return path
    model_path=verified(manifest['model_path'],manifest['model_sha256'])
    catalog_path=verified(manifest['catalog_path'],manifest['catalog_sha256'])
    for path,digest in manifest['sources'].items():verified(path,digest)
    if manifest['schema']=='ihm.lumbar-muscle-variant.v1':
        # A variant is a base model plus a ForceSet insert, so BOTH halves have to
        # be named and hash-verified here. Verifying only `model_path` would
        # certify a 98-muscle body while knowing where just the file came from and
        # not what was added to it. The base is required to be one of the already
        # verified `sources` rather than trusted from its own field, so the base
        # cannot be a path this function never hashed.
        insert_path=verified(manifest['insert_path'],manifest['insert_sha256'])
        base=manifest['base_model_path']
        if not isinstance(base,str) or manifest['sources'].get(base) is None:
            raise ValueError('Variant base model is not a verified registered source')
        if insert_path==model_path or (root/base).resolve()==model_path:
            raise ValueError('Variant model must differ from its base and insert')
    rows=json.loads(catalog_path.read_text())
    if not isinstance(rows,list) or len(rows)!=manifest['muscle_count'] or len({r['id'] for r in rows})!=len(rows):
        raise ValueError('Invalid materialized effector identities')
    model=ET.parse(model_path).getroot()
    effectors={m.get('name'):m for m in model.findall('.//ForceSet/objects/*') if 'Muscle' in m.tag}
    if set(effectors)!={row['id'] for row in rows}:raise ValueError('Catalog and actual plant effectors differ')
    import math
    for row in rows:
        if manifest['sources'].get(row['source_path'])!=row['source_sha256']:
            raise ValueError('Catalog source owner differs from verified manifest')
        muscle=effectors[row['id']]
        for field,xml_field in [('max_isometric_force_n','max_isometric_force'),('optimal_fiber_length_m','optimal_fiber_length')]:
            value=row[field]
            if isinstance(value,bool) or not math.isfinite(value) or value<=0 or not math.isclose(value,float(muscle.findtext(xml_field)),rel_tol=1e-12):
                raise ValueError('Catalog normalization differs from actual plant')
    return rows
