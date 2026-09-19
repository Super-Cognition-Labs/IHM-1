#!/usr/bin/env python3
"""Two things a seat measurement cannot tell you about a skin contact bundle.

`scripts/verify_skin_contact.py` says where the skin sits relative to the bone
and the floor.  Both of the failures measured on 2026-09-18 are invisible to it:

**1. THE SUPPORT REGION AGAINST THE CENTRE OF MASS.**  A bundle can be seated
perfectly and still be unable to hold the body up, because what decides that is
WHERE the contact is relative to the weight line.  Measured from one snapshot,
using the plant's own emitted per-body mass and mass centre -- never the model
file's, which is not what the engine scaled.  The control is the sphere arm: at
the stance pose its force-weighted centre of pressure must equal the whole-body
centre of mass in x, because that pose IS a static equilibrium under the spheres.
If it does not, the reconstruction is wrong and nothing else here is readable.

**2. THE SEAM.**  A hard partition cuts a continuous surface, so a skin vertex on
a boundary belongs to the triangles of BOTH neighbouring pieces and is carried
twice.  The distance between its two images is the step the bundle presents to
the world.  Two controls, and the first of them failed the first time it was
written, which is why both are stated:

* at the pose the map was FITTED at, a single global map must carry both copies
  to the same place -- exactly zero.  It does (5.9e-14 mm).
* at any OTHER pose it is NOT zero and must not be: the seam opens by the joint's
  own motion away from the reference pose.  `radius_l|hand_l` reads exactly zero
  at every pose because the wrist is a `WeldJoint`, and that is what separates
  the joint's contribution from the map's.

The plant is the 22-segment SCAFFOLD.  This is its foot, not a human foot.
"""
from pathlib import Path
import argparse,gzip,json,sys,time
import numpy as np
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from ihm.assembly.anatomy_pose import OsimKinematics

MODEL='data/models/engineering_stance_v1/model.osim'
STANCE='data/models/engineering_stance_v1/initial_pose.json'
REGISTRATION='data/models/engineering_stance_v1/registration.json'
TARGET_MASS_KG=77.6122029
SEAMS=(('calcn_l','toes_l'),('calcn_r','toes_r'),('tibia_l','calcn_l'),('tibia_r','calcn_r'),
       ('femur_l','tibia_l'),('pelvis','torso'),('radius_l','hand_l'))
BUNDLES=('skin','skin-per-segment')


def read_obj(path):
    return np.asarray([l.split()[1:4] for l in Path(path).read_text().splitlines() if l.startswith('v ')],float)


def skin_partition():
    """the canonical exterior skin, and which piece owns each triangle"""
    mechanics=json.loads((ROOT/'data/derived/canonical/mechanics.json').read_text())
    skin=next(e for e in mechanics['entities'] if e['role']=='skin')
    geometry=json.loads(gzip.decompress((ROOT/skin['reference_geometry']['path']).read_bytes()))
    canonical=np.asarray(geometry['positions'],float).reshape(-1,3)
    faces=np.asarray(geometry['indices'],np.int64).reshape(-1,3)
    exterior=np.asarray(json.loads((ROOT/'data/research/engineered_skin_territories/materialization.json').read_text())['contact_eligible_triangle_ids'],np.int64)
    binding=json.loads(gzip.decompress((ROOT/'data/derived/canonical/continuous_surface_binding.json.gz').read_bytes()))
    segments=[s['id'] for s in binding['segments']]
    weights=np.asarray(binding['weights'],np.float32)
    owner=((weights[faces[:,0]]+weights[faces[:,1]]+weights[faces[:,2]])/3).argmax(axis=1)
    return canonical,faces,exterior,segments,owner


def support(root,out):
    from ihm.native.mechanical_stream import NativeMechanicalStream
    pose=json.loads((root/STANCE).read_text())
    work=root/'data/derived/skin-contact-support'/('run-'+str(int(time.time())))
    work.parent.mkdir(parents=True,exist_ok=True)
    stream=NativeMechanicalStream(root,work,environment='upright',target_mass_kg=TARGET_MASS_KG,
                                  initial_pose=pose,augmented_registration=REGISTRATION)
    try:first=stream.snapshot();second=stream.snapshot()
    finally:stream.close()
    # call it twice at the same input
    same=all(np.array_equal(np.asarray(first['bodies'][k]['transform_ground']),
                            np.asarray(second['bodies'][k]['transform_ground'])) for k in first['bodies'])
    transforms={k:np.asarray(v['transform_ground'],float) for k,v in first['bodies'].items()}
    mass=0.;centre=np.zeros(3)
    for name,body in first['bodies'].items():
        m=float(body['mass_kg']);local=np.asarray(body['mass_center_local_m'],float)
        centre+=m*(transforms[name][:3,:3]@local+transforms[name][:3,3]);mass+=m
    centre/=mass
    total=0.;moment=0.;xs=[]
    for contact in first['contacts']:
        force=contact['force_n'][1]
        if force<=1e-9:continue
        point=np.asarray(contact['center_m'],float);xs.append(float(point[0]))
        total+=force;moment+=force*point[0]
    report=dict(idempotent=bool(same),mass_kg=float(mass),com_ground_m=[float(v) for v in centre],
                sphere_arm=dict(loaded=len(xs),vertical_n=float(total),cop_x_m=float(moment/total),
                                centre_x_min_m=min(xs),centre_x_max_m=max(xs),
                                com_minus_cop_mm=float((centre[0]-moment/total)*1e3),
                                behind_com_mm=float((centre[0]-min(xs))*1e3),
                                ahead_of_com_mm=float((max(xs)-centre[0])*1e3)),
                bundles={})
    for bundle in BUNDLES:
        path=root/'data/derived/segment-contact-meshes'/bundle
        if not (path/'manifest.json').exists():continue
        manifest=json.loads((path/'manifest.json').read_text());below=[];bodies=[]
        for record in manifest['records']:
            vertices=read_obj(path/'meshes'/record['mesh_file'])
            T=transforms[record['body']];ground=vertices@T[:3,:3].T+T[:3,3]
            mask=ground[:,1]<0.
            if mask.any():below.append(ground[mask]);bodies.append(record['body'])
        if not below:report['bundles'][bundle]=dict(vertices_below_floor=0);continue
        points=np.concatenate(below)
        report['bundles'][bundle]=dict(vertices_below_floor=int(len(points)),bodies=bodies,
            x_min_m=float(points[:,0].min()),x_max_m=float(points[:,0].max()),
            behind_com_mm=float((centre[0]-points[:,0].min())*1e3),
            ahead_of_com_mm=float((points[:,0].max()-centre[0])*1e3))
    out['support']=report
    print('snapshot called twice, identical:',report['idempotent'])
    print('whole body %.7f kg, CoM ground x %+.6f m'%(report['mass_kg'],report['com_ground_m'][0]))
    s=report['sphere_arm']
    print('CONTROL  the sphere arm\'s centre of pressure minus the CoM in x: %+.4f mm'%s['com_minus_cop_mm'])
    print('         (the stance pose is a static equilibrium UNDER THE SPHERES, so this is 0)')
    print('%-18s %10s %10s %14s %14s'%('contact set','x min','x max','behind CoM','ahead of CoM'))
    print('%-18s %10.4f %10.4f %11.1f mm %11.1f mm'%('source spheres',s['centre_x_min_m'],s['centre_x_max_m'],
          s['behind_com_mm'],s['ahead_of_com_mm']))
    for bundle,row in report['bundles'].items():
        if not row.get('vertices_below_floor'):print('%-18s  nothing below the floor'%bundle);continue
        print('%-18s %10.4f %10.4f %11.1f mm %11.1f mm   (%d vertices on %s)'%(
            bundle,row['x_min_m'],row['x_max_m'],row['behind_com_mm'],row['ahead_of_com_mm'],
            row['vertices_below_floor'],','.join(row['bodies'])))


def seams(root,out):
    binding=json.loads((root/'data/derived/anatomy-segment-binding/binding.json').read_text())
    G=np.linalg.inv(np.asarray(binding['similarity_atlas_from_opensim_ground'],float))
    per=json.loads((root/'data/derived/anatomy-segment-registration/registration.json').read_text())
    kinematics=OsimKinematics(root/MODEL)
    node=ET.parse(root/MODEL).getroot().find('Model')
    defaults={c.get('name'):float(c.findtext('default_value') or 0.) for c in node.iter('Coordinate')}
    reference=kinematics.forward(kinematics.complete(dict(binding['reference_pose_rad'])))
    stance=kinematics.forward(kinematics.complete(json.loads((root/STANCE).read_text()),fill=defaults))
    canonical,faces,exterior,segments,owner=skin_partition()
    sets={seg:set(np.unique(faces[exterior[owner[exterior]==i]]).tolist()) for i,seg in enumerate(segments)}

    def place(segment,index,M,transforms):
        ground=canonical[index]@M[:3,:3].T+M[:3,3]
        inverse=np.linalg.inv(np.asarray(reference[segment],float))
        local=ground@inverse[:3,:3].T+inverse[:3,3]
        T=np.asarray(transforms[segment],float)
        return local@T[:3,:3].T+T[:3,3]

    rows={}
    print('\nthe same canonical skin vertex carried by BOTH neighbouring pieces (mm)')
    for label,transforms in (('reference',reference),('stance',stance)):
        print('  == at the %s pose'%label)
        print('  %-20s %7s %12s %12s %12s'%('seam','shared','global max','per-seg med','per-seg max'))
        for a,b in SEAMS:
            shared=sorted(sets[a]&sets[b])
            if not shared:continue
            index=np.asarray(shared)
            g=np.linalg.norm(place(a,index,G,transforms)-place(b,index,G,transforms),axis=1)*1e3
            Ma=np.asarray(per['segments'][a]['atlas_to_ground'],float)
            Mb=np.asarray(per['segments'][b]['atlas_to_ground'],float)
            p=np.linalg.norm(place(a,index,Ma,transforms)-place(b,index,Mb,transforms),axis=1)*1e3
            rows[label+' '+a+'|'+b]=dict(shared=len(shared),global_max_mm=float(g.max()),
                per_segment_median_mm=float(np.median(p)),per_segment_p90_mm=float(np.percentile(p,90)),
                per_segment_max_mm=float(p.max()))
            print('  %-20s %7d %12.4g %12.2f %12.2f'%(a+'|'+b,len(shared),g.max(),np.median(p),p.max()))
    out['seams']=rows
    welded=rows['stance radius_l|hand_l']['global_max_mm']
    print('  CONTROLS: global at the reference pose is 0 (one map, segments at the fit pose);')
    print('            at the stance pose it is NOT 0 and must not be -- the seam opens by the')
    print('            JOINT\'s own motion; radius_l|hand_l reads %.4g mm at BOTH poses because'%welded)
    print('            the wrist is a WeldJoint, which separates the two causes.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default=None)
    parser.add_argument('--skip-engine',action='store_true')
    args=parser.parse_args()
    report=dict(schema='ihm.skin-contact-support.v1',pose=STANCE,model=MODEL)
    if not args.skip_engine:support(ROOT,report)
    seams(ROOT,report)
    if args.out:
        out=ROOT/args.out;out.parent.mkdir(parents=True,exist_ok=True)
        out.write_text(json.dumps(report,indent=2)+'\n');print('\nwrote',out)
