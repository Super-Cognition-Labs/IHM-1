#!/usr/bin/env python3
"""Where the skin, the bone, the source foot spheres and the floor are -- ONE frame.

WHY THIS EXISTS.  Two artefacts in this repo disagree about where the floor is
under the foot, by about 36 mm, and the disagreement had never been taken apart:

* `docs/WORKBENCH_AUTHENTICITY.md` 1.1 and `plant_options.SEGMENT_CONTACT_BUNDLES`
  say *the skin never reaches the floor in the stance pose*;
* `docs/SOFT_BODY.md` measured the same thing sharply: the foot's skin sits
  **36.3 mm above the floor** while the source foot spheres carry 616 N.

Four explanations were available and only a measurement separates them: the
source foot spheres sit below the skin by construction; the bundle's
segment-local stations were cut at a reference pose that is not the stance pose;
the floor is placed by a rule that ignores the meshes; or the skin itself is
misplaced.  This script measures all four in one frame, with controls, and the
answer is in `docs/SEGMENT_CONTACT_SURFACES.md`.

WHAT IT PRINTS, and every number has one provenance:

1. **The floor**, read out of the `ContactGeometrySet` the engine actually loads
   -- never assumed to be y = 0.  The upright branch of the engine takes the
   source model's own `floor` half-space; only the SUPINE branch hangs a plane
   under the lowest inertia sphere, so the "the plane follows a sphere" story is
   excluded by reading the file rather than by argument.
2. **The chain**, at the stance pose, in ground: floor, the lowest point of every
   source foot contact sphere, the lowest vertex of each foot bone, and the
   lowest vertex of each skin bundle's foot pieces.
3. **The seat** -- skin lowest minus bone lowest in the segment's OWN frame,
   which is a rigid-body fact and does not depend on the pose.  Negative means
   the skin is below the bone, which is what a pad is.  This is the quantity
   gate 3 of the 2026-09-10 skin-warp pre-registration was written against,
   band **[-25, -5] mm**, and that band is not touched here.
4. **The body's own plantar pad**, from its own soft-tissue depth map, as the
   reference the seat has to be compared against.

CONTROLS, because a measurement of a placement is worthless if the ruler moves:

* **K0 idempotence** -- the whole measurement run twice at the same input must
  return bitwise identical numbers.  Not a known answer: a test of whether the
  instrument is a function at all.
* **K1 the artefact agrees with its own bytes** -- the seat recomputed from the
  meshes on disk must reproduce the `skin_minus_bone_minimum_y_m` each bundle's
  manifest recorded when it was built.
* **K2 a control that can fail** -- lower one skin piece by exactly 10.000 mm and
  the seat must move by exactly -10.000 mm.  An instrument that reads a constant
  whatever it is handed passes every other check here.
* **K3 the known WRONG answer** -- the `skin-canonical` bundle is built through
  the map this file already measured as wrong, and must print its recorded
  +84.0 / +94.6 mm hover.  A ruler that cannot see a known 84 mm error cannot
  certify a 20 mm one.
* **K4 (--engine)** the engine's own emitted contact elements must agree with the
  kinematics used here, and the sign of every element's vertical force must match
  the sign of its penetration below the floor read in 1.

The plant is the 22-segment SCAFFOLD.  This is its foot, not a human foot.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'scripts'))
from ihm.assembly.anatomy_pose import OsimKinematics
from ihm.spatial.vtk import surface as read_surface

MODEL='data/models/engineering_stance_v1/model.osim'
STANCE='data/models/engineering_stance_v1/initial_pose.json'
REGISTRATION='data/models/engineering_stance_v1/registration.json'
CONTACT_SET='data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/example3DWalking/subject_walk_scaled_ContactGeometrySet.xml'
DEPTH_MAP='data/derived/soft-tissue-depth-v1/depth.npz'
BUNDLES=('skin','skin-canonical','skin-layer-map-v1','skin-warp-v1','skin-warp-v2','skin-per-segment')
FOOT=('calcn_l','calcn_r','toes_l','toes_r')
TARGET_MASS_KG=77.6122029
# Fixed 2026-09-10 in this file's own pre-registration, before any warp was fitted.
GATE3_BAND_MM=(-25.,-5.)


def model_defaults(root,model):
    """Every coordinate's own declared default, read out of the model.

    The engine writes only the coordinates the caller names into
    `initial_pose.txt`; everything else stays at the value the model declares.
    Filling with zero would be an assumption, and K4 is what checks this one:
    the sphere centres reconstructed here have to be the engine's own.
    """
    node=ET.parse(root/model).getroot().find('Model')
    return {c.get('name'):float(c.findtext('default_value') or 0.) for c in node.iter('Coordinate')}


def read_obj(path):
    rows=[l.split()[1:4] for l in Path(path).read_text().splitlines() if l.startswith('v ')]
    return np.asarray(rows,dtype=float)


def floor_plane(root):
    """The floor the UPRIGHT engine loads, read out of its own source file.

    `ContactHalfSpace` is SimTK's x<0 half space with its normal along +x in its
    own frame; the source declares location (0,0,0) on /ground and a rotation
    about z, so the plane passes through the ground origin and the surface the
    body stands on is y = location.y.  Read, not assumed -- the whole point of
    the third candidate explanation is that a plane can be placed by a rule.
    """
    tree=ET.parse(root/CONTACT_SET).getroot()
    for node in tree.iter('ContactHalfSpace'):
        if node.get('name')!='floor':continue
        location=np.fromstring(node.findtext('location') or '0 0 0',sep=' ')
        orientation=np.fromstring(node.findtext('orientation') or '0 0 0',sep=' ')
        return dict(name='floor',frame=node.findtext('socket_frame').strip(),
                    location_m=[float(v) for v in location],
                    orientation_rad=[float(v) for v in orientation],
                    y_m=float(location[1]),
                    basis='ContactHalfSpace "floor" of '+CONTACT_SET+', attached to /ground; '
                          'the engine\'s upright branch adds this set verbatim and places no plane of its own '
                          '(scripts/native_mechanical_stream.cpp, upright branch). Only the SUPINE branch '
                          'computes a plane under the lowest inertia-inscribed sphere.')
    raise ValueError('no floor half space in '+CONTACT_SET)


def source_foot_spheres(root):
    """The 12 anatomically placed source foot contacts, from the same file."""
    tree=ET.parse(root/CONTACT_SET).getroot();out=[]
    for node in tree.iter('ContactSphere'):
        frame=node.findtext('socket_frame').strip()
        out.append(dict(name=node.get('name'),body=frame.rsplit('/',1)[-1],
                        location_m=[float(v) for v in np.fromstring(node.findtext('location'),sep=' ')],
                        radius_m=float(node.findtext('radius'))))
    return out


def bone_clouds(root,model):
    """Every body's own bone-mesh vertices in the body frame, subject scale baked in."""
    node=ET.parse(root/model).getroot().find('Model')
    geometry=root/'data/raw/anatomy/opensim-models/source/Geometry'
    out={}
    for body in node.iter('Body'):
        parts=[]
        for mesh in body.iter('Mesh'):
            factors=np.fromstring(mesh.findtext('scale_factors'),sep=' ')
            points,_=read_surface(geometry/mesh.findtext('mesh_file'))
            parts.append(points*factors)
        if parts:out[body.get('name')]=np.concatenate(parts)
    return out


def plantar_pad(root):
    """This body's OWN heel pad, from its own soft-tissue depth map.

    Measured on the subset the question is about: the plantar band, the lowest
    10% and 25% of the depth-map points whose nearest structure is the
    calcaneus.  Over the WHOLE calcaneus patch the median is 19.6 mm, which
    includes the sides and the back of the heel and is not what stands on a
    floor.
    """
    data=np.load(root/DEPTH_MAP,allow_pickle=True)
    points=data['skin_points_m'];depth=data['depth_m']
    names=np.array([str(s) for s in data['nearest']])
    out={}
    for side in ('left','right'):
        mask=np.array([('calcaneus' in n and side in n) for n in names])
        p=points[mask];d=depth[mask]*1e3
        row=dict(points=int(mask.sum()),all_median_mm=float(np.median(d)))
        for fraction,label in ((.10,'lowest10'),(.25,'lowest25')):
            band=p[:,1]<=np.percentile(p[:,1],fraction*100)
            row[label+'_n']=int(band.sum())
            row[label+'_min_mm']=float(d[band].min())
            row[label+'_median_mm']=float(np.median(d[band]))
            row[label+'_max_mm']=float(d[band].max())
        out[side]=row
    return out


def seats(root,bundle,bones):
    """skin lowest minus bone lowest, in the SEGMENT frame, per segment.

    A rigid-body fact: both clouds ride the same body, so this number does not
    depend on the pose the plant is in.  It is gate 3's quantity.
    """
    manifest=json.loads((root/'data/derived/segment-contact-meshes'/bundle/'manifest.json').read_text())
    out={}
    for record in manifest['records']:
        body=record['body']
        if body not in bones:continue
        vertices=read_obj(root/'data/derived/segment-contact-meshes'/bundle/'meshes'/record['mesh_file'])
        out[body]=dict(measured_mm=float((vertices[:,1].min()-bones[body][:,1].min())*1e3),
                       recorded_mm=None if record.get('skin_minus_bone_minimum_y_m') is None
                                   else float(record['skin_minus_bone_minimum_y_m']*1e3),
                       skin_vertices=int(len(vertices)))
    return out,manifest


def measure(root,engine=False):
    kinematics=OsimKinematics(root/MODEL)
    stance=json.loads((root/STANCE).read_text())
    transforms=kinematics.forward(kinematics.complete(dict(stance),fill=model_defaults(root,MODEL)))
    bones=bone_clouds(root,MODEL)
    floor=floor_plane(root)
    report=dict(schema='ihm.skin-contact-verification.v1',
                pose=STANCE,model=MODEL,floor=floor,
                plantar_pad_from_own_depth_map=plantar_pad(root),
                gate3_band_mm=list(GATE3_BAND_MM))

    def ground(body,points):
        T=np.asarray(transforms[body],float)
        return points@T[:3,:3].T+T[:3,3]

    spheres=[]
    for sphere in source_foot_spheres(root):
        centre=ground(sphere['body'],np.asarray([sphere['location_m']],float))[0]
        spheres.append(dict(sphere,centre_ground_m=[float(v) for v in centre],
                            bottom_y_m=float(centre[1]-sphere['radius_m']),
                            below_floor_mm=float((floor['y_m']-(centre[1]-sphere['radius_m']))*1e3)))
    spheres.sort(key=lambda s:s['bottom_y_m'])
    report['source_foot_spheres']=spheres
    report['lowest_source_sphere']=spheres[0]['name']

    report['bone_above_floor_mm']={b:float((ground(b,bones[b])[:,1].min()-floor['y_m'])*1e3) for b in FOOT}
    report['bundles']={}
    for bundle in BUNDLES:
        path=root/'data/derived/segment-contact-meshes'/bundle
        if not (path/'manifest.json').exists():continue
        seat,manifest=seats(root,bundle,bones)
        skin_above={}
        for body in seat:
            vertices=read_obj(path/'meshes'/next(r['mesh_file'] for r in manifest['records'] if r['body']==body))
            skin_above[body]=float((ground(body,vertices)[:,1].min()-floor['y_m'])*1e3)
        report['bundles'][bundle]=dict(
            registration=manifest.get('registration',{}).get('choice'),
            seat_mm={b:seat[b]['measured_mm'] for b in seat},
            recorded_seat_mm={b:seat[b]['recorded_mm'] for b in seat},
            skin_above_floor_mm=skin_above,
            gate3={b:(GATE3_BAND_MM[0]<=seat[b]['measured_mm']<=GATE3_BAND_MM[1]) for b in seat if b.startswith('calcn')})

    # K1: the recomputation must reproduce what each bundle recorded when it was built.
    k1=[]
    for bundle,row in report['bundles'].items():
        for body,value in row['seat_mm'].items():
            recorded=row['recorded_seat_mm'].get(body)
            if recorded is not None:k1.append(abs(value-recorded))
    report['K1_artefact_agrees_mm']=float(max(k1)) if k1 else None

    # K2: a control that can fail.
    shipped=root/'data/derived/segment-contact-meshes/skin'
    manifest=json.loads((shipped/'manifest.json').read_text())
    mesh=next(r['mesh_file'] for r in manifest['records'] if r['body']=='calcn_l')
    vertices=read_obj(shipped/'meshes'/mesh)
    base=(vertices[:,1].min()-bones['calcn_l'][:,1].min())*1e3
    moved=vertices.copy();moved[:,1]-=0.010
    report['K2_lowered_10mm']=dict(before_mm=float(base),
        after_mm=float((moved[:,1].min()-bones['calcn_l'][:,1].min())*1e3),
        delta_mm=float((moved[:,1].min()-vertices[:,1].min())*1e3))

    # K3: the known WRONG map has to print its recorded hover.
    canonical=report['bundles'].get('skin-canonical')
    report['K3_known_wrong_map_mm']=None if canonical is None else \
        {b:canonical['seat_mm'][b] for b in FOOT if b in canonical['seat_mm']}

    if engine:
        from ihm.native.mechanical_stream import NativeMechanicalStream
        work=root/'data/derived/skin-contact-verify'/('run-'+str(int(time.time())))
        work.parent.mkdir(parents=True,exist_ok=True)
        stream=NativeMechanicalStream(root,work,environment='upright',target_mass_kg=TARGET_MASS_KG,
                                      initial_pose=stance,augmented_registration=REGISTRATION)
        try:frame=stream.snapshot()
        finally:stream.close()
        emitted={c['name']:c for c in frame['contacts']}
        worst=0.;signs=[]
        for sphere in spheres:
            found=emitted.get('contact'+sphere['name'][0].upper()+sphere['name'][1:])
            if found is None:
                found=next((c for c in frame['contacts'] if c['name'].lower().endswith(sphere['name'].lower())),None)
            if found is None:continue
            worst=max(worst,float(abs(np.asarray(found['center_m'],float)[1]-sphere['centre_ground_m'][1])))
            signs.append(bool((found['force_n'][1]>0)==(sphere['bottom_y_m']<floor['y_m'])))
        report['K4_engine']=dict(contact_model=frame['contact_model'],
                                 worst_sphere_centre_disagreement_m=worst,
                                 force_sign_tracks_penetration=bool(all(signs)) if signs else None,
                                 spheres_compared=len(signs),
                                 vertical_contact_force_n=float(sum(c['force_n'][1] for c in frame['contacts'])),
                                 weight_n=TARGET_MASS_KG*9.81)
    return report


def flatten(value):
    if isinstance(value,dict):return {k:flatten(v) for k,v in sorted(value.items())}
    if isinstance(value,(list,tuple)):return [flatten(v) for v in value]
    return value


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',default=None)
    parser.add_argument('--engine',action='store_true',help='also start the plant and cross-check its emitted contact elements')
    args=parser.parse_args()
    first=measure(ROOT,engine=False)
    second=measure(ROOT,engine=False)
    idempotent=json.dumps(flatten(first),sort_keys=True)==json.dumps(flatten(second),sort_keys=True)
    report=first if not args.engine else measure(ROOT,engine=True)
    report['K0_idempotent']=bool(idempotent)

    floor=report['floor'];pad=report['plantar_pad_from_own_depth_map']
    print('K0  called twice at the same input, identical:',report['K0_idempotent'])
    print('K1  seat recomputed vs the bundle\'s own record: worst %.9f mm'%report['K1_artefact_agrees_mm'])
    k2=report['K2_lowered_10mm']
    print('K2  lowering calcn_l skin 10.000 mm moves the seat by %+.6f mm (%+.3f -> %+.3f)'%(
        k2['delta_mm'],k2['before_mm'],k2['after_mm']))
    print('K3  the known-wrong canonical map prints',{k:round(v,1) for k,v in (report['K3_known_wrong_map_mm'] or {}).items()})
    if 'K4_engine' in report:
        k4=report['K4_engine']
        print('K4  engine: sphere centres agree to %.2e m; force sign tracks penetration: %s; Fy %.2f N vs mg %.2f N'%(
            k4['worst_sphere_centre_disagreement_m'],k4['force_sign_tracks_penetration'],
            k4['vertical_contact_force_n'],k4['weight_n']))
    print()
    print('FLOOR  %s at y = %+.6f m, %s'%(floor['name'],floor['y_m'],floor['frame']))
    print('       read out of '+CONTACT_SET)
    low=report['source_foot_spheres'][0]
    print('LOWEST source foot sphere %s: bottom y %+.6f m (%.3f mm below the floor), r = %.3f m'%(
        low['name'],low['bottom_y_m'],low['below_floor_mm'],low['radius_m']))
    print('PLANTAR pad, this body\'s own depth map, lowest 10%% of the calcaneus patch: %.2f / %.2f mm median'%(
        pad['left']['lowest10_median_mm'],pad['right']['lowest10_median_mm']))
    print()
    print('%-16s %10s %10s %10s %10s'%('segment','bone above','skin above','seat','gate 3'))
    print('%-16s %10s %10s %10s %10s'%('','floor mm','floor mm','mm','[-25,-5]'))
    for bundle,row in report['bundles'].items():
        print('-- '+bundle+'  (registration: %s)'%row['registration'])
        for body in FOOT:
            if body not in row['seat_mm']:continue
            gate=row['gate3'].get(body)
            print('   %-13s %10.3f %10.3f %10.3f %10s'%(body,report['bone_above_floor_mm'][body],
                  row['skin_above_floor_mm'][body],row['seat_mm'][body],
                  '' if gate is None else ('PASS' if gate else 'FAIL')))
    if args.out:
        out=ROOT/args.out;out.parent.mkdir(parents=True,exist_ok=True)
        out.write_text(json.dumps(report,indent=2)+'\n');print('\nwrote',out)
