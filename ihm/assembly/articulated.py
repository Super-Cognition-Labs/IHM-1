"""Native articulated inertia with explicit canonical material registration.

Every native segment owns mass/inertia once. Canonical tissue meshes are attached
geometric materializations; this wrapper does not create 2408 independent masses
or pretend that attachment geometry is an all-organ volumetric contact solve.
"""
from pathlib import Path
import copy,hashlib,json
import numpy as np
from ihm.native.mechanical_stream import NativeMechanicalStream,finite,vec
from .sensorimotor_catalog import native_muscle_catalog
from .snapshot_data import clone_snapshot_data
from .continuous_surface_binding import ContinuousSurfaceBinding

BASIS=np.array([[0.,0.,-1.],[0.,1.,0.],[1.,0.,0.]])

class CanonicalRegistration:
    def __init__(self,payload,native_reference):
        # A requested initial pose changes geometry, not the material embedding.
        # Retain the source-default body transforms emitted before posing.
        if 'registration_reference_bodies' in native_reference:
            native_reference=copy.deepcopy(native_reference)
            native_reference['bodies']=native_reference['registration_reference_bodies']
        self.specs={e['id']:e for e in payload['entities']};self.ids=list(self.specs);self.reference=copy.deepcopy(native_reference)
        self.bodies=list(native_reference['bodies']);self.groups={};self.named={};self.maps={};self.rows={};source_landmarks=[];canonical_landmarks=[]
        for body in self.bodies:
            entry=payload.get('registration',{}).get(body,{})
            bones=[self.specs[k] for k in entry.get('canonical_bones',[]) if k in self.specs]
            if not bones:raise ValueError('Native segment has no retained canonical bone anchors: '+body)
            lower=np.min([e['bounds_m']['min'] for e in bones],axis=0);upper=np.max([e['bounds_m']['max'] for e in bones],axis=0);anchor=(lower+upper)/2
            for e in bones:
                if e['id'] in self.named and self.named[e['id']]!=body:raise ValueError('Canonical bone has multiple segment owners')
                self.named[e['id']]=body
            record=native_reference['bodies'][body];t=np.array(record['transform_ground']);com=t[:3,:3]@record['mass_center_local_m']+t[:3,3]
            source_landmarks.append(com);canonical_landmarks.append(anchor)
            alignment=np.eye(4)
            self.maps[body]=alignment;self.groups[body]={'bounds_min_m':lower,'bounds_max_m':upper,'anchor_m':anchor,'canonical_bones':[e['id'] for e in bones],
                'source_to_canonical_ground':alignment,'native_reference_transform':t,
                'basis':'Named retained bone group provides an approximate canonical anchor for one shared global rigid fit; COM and bounding-box center are not measured homologous joint landmarks.'}
        a=np.array(source_landmarks);b=np.array(canonical_landmarks);ac=a.mean(0);bc=b.mean(0)
        u,_,vh=np.linalg.svd((a-ac).T@(b-bc));r=vh.T@np.diag([1.,1.,np.linalg.det(vh.T@u.T)])@u.T
        self.basis=r;self.global_map=np.eye(4);self.global_map[:3,:3]=r;self.global_map[:3,3]=bc-r@ac
        residual=(a@r.T+self.global_map[:3,3])-b
        self.global_fit={'basis':'Unweighted proper-rigid least-squares fit of 22 approximate source COM/canonical bone-envelope center correspondences; no scale fit',
                         'rms_landmark_residual_m':float(np.sqrt(np.mean(np.sum(residual**2,axis=1)))),'maximum_landmark_residual_m':float(np.linalg.norm(residual,axis=1).max()),
                         'source_to_canonical_ground':self.global_map.tolist(),'landmark_residuals_m':dict(zip(self.bodies,residual.tolist()))}
        for body in self.bodies:
            self.maps[body]=self.global_map.copy();self.groups[body]['source_to_canonical_ground']=self.global_map.copy()
            self.groups[body]['canonical_reference_embedding_in_source_body']=np.linalg.inv(self.groups[body]['native_reference_transform'])@np.linalg.inv(self.global_map)
        for ident,e in self.specs.items():
            point=np.array(e['centroid_m']);rank=self._ranking(point)
            owner=self.named.get(ident,rank[0][1]);self.rows[ident]={'body':owner,'basis':'retained_named_bone_group' if ident in self.named else 'nearest_named_bone_envelope_prior',
                    'candidate_distances_m':[[b,d] for d,b in rank[:2]],'registration_error_m':None,
                    'unresolved':['inter-segment anatomical joint-gap/soft-continuity error not calibrated'] if ident in self.named else ['tissue compartment/segment assignment is inferred','independent internal tissue deformation unresolved']}
    def _ranking(self,point):
        return sorted((float(np.linalg.norm(np.maximum(np.maximum(g['bounds_min_m']-point,point-g['bounds_max_m']),0))),b) for b,g in self.groups.items())
    def transforms(self,native):
        out={}
        for body in self.bodies:
            c=self.maps[body];t=np.array(native['bodies'][body]['transform_ground']);reference=self.groups[body]['native_reference_transform']
            out[body]=c@t@np.linalg.inv(reference)@np.linalg.inv(c)
        return out
    def project(self,native):
        transforms=self.transforms(native);entities={}
        for ident,e in self.specs.items():
            t=transforms[self.rows[ident]['body']];center=np.array(e['centroid_m']);current=t[:3,:3]@center+t[:3,3]
            entities[ident]={'centroid_m':current.tolist(),'translation_m':(current-center).tolist(),'rotation_matrix':t[:3,:3].tolist(),'deformation_gradient':np.eye(3).tolist()}
        return entities
    def force(self,ident,point,force,native):
        if ident not in self.rows:raise ValueError('Unknown canonical force entity')
        point=vec(point);force=vec(force);owner=self.rows[ident]['body']
        # Extensive source surfaces require spatial support, not the segment
        # nearest their single whole-surface centroid. Preserve the chosen prior.
        if ident not in self.named:
            candidates=[]
            for b,t in self.transforms(native).items():
                reference=t[:3,:3].T@(point-t[:3,3]);g=self.groups[b]
                distance=float(np.linalg.norm(np.maximum(np.maximum(g['bounds_min_m']-reference,reference-g['bounds_max_m']),0)))
                candidates.append((distance,b))
            owner=min(candidates)[1]
        c=self.maps[owner];source_point=self.basis.T@(point-c[:3,3]);source_force=self.basis.T@force
        return {'body':owner,'point_m':source_point.tolist(),'force_n':source_force.tolist()}
    def contact_wrenches(self,native):
        records=[]
        for contact in native['contacts']:
            body=contact['body_frame']
            if body not in self.maps:raise ValueError('Native contact frame has no canonical registration')
            c=self.maps[body];origin=np.array(native['bodies'][body]['transform_ground'])[:3,3]
            records.append({'body':body,'canonical_bone_ids':self.groups[body]['canonical_bones'],
                'point_m':(self.basis@origin+c[:3,3]).tolist(),'force_n':(self.basis@contact['force_n']).tolist(),
                'moment_nm':(self.basis@contact['moment_nm']).tolist(),'name':contact['name'],
                'kind':'resultant_wrench_about_body_origin','frame':'canonical_current_world',
                'scope':'Rigid-segment contact resultant; distribution onto deforming thorax/skin is unresolved; proxy area is not anatomical pressure area'})
        return records
    def cutaneous_contacts(self,native):
        surface=native.get('surface_foundation')
        if surface is None:return None
        result=[]
        for point in surface['sensor_points']:
            row=copy.deepcopy(point)
            row['point_m']=(self.basis@np.array(point['point_source_m'])+self.global_map[:3,3]).tolist()
            row['normal']=(self.basis@np.array(point['normal_source'])).tolist()
            row['force_n']=(self.basis@np.array(point['force_n'])).tolist()
            row['coordinate_frame']='canonical_current_world'
            result.append(row)
        return result

    def manifest(self):
        def plain(value):
            if isinstance(value,np.ndarray):return value.tolist()
            if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
            if isinstance(value,list):return [plain(v) for v in value]
            return value
        return plain({'schema':'ihm.canonical-articulation.v1','groups':self.groups,'entities':self.rows,
            'coordinate_basis_source_to_canonical':self.basis,'global_rigid_fit':self.global_fit,
            'scope':'One global ground frame preserves native joint and environment consistency. Canonical source meshes have fixed segment-local embeddings; their anatomical boundaries are not guaranteed to coincide with native joint locations. No anatomical registration precision is fabricated.'})

class ArticulatedBodyPlant:
    # `mechanical_fidelity` carries the joint stops, the real segment contact
    # surfaces and the derived tissue force elements that NativeMechanicalStream has
    # always accepted and that nothing serving a live body could previously ask for
    # (docs/WORKBENCH_AUTHENTICITY.md, Tier 1).  It is resolved from server-owned
    # identities by ihm/assembly/plant_options.py; a caller never supplies a path.
    # None is the historical plant to the float, so no existing measurement moves.
    def __init__(self,root,output,*,environment='supine',target_mass_kg=None,augmented_registration=None,enable_garments=False,surface_contact_manifest=None,surface_sensor_indices=(),bed_material=None,instance_mass_variant=None,initial_pose=None,mechanical_fidelity=None,display_pose=None):
        self.root=Path(root).resolve();self.output=Path(output).resolve()
        if self.output.exists() or not self.output.is_relative_to(self.root):raise ValueError('Fresh owned articulated output required')
        path=self.root/'data/derived/canonical/mechanics.json';raw=path.read_bytes();payload=json.loads(raw)
        canonical_mass=sum(e['mass_kg'] for e in payload['entities']);target=canonical_mass if target_mass_kg is None else finite(target_mass_kg)
        self.output.mkdir(parents=True);(self.output/'canonical_mechanics.json').write_bytes(raw)
        from .plant_options import resolve_fidelity
        fidelity_kwargs,self.mechanical_fidelity=resolve_fidelity(self.root,mechanical_fidelity,environment=environment)
        self.native=NativeMechanicalStream(self.root,self.output/'native',environment=environment,target_mass_kg=target,augmented_registration=augmented_registration,surface_contact_manifest=surface_contact_manifest,surface_sensor_indices=surface_sensor_indices,bed_material=bed_material,instance_mass_variant=instance_mass_variant,**fidelity_kwargs,**({'initial_pose':initial_pose} if initial_pose is not None else {}))
        try:
            self.registration=CanonicalRegistration(payload,self.native.snapshot());self.muscle_catalog=self.native.muscle_catalog or native_muscle_catalog(self.root)
            self.surface_binding=ContinuousSurfaceBinding.from_root(self.root,self.registration)
            (self.output/'surface_binding.json').write_text(json.dumps(self.surface_binding.manifest(),indent=2)+'\n')
            self.source_registration_manifest=self.registration.manifest();(self.output/'registration.json').write_text(json.dumps(self.source_registration_manifest,indent=2)+'\n')
            self.identity={'canonical_mechanics_sha256':hashlib.sha256(raw).hexdigest(),'canonical_reference_mass_kg':canonical_mass,'effective_native_body_mass_kg':target,
                          'mass_change_basis':'Explicit uniform scaling of source segment masses and inertias; source proportions retained; no second canonical inertial owner'}
            (self.output/'identity.json').write_text(json.dumps(self.identity,indent=2)+'\n')
            (self.output/'mechanical_fidelity.json').write_text(json.dumps(self.mechanical_fidelity,indent=2)+'\n')
            # THE DISPLAY POSE, and why it is not `entities`.
            # `registration.project()` is a MATERIAL EMBEDDING: its transforms are
            # relative to the plant's own reference bodies, so force routing, the
            # surface binding and the world exchange can all apply loads in a frame
            # that does not move when a caller asks for an initial pose. Measured, at
            # the plant's own first frame with pelvis_ty=0.93: every entity's
            # `translation_m` is EXACTLY 0.00 mm while the mechanical body sits 88.5 mm
            # below the pose the atlas was registered at. That is correct for a force
            # frame and wrong for a picture -- the anatomy is drawn where the atlas left
            # it while the scaffold walks away underneath.
            # AnatomyPoser is the pose (FK verified against Simbody to 7.8e-16,
            # scripts/verify_anatomy_pose.py). It costs 0.94 ms median a frame. It is
            # carried BESIDE `entities`, never instead of it, because swapping it in
            # would move the frame that forces are applied in.
            self.display_poser=None;self.display_pose_basis=None;self.display_pose_map=None
            if display_pose is not None:
                if display_pose not in ('opensim','anatomical'):raise ValueError('Unknown display pivot')
                from .anatomy_pose import AnatomyPoser
                self.display_poser=AnatomyPoser.from_workspace(self.root,pivot=display_pose)
                # static for the session: which segment carries each entity, and each
                # entity's rest centroid. A client applies segment_motion[segment_of[i]]
                # to entity i and gets the per-entity pose back exactly.
                self.display_pose_map={'schema':'ihm.display-pose-map.v1',
                    'frame':'bodyparts3d-display-m','pivot':display_pose,
                    'segments':list(self.display_poser.segments),
                    'entity_ids':list(self.display_poser.entity_ids),
                    'segment_of':[int(v) for v in self.display_poser.segment_of],
                    'rest_centroid_m':self.display_poser.rest_centroid.tolist(),
                    'apply':'x_now = R[segment_of[i]] @ x_rest + t[segment_of[i]], with R and t '
                            'the 3x3 and 3x1 blocks of segment_motion[segment_of[i]]'}
                (self.output/'display_pose_map.json').write_text(json.dumps(self.display_pose_map)+'\n')
                self.display_pose_basis={'registration':'data/derived/anatomy-segment-binding/binding.json',
                    'pivot':display_pose,'frame':'bodyparts3d-display-m',
                    'verification':'scripts/verify_anatomy_pose.py -- forward kinematics to 7.8e-16 of '
                                   "Simbody's own transform_ground over 12 native frames, rest pose to "
                                   '5e-16, distal-only motion with 0 violations, bit-identical idempotence',
                    'basis':'The pose of the anatomy, for display. NOT the frame forces are applied in: '
                            "`entities` remains the plant's material embedding and is what every force, "
                            'the surface binding and the world exchange use.',
                    'pivot_basis':('OpenSim joint centres carried into the atlas frame; the anatomy follows '
                                   'the simulated segment exactly and the joint it turns about misses the '
                                   "atlas's own by 38 mm median" if display_pose=='opensim' else
                                   'Each joint re-seated where the two segments\' bone surfaces are closest: '
                                   'orientations stay exact and joints open less (knee 32 mm against 108 mm), '
                                   'at the price of entities sitting up to 50 mm off the simulated segment'),
                    'not_posed':'Entities with no rigid pose are listed on every frame; nothing is silently '
                                'left at rest. One rigid segment per entity: no soft-tissue deformation, no '
                                'volume preservation, no sliding, and entities straddling a joint tear at it.'}
            self.garments=None
            if enable_garments:
                from .garment_feedback import GarmentFeedback
                self.garments=GarmentFeedback.from_root(self.root,self.registration)
                (self.output/'garment_identity.json').write_text(json.dumps(self.garments.identity,indent=2)+'\n')
            self.state=self._project(self.native.snapshot(),0.)
        except BaseException:self.native.close();raise
    def _display_pose(self,native):
        """The anatomy's pose this frame, as 22 segment motions rather than 3,995 poses.

        A rigid binding gives every entity on a segment the SAME motion, so the
        per-entity form was 3,995 copies of 22 matrices -- 1.16 MB of JSON a frame,
        which is not a per-frame payload at interactive rates. This sends the 22, and
        the static entity->segment map goes once as `display_pose_map` beside the
        geometry. Reconstruction is exact, not approximate, and
        scripts/verify_display_pose.py checks it against the per-entity pose.
        """
        M=self.display_poser.segment_motions({b:r['transform_ground'] for b,r in native['bodies'].items()})
        return {**self.display_pose_basis,'segments':list(self.display_poser.segments),
                'segment_motion':[m.tolist() for m in M],
                'entity_segment_map':'display_pose_map.json',
                'unbound':copy.deepcopy(self.display_poser.unbound)}

    def _contact_scope(self):
        """What the body is actually touching the world with, on THIS plant.

        The sentence used to be a constant asserting inertia-derived spheres. With a
        segment contact bundle loaded that assertion is false, and a scope field that
        keeps saying `engineering proxies, not anatomical surface` about a run that
        stood on its skin is exactly the stale declaration this repo keeps catching.
        """
        contact=self.mechanical_fidelity.get('segment_contact')
        if not contact:
            return ('Unilateral native contact; source feet and disclosed environment-specific '
                    'inertia-derived spheres are engineering proxies, not anatomical surface or '
                    'calibrated contact')
        if contact['layer']=='skin':
            return ('Unilateral native contact over the real per-segment SKIN exterior as '
                    'ContactMesh under an elastic foundation. The skin is carried rigidly by its '
                    'segment: no in-plane stretch, no sliding, no deformable continuum, and every '
                    'segment boundary is a seam the real body does not have. Not calibrated contact.')
        return ('Unilateral native contact over real per-segment BONE surfaces. Bone is a collider '
                'against other bones and its own soft tissue; it is NOT what meets the floor in a '
                'real body. Not anatomical skin contact and not calibrated contact.')

    def _project(self,native,positive_work):
        frame={'model_id':'ihm-body','effective_native_body_mass_kg':native['mass_kg'],'mass_transfer':copy.deepcopy(native['mass_transfer']),'time_s':native['time_s'],'entities':self.registration.project(native),'muscles':clone_snapshot_data(native['muscles']),
                'foot_contact_force_n':copy.deepcopy(native['foot_contact_force_n']),'positive_muscle_work_j':float(positive_work),
                'total_muscle_metabolic_w':native.get('total_muscle_metabolic_w'),'muscle_metabolic_energy_j':native.get('muscle_metabolic_energy_j'),
                'signed_active_fiber_work_j':native['signed_active_fiber_work_j'],'muscle_heat_energy_j':native['muscle_heat_energy_j'],
                'signed_active_fiber_power_w':native['signed_active_fiber_power_w'],'metabolic_reference':copy.deepcopy(native['metabolic_reference']),
                'positive_active_fiber_work_j':native['positive_active_fiber_work_j'],
                'muscle_energy_ledger_basis':'Cumulative chemical and signed active-fiber energy share native endpoint quadrature; heat is their difference; reference is frozen at initialization',
                'muscle_heat_w':native.get('muscle_heat_w'),'metabolic_analysis_mass_kg':native.get('metabolic_analysis_mass_kg'),
                'muscle_metabolic_basis':'Native Umberger/Uchida model with explicit generic composition/mass priors; analysis mass adds no inertia; total demand must not be added directly atop an overlapping basal budget',
                'positive_muscle_work_basis':'Native active-fiber positive mechanical work, endpoint power quadrature over interval; metabolic cost is separately modeled',
                'joints':clone_snapshot_data(native['coordinates']),'native_bodies':clone_snapshot_data(native['bodies']),
                'body_environment':{'kind':native['environment'],'contact_model':copy.deepcopy(native.get('contact_model')),'contacts':copy.deepcopy(native['contacts']),'contact_force_source_n':native['contact_force_n'],'canonical_wrenches':self.registration.contact_wrenches(native),
                    'scope':self._contact_scope(),'segment_contact':self.mechanical_fidelity.get('segment_contact')},
                'audit':{'momentum_balance_residual_source_n':native['momentum_balance_residual_n'],'constraint_position_error':native['constraint_position_error'],
                         'constraint_velocity_error':native['constraint_velocity_error'],'external_work_j':native['external_work_j'],'external_power_w':native['external_power_w']},
                'ownership':{'inertia':'Native source segments exclusively own body mass/inertia','tissues':'Canonical meshes attached to inferred native segment supports; no additional tissue mass integrated','muscle_activation':'Native muscle activation/tendon/fiber states'},
                **({} if self.display_poser is None else {'display_pose':self._display_pose(native)}),
                'mechanical_fidelity':copy.deepcopy(self.mechanical_fidelity),
                'limitations':['Canonical anatomical joint locations and surface continuity remain uncalibrated',*(['Garment partitioned face contact lacks CCD/self/edge/full containment validation'] if self.garments else ['Whole-garment feedback disabled for this plant']),'All-organ volumetric deformation/contact is not implemented by this registration',*([] if self.mechanical_fidelity.get('joint_stops') else ['No coordinate limits are enforced: the model declares ranges, holds zero CoordinateLimitForce, and OpenSim does not clamp during forward dynamics']),*([] if self.mechanical_fidelity.get('tissue_ligaments') else ['No tissue force elements: 645 classified tissue structures exist as geometry and carry no force in this plant'])]}
        if native['environment']=='supine':
            p=np.array([native['support_plane_source_x_m'],0.,0.]);frame['body_environment']['plane']={'point_m':(self.registration.basis@p+self.registration.global_map[:3,3]).tolist(),'normal':self.registration.basis[:,0].tolist(),'basis':'One native ideal plane shared by all contact proxies'}
        contacts=self.registration.cutaneous_contacts(native)
        if contacts is not None:
            frame['cutaneous_contacts']=contacts
            frame['body_environment']['surface_foundation']=copy.deepcopy(native['surface_foundation'])
            frame['body_environment']['scope']='Retained skin quadrature with declared nonlinear material law; selected material sensors preserve source identity. Registration and bed material scope remain explicit in native receipts.'
        frame['gravity_m_s2']=(self.registration.basis@native['gravity_m_s2']).tolist()
        binding=getattr(self,'surface_binding',None)
        # Internal selective observations need no display metadata. Public
        # endpoints carry all supports and expose the exact shared attachment.
        if binding is not None and len(frame['entities'])==len(self.registration.rows):
            frame['surface_binding']=binding.manifest()
            frame['surface_transforms']=binding.frame(frame['entities'])
        if self.garments is not None:frame['garment_mechanics']=self.garments.frame()
        return frame
    def snapshot(self):return clone_snapshot_data(self.state)
    def advance(self,dt_s,forces=(),actuation=None):
        old=self.native.snapshot();mapped=[]
        for f in forces:
            if set(f)!={'id','point_m','force_n'}:raise ValueError('Canonical force requires id, point_m, force_n')
            mapped.append(self.registration.force(f['id'],f['point_m'],f['force_n'],old))
        checkpoint=self.native.checkpoint();garment_checkpoint=None if self.garments is None else self.garments.checkpoint()
        try:
            result=self.native.advance(dt_s,mapped,actuation) if self.garments is None else self.garments.advance(self.native,dt_s,mapped,actuation)[0]
            work=result['positive_active_fiber_work_j']-old['positive_active_fiber_work_j']
            self.state=self._project(result,work);return self.snapshot()
        except BaseException:
            self.native.restore(checkpoint)
            if self.garments is not None:self.garments.restore(garment_checkpoint)
            raise
        finally:self.native.release(checkpoint)
    def body_point(self,**query):
        """Native source frame station query; no canonical-frame conversion."""
        return self.native.body_point(**query)
    def transfer_mass(self,**payload):
        """Explicit native-frame endpoint material port with projected rollback."""
        checkpoint=self.checkpoint()
        try:
            result=self.native.transfer_mass(**payload)
            self.state=self._project(result,self.state['positive_muscle_work_j']);return self.snapshot()
        except BaseException:
            if not self.native.closed:self.restore(checkpoint)
            raise
        finally:
            if not self.native.closed:self.release(checkpoint)
    def checkpoint(self):return {'native':self.native.checkpoint(),'frame':self.snapshot(),'garments':None if self.garments is None else self.garments.checkpoint()}
    def restore(self,checkpoint):
        self.native.restore(checkpoint['native']);self.state=clone_snapshot_data(checkpoint['frame'])
        if self.garments is not None:self.garments.restore(checkpoint['garments'])
    def release(self,checkpoint):self.native.release(checkpoint['native'])
    def close(self):self.native.close()
