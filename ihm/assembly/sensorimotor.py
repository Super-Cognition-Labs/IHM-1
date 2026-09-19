"""Source ankle reflexes plus an explicitly engineered pinned-IBM descending bridge.

Excitation is dimensionless. Native mechanics owns activation, fiber dynamics,
force and work. No cortical label is treated as an identified motor policy.

Two selectable axes, both defaulting to the behaviour this module always had:

  decoder   None / 'rate_gain'  engineered regional-rate gain on the reflex and an
                                effector-gated drive (the original, a dial).
            'recruitment'       a size-principle motor-unit pool (ihm/assembly/
                                recruitment.py, Potvin & Fuglevand 2017 after
                                Fuglevand et al. 1993): segmental and descending
                                input sum at the pool; excitation is the pool's
                                rate-averaged force fraction. No rate gain.
  reflexes  None / 'ankle'      the eight Geyer & Herr 2010 ankle primitives.
            'geyer_herr_2010'   the source's whole stance/swing set: SOL, TA, GAS,
                                VAS, HAM, GLU, HFL (reflexes.py, Table I and
                                Appendix I of the hash-pinned PDF).

`None`/`None` reproduces the original excitations exactly (scripts/
verify_recruitment.py compares against the committed implementation). Only the
model hash differs, because it hashes this file.
"""
from ihm.brain.active_source import DEFAULT_SOURCE,resolve_source
from copy import deepcopy
from dataclasses import asdict, dataclass
import hashlib
import heapq
import json
import math
from pathlib import Path
import numpy as np
from .brain import BodyBrain
from .reflexes import (finite, GEYER_HERR_2010_DOI, GEYER_HERR_2010_TABLE_I, GEYER_HERR_2010_DELAYS_S,
    GEYER_HERR_2010_SOURCE_NOTES, GEYER_HERR_2010_GROUPS, GEYER_HERR_2010_LUMPING, geyer_herr_2010_stimulation)
from .recruitment import MotorUnitPool, PoolParameters, POTVIN_FUGLEVAND_DOI, FUGLEVAND_1993_DOI
from .sensorimotor_catalog import native_muscle_catalog

MUSCLES=tuple(f'{name}_{side}' for side in ('r','l') for name in ('tibant','soleus','gasmed','gaslat'))
SOURCE_URL='https://www.cs.cmu.edu/~cga/tmp-public/song.pdf'
SOURCE_SHA256='8a60efec9ceba5120c1679610f676ab297e8dce887ed4ae6f0890f7e17e13d94'
# The file at SOURCE_URL is Geyer & Herr 2010 (IEEE TNSRE 18:263), despite the URL's
# name; its Table I is where 1.1/.71/1.2/.3 come from.
DECODERS=('rate_gain','recruitment')
REFLEX_SETS=('ankle','geyer_herr_2010')
STANDARD_GRAVITY_M_S2=9.80665
POSTURE_KEY='__posture__'
DECODER_BASIS={
    'rate_gain':'engineered regional-rate gain and effector-gated drive; not identified motor recruitment',
    'recruitment':('size-principle motor-unit pool (Potvin & Fuglevand 2017 rested pool, after Fuglevand et al. 1993): '
        'segmental reflex stimulation and descending drive sum as pool input; excitation is the pool rate-averaged '
        'isometric force fraction; one representative pool for every muscle, not muscle-specific; descending drive '
        'magnitude still carries the engineered cortical_drive_per_hz prior; not identified motor recruitment of this body'),
}
REFLEX_SCOPE={
    'ankle':'all catalog effectors have engineered sensory/descending ports; eight ankle effectors have source reflex primitives; no autonomous walking policy',
    'geyer_herr_2010':('all catalog effectors have engineered sensory/descending ports; Geyer & Herr 2010 stance/swing reflexes '
        'on SOL, TA, GAS, VAS, HAM, GLU, HFL group members (lumped groups split onto plant muscles by declared membership); '
        'no autonomous walking policy'),
}

@dataclass(frozen=True)
class SensorimotorParameters:
    afferent_delay_s: float=.01
    efferent_delay_s: float=.01
    stance_threshold_n: float=5.
    cortical_afferent_hz_per_strain: float=200.
    cortical_afferent_hz_per_normalized_force: float=40.
    cortical_command_hz: float=100.
    cortical_gain_per_hz: float=.005
    cortical_drive_per_hz: float=.002

    def __post_init__(self):
        for key,value in asdict(self).items():finite(value,0,1000,key)
        for key in ('afferent_delay_s','efferent_delay_s'):
            finite(getattr(self,key),0,1,key)

class SensorimotorController:
    """Source catalog effectors, eight ankle reflexes, transactional neural state.

    At t observe the plant, deliver sensory samples after afferent delay, execute
    the preserved IBM ODE over this exchange interval, derive spinal excitations
    and queue them after efferent delay. Outputs at t+dt apply on the NEXT plant
    interval. The exchange discretization adds latency; it is reported, not
    subtracted from measured time. Blocks discard in-flight signals immediately.
    ``descending`` contains per-effector requested drive fractions: these enter
    the corresponding precentral population and gate an explicit rate decoder.
    Additional sensory inputs are already-delayed external receptor endpoints
    from the previous exchange; they share this brain integration and saturation.
    """
    @classmethod
    def from_root(cls,root,*,source_pin=DEFAULT_SOURCE,**kwargs):
        root=Path(root)
        source_pin=resolve_source(root,source_pin)
        paper=root/'data/raw/sensorimotor/geyer_herr_2010.pdf'
        if hashlib.sha256(paper.read_bytes()).hexdigest()!=SOURCE_SHA256:
            raise ValueError('Reflex primary source hash mismatch')
        kwargs.setdefault('muscle_catalog',native_muscle_catalog(root))
        return cls(BodyBrain(json.loads((root/'data/derived/canonical/brain.json').read_text()),root=root,source_pin=source_pin),**kwargs)

    def __init__(self,brain,parameters=None,*,muscle_catalog=None,decoder=None,reflexes=None,pool_parameters=None):
        if not isinstance(brain,BodyBrain) or brain.time_s!=0:raise ValueError('Fresh pinned BodyBrain required')
        self.brain=brain;self.parameters=parameters or SensorimotorParameters()
        if not isinstance(self.parameters,SensorimotorParameters):raise ValueError('Invalid controller parameters')
        for side in ('lh','rh'):
            for region in ('postcentral','precentral'):
                if f'brain-{side}-{region}' not in brain.ids:raise ValueError('Required sensorimotor population missing')
        self.catalog=deepcopy(muscle_catalog or [{'id':k,'side':k[-1],'body_group':'ankle_foot',
            'motor_region':'brain-'+('lh' if k[-1]=='r' else 'rh')+'-precentral',
            'sensory_region':'brain-'+('lh' if k[-1]=='r' else 'rh')+'-postcentral',
            'assignment_basis':'explicit ankle-only engineering fixture'} for k in MUSCLES])
        self.muscles=tuple(row['id'] for row in self.catalog)
        if len(set(self.muscles))!=len(self.muscles) or not set(MUSCLES)<=set(self.muscles):raise ValueError('Catalog must contain unique ankle source effectors')
        self.bindings={row['id']:row for row in self.catalog}
        for row in self.catalog:
            if row['side'] not in ('r','l') or row['motor_region'] not in brain.ids or row['sensory_region'] not in brain.ids:raise ValueError('Invalid sensorimotor assignment')
        self.reference_rates=dict(zip(brain.ids,brain.state[:,1].tolist()))
        self.decoder='rate_gain' if decoder is None else decoder
        self.reflex_set='ankle' if reflexes is None else reflexes
        if self.decoder not in DECODERS:raise ValueError(f'Unknown decoder {decoder!r}')
        if self.reflex_set not in REFLEX_SETS:raise ValueError(f'Unknown reflex set {reflexes!r}')
        if pool_parameters is not None and self.decoder!='recruitment':raise ValueError('pool_parameters require the recruitment decoder')
        # The original behaviour runs through the original code path, untouched.
        self.legacy=self.decoder=='rate_gain' and self.reflex_set=='ankle'
        self.pool=MotorUnitPool(pool_parameters) if self.decoder=='recruitment' else None
        groups=({g:GEYER_HERR_2010_GROUPS[g] for g in ('SOL','TA','GAS')} if self.reflex_set=='ankle' else dict(GEYER_HERR_2010_GROUPS))
        self.reflex_groups={}
        for side in ('r','l'):
            for group,names in groups.items():
                members=tuple(f'{n}_{side}' for n in names if f'{n}_{side}' in self.bindings)
                if not members:raise ValueError(f'Catalog has no member of source reflex group {group} on side {side}')
                self.reflex_groups[(side,group)]=members
        self.reflex_member={m:key for key,members in self.reflex_groups.items() for m in members}
        self.posture=None
        self.reflex_memory={s:{'stance':None,'touchdown_s':None,'theta_takeoff':None} for s in ('r','l')}
        identity={'brain':brain.data,'brain_source_identity':brain.source_identity,'brain_step_s':brain.max_step_s,'parameters':asdict(self.parameters),
                  'catalog':self.catalog,'reflex_source_sha256':SOURCE_SHA256,'implementation_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  'brain_implementation_sha256':hashlib.sha256(Path(__file__).with_name('brain.py').read_bytes()).hexdigest()}
        if not self.legacy:
            identity.update(decoder=self.decoder,reflex_set=self.reflex_set,reflex_groups=sorted([list(k),list(v)] for k,v in self.reflex_groups.items()),
                pool_model_sha256=None if self.pool is None else self.pool.model_sha256,
                reflex_implementation_sha256=hashlib.sha256(Path(__file__).with_name('reflexes.py').read_bytes()).hexdigest(),
                recruitment_implementation_sha256=hashlib.sha256(Path(__file__).with_name('recruitment.py').read_bytes()).hexdigest())
        self.model_sha256=hashlib.sha256(json.dumps(identity,sort_keys=True,allow_nan=False).encode()).hexdigest()
        self.time_s=0.;self.events=[];self.serial=0
        self.arrived={};self.excitations={k:0. for k in self.muscles}

    def _blocks(self,value,label):
        if not isinstance(value,(list,tuple,set)) or any(not isinstance(k,str) for k in value) or set(value)-set(self.muscles):
            raise ValueError(f'Unknown {label} muscle')
        return set(value)

    def _validate_observation(self,observation):
        if not isinstance(observation,dict):raise ValueError('Mechanical observation required')
        time=finite(observation.get('time_s'),0,1e9,'mechanical time')
        if abs(time-self.time_s)>1e-8:raise ValueError('Mechanics and neural clocks differ')
        muscles=observation.get('muscles');contact=observation.get('foot_contact_force_n')
        if not isinstance(muscles,dict) or not isinstance(contact,dict):raise ValueError('Muscles and actual foot contact required')
        loads={s:finite(contact.get(s),0,1e9,'foot normal force') for s in ('r','l')}
        sensors={}
        for key in self.muscles:
            m=muscles.get(key)
            if not isinstance(m,dict) or not isinstance(m.get('sensor_basis'),str) or not m['sensor_basis']:
                raise ValueError(f'Missing sensor provenance: {key}')
            opt=finite(m.get('optimal_fiber_length_m'),1e-9,100,'optimal fiber length')
            fmax=finite(m.get('max_isometric_force_n'),1e-9,1e9,'maximum force')
            # A proxy must remain explicitly named. Never reconstruct CE from path silently.
            field='fiber_length_m' if 'fiber_length_m' in m else 'fiber_length_proxy_m'
            length=finite(m.get(field),1e-9,100,'fiber length')
            force=finite(m.get('tendon_force_n'),0,1e9,'tendon force')
            sensors[key]={'length':length/opt,'force':force/fmax,'force_n':force,'fmax_n':fmax,
                          'stance':loads[self.bindings[key]['side']]>self.parameters.stance_threshold_n,
                          'sensor_basis':m['sensor_basis'],'length_field':field}
        return sensors

    def _posture(self,observation):
        """Non-muscle afference the Geyer & Herr hip/knee laws need, read from the frame.

        Trunk pitch is the torso long axis from vertical, signed along the pelvis
        heading (+ = forward lean), and its rate is the torso angular velocity about
        the heading's lateral axis: the vestibular/trunk signal of the source's
        theta. Knee phi_k = pi - knee_angle (source: 180 deg = straight). Leg load is
        the native foot NORMAL force over body weight (the source's |F_leg| is the
        GRF magnitude; the tangential part is not observed here). Ground +y is up.
        These signals are not blockable by muscle sensory blocks.
        """
        bodies=observation.get('native_bodies',observation.get('bodies'))
        joints=observation.get('joints',observation.get('coordinates'))
        mass=observation.get('effective_native_body_mass_kg',observation.get('mass_kg'))
        if not isinstance(bodies,dict) or not isinstance(joints,dict):raise ValueError('Geyer-Herr reflexes need native bodies and joints')
        weight=finite(mass,1,1000,'body mass')*STANDARD_GRAVITY_M_S2
        def rotation(name):
            try:r=np.asarray(bodies[name]['transform_ground'],float)[:3,:3]
            except (KeyError,TypeError,ValueError,IndexError):raise ValueError(f'Missing {name} transform') from None
            if r.shape!=(3,3) or not np.isfinite(r).all() or not np.allclose(r.T@r,np.eye(3),atol=1e-6):raise ValueError(f'Invalid {name} rotation')
            return r
        torso=rotation('torso');pelvis=rotation('pelvis')
        omega=np.asarray(bodies['torso'].get('angular_velocity_rad_s'),float)
        if omega.shape!=(3,) or not np.isfinite(omega).all():raise ValueError('Missing torso angular velocity')
        up=np.array([0.,1.,0.]);heading=pelvis[:,0]-up*pelvis[1,0];norm=np.linalg.norm(heading)
        if norm<1e-6:raise ValueError('Pelvis heading undefined (not upright)')
        heading/=norm;lateral=np.cross(heading,up)
        axis=torso[:,1]
        theta=float(math.atan2(axis@heading,axis@up));theta_rate=float(-(omega@lateral))
        contact=observation['foot_contact_force_n']
        sample={'time_s':self.time_s,'theta':theta,'theta_rate':theta_rate,'body_weight_n':weight,'stance':{},'load_bw':{},'knee_phi':{},'knee_phi_rate':{},
                'basis':'native torso/pelvis transforms and angular velocity; knee coordinate; native foot normal force / (mass x 9.80665)'}
        for side in ('r','l'):
            knee=joints.get(f'knee_angle_{side}')
            if not isinstance(knee,dict):raise ValueError(f'Missing knee_angle_{side}')
            load=finite(contact.get(side),0,1e9,'foot normal force')
            sample['stance'][side]=load>self.parameters.stance_threshold_n
            sample['load_bw'][side]=load/weight
            sample['knee_phi'][side]=math.pi-finite(knee.get('value'),-10,10,'knee angle')
            sample['knee_phi_rate'][side]=-finite(knee.get('speed'),-1e4,1e4,'knee speed')
        return sample

    def _arrive_posture(self,sample):
        """Delayed posture arrives: update touchdown/take-off memory, in arrival order."""
        first=self.posture is None
        for side in ('r','l'):
            m=self.reflex_memory[side];now=sample['stance'][side]
            if not first:
                if now and not m['stance']:m['touchdown_s']=sample['time_s']
                if m['stance'] and not now:m['theta_takeoff']=sample['theta']
            m['stance']=now
        self.posture=sample

    def _trailing(self,side):
        """DSup: this leg is the TRAILING leg of a double support (stance began first)."""
        other='l' if side=='r' else 'r';q=self.posture
        if q is None or not(q['stance'][side] and q['stance'][other]):return False
        mine=self.reflex_memory[side]['touchdown_s'];theirs=self.reflex_memory[other]['touchdown_s']
        if mine is None and theirs is None:return False  # order unknown: both in stance since before observation
        if mine is None:return True
        if theirs is None:return False
        return mine<theirs

    def _source_stimulations(self):
        """Source stimulation S (unlimited) and S0 per reflex muscle, from delayed afference only."""
        out={}
        for side in ('r','l'):
            force={};length={}
            for (s_,group),members in self.reflex_groups.items():
                if s_!=side:continue
                samples=[self.arrived.get(m) for m in members]
                if all(samples):
                    fmax=sum(x['fmax_n'] for x in samples)
                    force[group]=sum(x['force_n'] for x in samples)/fmax
                    length[group]=sum(x['fmax_n']*x['length'] for x in samples)/fmax
            q=self.posture
            state={}
            if q is not None:
                other='l' if side=='r' else 'r'
                state=dict(load_bw=q['load_bw'][side],contra_load_bw=q['load_bw'][other],dsup=self._trailing(side),
                    knee_phi=q['knee_phi'][side],knee_phi_rate=q['knee_phi_rate'][side],theta=q['theta'],theta_rate=q['theta_rate'],
                    theta_takeoff=self.reflex_memory[side]['theta_takeoff'])
            for (s_,group),members in self.reflex_groups.items():
                if s_!=side:continue
                if group in ('SOL','TA','GAS'):
                    # Ankle laws switch on the muscle sample's own delayed stance flag,
                    # exactly as the retained eight primitives always did.
                    own=[self.arrived.get(m) for m in members]
                    if not all(own):result=None
                    else:result=geyer_herr_2010_stimulation(group,stance=own[0]['stance'],force=force,length=length,side_state=state)
                elif q is None:result=None
                else:result=geyer_herr_2010_stimulation(group,stance=q['stance'][side],force=force,length=length,side_state=state)
                for m in members:out[m]=result
        return out

    def _queue(self,arrival,kind,key,value):
        self.serial+=1;heapq.heappush(self.events,(arrival,self.serial,kind,key,deepcopy(value)))

    def step(self,dt_s,mechanical_observation,*,descending=None,sensory_blocks=(),motor_blocks=(),physiology=None,additional_sensory_inputs_hz=None):
        dt=finite(dt_s,1e-9,.1,'dt_s');sensors=self._validate_observation(mechanical_observation)
        sb=self._blocks(sensory_blocks,'sensory block');mb=self._blocks(motor_blocks,'motor block')
        descending={} if descending is None else descending
        if not isinstance(descending,dict) or set(descending)-set(self.muscles):raise ValueError('Unknown descending effector')
        descending={k:finite(v,0,1,'descending drive fraction') for k,v in descending.items()}
        additional={} if additional_sensory_inputs_hz is None else additional_sensory_inputs_hz
        if not isinstance(additional,dict) or set(additional)-set(self.brain.ids):raise ValueError('Unknown additional sensory population')
        additional={k:finite(v,0,1000,'additional sensory rate Hz') for k,v in additional.items()}
        posture=self._posture(mechanical_observation) if getattr(self,'reflex_set','ankle')=='geyer_herr_2010' else None
        saved=self.checkpoint()
        try:return self._step(dt,sensors,descending,sb,mb,physiology,additional,posture)
        except Exception:
            self.restore(saved)
            raise

    def _step(self,dt,sensors,descending,sb,mb,physiology,additional,posture=None):
        p=self.parameters;end=self.time_s+dt
        if end<=self.time_s:raise ValueError('Clock must advance')
        self.events=[e for e in self.events if not(e[2]=='sensor' and e[3] in sb or e[2]=='motor' and e[3] in mb)]
        heapq.heapify(self.events)
        for k in sb:self.arrived.pop(k,None)
        for k in mb:self.excitations[k]=0.
        for k,v in sensors.items():
            if k not in sb:self._queue(self.time_s+p.afferent_delay_s,'sensor',k,v)
        if posture is not None:self._queue(self.time_s+p.afferent_delay_s,'posture',POSTURE_KEY,posture)
        # Brain integrates the previously arrived afferents (explicit causal exchange).
        inputs={}
        for key,s in self.arrived.items():
            region=self.bindings[key]['sensory_region']
            rate=p.cortical_afferent_hz_per_strain*max(0,s['length']-1)+p.cortical_afferent_hz_per_normalized_force*s['force']
            inputs[region]=min(1000.,inputs.get(region,0)+rate)
        for key,value in descending.items():
            region=self.bindings[key]['motor_region']
            inputs[region]=min(1000.,inputs.get(region,0)+p.cortical_command_hz*value)
        # Caller supplies already-delayed receptor output from the previous exchange.
        # Pool here so the existing shared brain advances exactly once.
        for region,rate in additional.items():
            inputs[region]=min(1000.,inputs.get(region,0)+rate)
        neural=self.brain.step(dt,physiology=physiology,sensory_inputs_hz=inputs)
        while self.events and self.events[0][0]<=end+1e-12:
            _,_,kind,key,value=heapq.heappop(self.events)
            if kind=='sensor':self.arrived[key]=value
            elif kind=='posture':self._arrive_posture(value)
            else:self.excitations[key]=value
        rates=dict(zip(neural['regional_state']['node_ids'],neural['regional_state']['activity_hz']))
        targets={};gains={};drives={}
        if not self.legacy:return self._decode(end,dt,sensors,descending,sb,mb,additional,inputs,neural,rates)
        for key in self.muscles:
            side=self.bindings[key]['side'];region=self.bindings[key]['motor_region']
            gains[key]=float(np.clip(1+p.cortical_gain_per_hz*(rates[region]-self.reference_rates[region]),0,2))
            drives[key]=descending.get(key,0)*float(np.clip(p.cortical_drive_per_hz*rates[region],0,1))
            sensor=self.arrived.get(key);reflex=0.
            if sensor and key in MUSCLES:
                if key.startswith('tibant'):
                    reflex=1.1*(sensor['length']-.71)
                    soleus=self.arrived.get('soleus_'+side)
                    if soleus and sensor['stance']:reflex-=.3*soleus['force']
                elif sensor['stance']:
                    if key.startswith('soleus'):reflex=1.2*sensor['force']
                    else:
                        group=[self.arrived.get(n+'_'+side) for n in ('gasmed','gaslat')]
                        # Both heads required for the explicitly lumped GAS transfer.
                        if all(group):reflex=1.1*sum(s['force_n'] for s in group)/sum(s['fmax_n'] for s in group)
            baseline=.01 if key in MUSCLES else 0.
            target=float(np.clip(baseline+gains[key]*reflex+drives[key],baseline,1))
            targets[key]=0. if key in mb else target
            if key not in mb:self._queue(end+p.efferent_delay_s,'motor',key,target)
        # Zero-delay motor commands are available at this endpoint only.
        while self.events and self.events[0][0]<=end+1e-12:
            _,_,kind,key,value=heapq.heappop(self.events)
            if kind=='motor':self.excitations[key]=value
            elif kind=='posture':self._arrive_posture(value)
            else:self.arrived[key]=value
        self.time_s=end
        return {'schema':'ihm.sensorimotor.v1','time_s':end,'motor_excitations':dict(self.excitations),
            'requested_excitations':targets,'sensors':deepcopy(sensors),'delayed_sensors':deepcopy(self.arrived),
            'brain':neural,'additional_sensory_inputs_hz':dict(additional),'brain_sensory_inputs_hz':dict(inputs),'descending_gain':gains,'descending_drive':drives,'pending_events':len(self.events),
            'sensory_blocks':sorted(sb),'motor_blocks':sorted(mb),'exchange_interval_s':dt,
            'neural_delay_s':p.afferent_delay_s+p.efferent_delay_s,'activation_owner':'mechanical_plant',
            'biological_validation':False,'decoder_basis':'engineered regional-rate gain and effector-gated drive; not identified motor recruitment',
            'scope':'all catalog effectors have engineered sensory/descending ports; eight ankle effectors have source reflex primitives; no autonomous walking policy',
            'spinal_reflex_effectors':list(MUSCLES),'descending_effector_count':len(self.muscles),
            'model_sha256':self.model_sha256}

    def _decode(self,end,dt,sensors,descending,sb,mb,additional,inputs,neural,rates):
        """Non-legacy decode: source stimulation per reflex muscle, then the selected decoder."""
        p=self.parameters;floor=GEYER_HERR_2010_TABLE_I['S_min']
        stimulations=self._source_stimulations()
        targets={};gains={};drives={};segmental={};total={};recruited={};source={}
        for key in self.muscles:
            region=self.bindings[key]['motor_region']
            drives[key]=descending.get(key,0)*float(np.clip(p.cortical_drive_per_hz*rates[region],0,1))
            reflex_muscle=key in self.reflex_member
            result=stimulations.get(key)
            source[key]=None if result is None else float(result[0])
            if self.decoder=='rate_gain':
                gains[key]=float(np.clip(1+p.cortical_gain_per_hz*(rates[region]-self.reference_rates[region]),0,2))
                if reflex_muscle:
                    # Same shape as the original ankle rule: the dial scales the reflex
                    # term above its prestimulation, the result is limited to 0.01..1.
                    s0=floor if result is None else result[1];reflex=0. if result is None else result[0]-result[1]
                    target=float(np.clip(s0+gains[key]*reflex+drives[key],floor,1))
                else:target=float(np.clip(drives[key],0,1))
            else:
                gains[key]=1.
                s=0. if not reflex_muscle else (floor if result is None else float(np.clip(result[0],floor,1)))
                excitation,base,summed=self.pool.combine(s,drives[key])
                segmental[key]=float(base);total[key]=float(summed);recruited[key]=int(self.pool.recruited(summed))
                target=float(np.clip(excitation,floor if reflex_muscle else 0.,1))
            targets[key]=0. if key in mb else target
            if key not in mb:self._queue(end+p.efferent_delay_s,'motor',key,target)
        while self.events and self.events[0][0]<=end+1e-12:
            _,_,kind,key,value=heapq.heappop(self.events)
            if kind=='motor':self.excitations[key]=value
            elif kind=='posture':self._arrive_posture(value)
            else:self.arrived[key]=value
        self.time_s=end
        reflex_muscles=sorted(self.reflex_member)
        frame={'schema':'ihm.sensorimotor.v1','time_s':end,'motor_excitations':dict(self.excitations),
            'requested_excitations':targets,'sensors':deepcopy(sensors),'delayed_sensors':deepcopy(self.arrived),
            'brain':neural,'additional_sensory_inputs_hz':dict(additional),'brain_sensory_inputs_hz':dict(inputs),'descending_gain':gains,'descending_drive':drives,'pending_events':len(self.events),
            'sensory_blocks':sorted(sb),'motor_blocks':sorted(mb),'exchange_interval_s':dt,
            'neural_delay_s':p.afferent_delay_s+p.efferent_delay_s,'activation_owner':'mechanical_plant',
            'excitation_owner':'this controller produces excitation only; native mechanics owns activation, fibre dynamics and force',
            'biological_validation':False,'decoder':self.decoder,'reflex_set':self.reflex_set,
            'decoder_basis':DECODER_BASIS[self.decoder],'scope':REFLEX_SCOPE[self.reflex_set],
            'spinal_reflex_effectors':reflex_muscles,'descending_effector_count':len(self.muscles),
            'source_stimulations':source,
            'reflex_source':{'doi':GEYER_HERR_2010_DOI,'sha256':SOURCE_SHA256,'values':'Table I p.269','laws':'Appendix I p.270',
                'groups':{f'{g}_{s_}':list(m) for (s_,g),m in self.reflex_groups.items()},'lumping':GEYER_HERR_2010_LUMPING,
                'source_delay_s':{g:GEYER_HERR_2010_DELAYS_S[g] for g in sorted({g for _,g in self.reflex_groups})},
                'realised_loop_delay_s':p.afferent_delay_s+p.efferent_delay_s,
                'delay_basis':'one afferent+efferent split for every group, plus up to one exchange interval; source VAS 10 ms and hip 5 ms are NOT honoured separately',
                'notes':list(GEYER_HERR_2010_SOURCE_NOTES)},
            'model_sha256':self.model_sha256}
        if self.reflex_set=='geyer_herr_2010':
            frame['reflex_state']={'memory':deepcopy(self.reflex_memory),'trailing':{s_:self._trailing(s_) for s_ in ('r','l')},
                'posture':deepcopy(self.posture)}
        if self.pool is not None:
            frame['recruitment']={'model_sha256':self.pool.model_sha256,'sources':{'parameterisation':POTVIN_FUGLEVAND_DOI,'original':FUGLEVAND_1993_DOI},
                'segmental_input_fraction':segmental,'total_input_fraction':total,'recruited_units':recruited,'units_per_pool':self.pool.parameters.units,
                'input_basis':'fractions of Emax; descending input = request x clip(cortical_drive_per_hz x regional rate), the same engineered factor the rate decoder uses'}
        return frame

    def checkpoint(self):
        if not getattr(self,'legacy',True):
            return {'schema':'ihm.sensorimotor-state.v1','model_sha256':self.model_sha256,'time_s':self.time_s,
                'brain_state':self.brain.state.tolist(),'arrived':deepcopy(self.arrived),
                'excitations':dict(self.excitations),'events':deepcopy([list(e) for e in sorted(self.events)]),'serial':self.serial,
                'posture':deepcopy(self.posture),'reflex_memory':deepcopy(self.reflex_memory)}
        return {'schema':'ihm.sensorimotor-state.v1','model_sha256':self.model_sha256,'time_s':self.time_s,
                'brain_state':self.brain.state.tolist(),'arrived':deepcopy(self.arrived),
                'excitations':dict(self.excitations),'events':deepcopy([list(e) for e in sorted(self.events)]),'serial':self.serial}

    def restore(self,checkpoint):
        c=deepcopy(checkpoint)
        if not isinstance(c,dict) or c.get('schema')!='ihm.sensorimotor-state.v1' or c.get('model_sha256')!=self.model_sha256:
            raise ValueError('Checkpoint identity mismatch')
        time=finite(c.get('time_s'),0,1e9,'checkpoint time');brain=np.asarray(c.get('brain_state'),float)
        if brain.shape!=self.brain.state.shape or not np.isfinite(brain).all():raise ValueError('Invalid brain checkpoint')
        extended_posture=getattr(self,'reflex_set','ankle')=='geyer_herr_2010'
        def posture(q):
            if not isinstance(q,dict):raise ValueError('Invalid posture checkpoint')
            finite(q.get('time_s'),0,1e9,'posture time');finite(q.get('theta'),-10,10,'theta');finite(q.get('theta_rate'),-1e4,1e4,'theta rate')
            finite(q.get('body_weight_n'),1,1e5,'body weight')
            for field,lo,hi in (('load_bw',0,1e6),('knee_phi',-20,20),('knee_phi_rate',-1e4,1e4)):
                if not isinstance(q.get(field),dict) or set(q[field])!={'r','l'}:raise ValueError('Invalid posture sides')
                for v in q[field].values():finite(v,lo,hi,field)
            if not isinstance(q.get('stance'),dict) or set(q['stance'])!={'r','l'} or any(type(v) is not bool for v in q['stance'].values()):raise ValueError('Invalid posture stance')
        def sensor(s):
            if not isinstance(s,dict) or type(s.get('stance')) is not bool:raise ValueError('Invalid sensor checkpoint')
            for k in ('length','fmax_n'):finite(s.get(k),1e-12,1e12,k)
            for k in ('force','force_n'):finite(s.get(k),0,1e12,k)
            if not isinstance(s.get('sensor_basis'),str) or s.get('length_field') not in ('fiber_length_m','fiber_length_proxy_m'):
                raise ValueError('Invalid sensor provenance')
        arrived=c.get('arrived');exc=c.get('excitations')
        if not isinstance(arrived,dict) or set(arrived)-set(self.muscles) or not isinstance(exc,dict) or set(exc)!=set(self.muscles):raise ValueError('Invalid effector state')
        for s in arrived.values():sensor(s)
        for v in exc.values():finite(v,0,1,'excitation')
        serial=c.get('serial');events=c.get('events')
        if type(serial) is not int or serial<0 or not isinstance(events,list):raise ValueError('Invalid event state')
        previous=(time,-1);seen=set();validated=[]
        for e in events:
            if not isinstance(e,list) or len(e)!=5:raise ValueError('Invalid event')
            t,n,kind,key,value=e;finite(t,time,1e9,'arrival')
            posture_event=extended_posture and kind=='posture' and key==POSTURE_KEY
            if type(n) is not int or not 0<n<=serial or n in seen or (t,n)<=previous or not posture_event and (key not in self.muscles or kind not in ('sensor','motor')):raise ValueError('Invalid event order/owner')
            seen.add(n);previous=(t,n)
            if kind=='sensor':sensor(value)
            elif posture_event:posture(value)
            else:finite(value,0,1,'queued excitation')
            validated.append(tuple(e))
        extended=not getattr(self,'legacy',True)
        if extended:
            if set(c)!={'schema','model_sha256','time_s','brain_state','arrived','excitations','events','serial','posture','reflex_memory'}:raise ValueError('Extended checkpoint fields differ')
            if c['posture'] is not None:
                if not extended_posture:raise ValueError('Posture state without posture reflexes')
                posture(c['posture'])
            memory=c['reflex_memory']
            if not isinstance(memory,dict) or set(memory)!={'r','l'}:raise ValueError('Invalid reflex memory')
            for m in memory.values():
                if not isinstance(m,dict) or set(m)!={'stance','touchdown_s','theta_takeoff'}:raise ValueError('Invalid reflex memory')
                if m['stance'] is not None and type(m['stance']) is not bool:raise ValueError('Invalid reflex memory stance')
                if m['touchdown_s'] is not None:finite(m['touchdown_s'],0,1e9,'touchdown')
                if m['theta_takeoff'] is not None:finite(m['theta_takeoff'],-10,10,'take-off lean')
        self.time_s=time;self.brain.time_s=time;self.brain.state=brain.copy()
        self.arrived=arrived;self.excitations=exc;self.serial=serial;self.events=validated
        if extended:self.posture=c['posture'];self.reflex_memory=c['reflex_memory']
