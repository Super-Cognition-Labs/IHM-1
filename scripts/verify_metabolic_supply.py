"""Known answers for metabolic supply as a limit on excitation (docs/METABOLIC_SUPPLY.md).

Default: local fixtures through the real EmbodiedRuntime.step, no native process.
``--native``: additionally drives the real native mechanical plant (one process,
prlimit-capped by its own launcher) to check that the cap reaches the engine.

Known answers:
  K1  supply >= demand  -> delivered excitation identical to requested, bit for bit,
      and every frame identical to the pre-change embodied.py (pinned blob).
  K2  supply = 0        -> every muscle's delivered excitation is exactly 0 from the
      next interval on, and that interval's active work increment is exactly 0.
  K3  the ledger raw = drawn + unsupplied + pending closes to round-off on EVERY
      exchange, and chemical = heat + work on the drawn and unsupplied triples.
  K4  idempotence: each stateful function called twice at the same input.
  K5  controls that CAN fail, made to fail on purpose.
  K6  None still aborts on unmet energy, with the original message.
"""
from copy import deepcopy
from pathlib import Path
import argparse,json,math,subprocess,sys,types,unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from ihm.assembly.embodied import EmbodiedRuntime
from ihm.assembly.metabolic_supply import (J_PER_KCAL,SupplyLimiter,IntervalCap,CappedActuationPlant,
    assess_exchange,ledger_closes)

# git blob of ihm/assembly/embodied.py immediately before metabolic supply existed
# (commit eea04112ee3ead067790a15e65d1b10f6e06cc61). A blob id never moves.
PRE_SUPPLY_EMBODIED_BLOB='542f566aec1d1fdfdd04fe9f754c0402ae9743cd'
DT=.02
REPORT={}


def pre_supply_runtime_class():
    source=subprocess.check_output(['git','-C',str(ROOT),'cat-file','blob',PRE_SUPPLY_EMBODIED_BLOB])
    module=types.ModuleType('ihm.assembly._embodied_pre_supply')
    module.__package__='ihm.assembly';module.__file__='embodied@'+PRE_SUPPLY_EMBODIED_BLOB
    exec(compile(source,module.__file__,'exec'),module.__dict__)
    return module.EmbodiedRuntime


class Plant:
    """Holds each muscle's last command, as the native engine does.

    Active work and heat above the zero-excitation reference are proportional to
    total held excitation, so zero excitation gives exactly the reference.
    Chemical = work + heat exactly, as the Umberger ledger does.
    """
    def __init__(self,held=None):
        self.t=0.;self.held=dict(held or {'a':.2,'b':.2,'c':.2});self.m=0.;self.w=0.;self.h=0.;self.commands=[]
    def snapshot(self):
        return {'time_s':self.t,'entities':{},'muscles':{k:{'excitation':v} for k,v in self.held.items()},
                'muscle_metabolic_energy_j':self.m,'signed_active_fiber_work_j':self.w,'muscle_heat_energy_j':self.h,
                'metabolic_reference':{'M0_w':100.,'W0_w':0.,'H0_w':100.}}
    def checkpoint(self):return deepcopy(self.__dict__)
    def restore(self,state):self.__dict__=deepcopy(state)
    def advance(self,dt_s,forces=(),actuation=None):
        self.commands.append(None if actuation is None else dict(actuation))
        for k,v in (actuation or {}).items():
            if k not in self.held or not 0<=v<=1:raise ValueError('bad excitation')
            self.held[k]=v
        total=sum(self.held.values())
        self.t+=dt_s;self.w+=40.*total*dt_s;self.h+=(100.+60.*total)*dt_s;self.m=self.w+self.h
        return {**self.snapshot(),'positive_muscle_work_j':max(0.,40.*total*dt_s)}
    def close(self):pass


class Neural:
    """Commands a and b only; c is never commanded, so the plant holds it."""
    def __init__(self,excitation=.5):self.t=0.;self.excitation=excitation
    def checkpoint(self):return self.t
    def restore(self,state):self.t=state
    def step(self,dt,observation,**inputs):
        self.t+=dt;return {'time_s':self.t,'motor_excitations':{'a':self.excitation,'b':self.excitation}}


class Native:
    """Muscle request = basal share + signed increment; supply capacity is scripted."""
    BASAL_W=20.
    def __init__(self,capacity_w=None):
        self.t=0.;self.capacity_w=capacity_w or (lambda step:1e9);self.demands=[];self.last={}
    def snapshot(self):
        return {'elapsed_s':self.t,'time_s':100+self.t,'values':{'mean_arterial_pressure_mmhg':90,'oxygen_saturation':.98,
            'core_temperature_c':37,'lung_volume_ml':3000,**self.last}}
    def respiratory_load(self,p):pass
    def signed_step(self,reference,m,h,w):
        self.demands.append((m,h,w))
        requested=(self.BASAL_W+m)*DT/J_PER_KCAL
        capacity=self.capacity_w(len(self.demands))*DT/J_PER_KCAL
        self.last={'coupling.muscle_requested_kcal':requested,'coupling.muscle_unmet_kcal':max(0.,requested-capacity)}
        self.t+=DT;return self.snapshot()
    def close(self,graceful=True):pass


class Load:
    def project_load(self,forces,entities,volume):return {'external_pressure_pa':0.,'force_ports':[],'ignored_nonrespiratory_ids':[]}
    def geometry(self,volume,entities,time):return {'entities':{},'skin_field':{},'time_s':time}


class Exchange:
    def observe(self,snapshot):return {'time_s':snapshot['time_s']}


def body(mode=None,capacity_w=None,cls=EmbodiedRuntime,excitation=.5):
    kwargs={} if mode is None else {'metabolic_supply':mode}
    return cls(Plant(),Neural(excitation),Native(capacity_w),Exchange(),Load(),**kwargs)


def run(runtime,steps):
    return [runtime.step({}) for _ in range(steps)]


class LawTests(unittest.TestCase):
    def test_full_supply_is_coverage_one_exactly(self):
        a=assess_exchange(30.,20.,10.,DT,requested_kcal=(20+30)*DT/J_PER_KCAL,unmet_kcal=0.)
        self.assertEqual(a['coverage'],1.)
        self.assertEqual(a['drawn_j'],a['sent_j']);self.assertEqual(a['unsupplied_j'],{'m_j':0.,'h_j':0.,'w_j':0.})
    def test_zero_supply_is_coverage_zero_exactly(self):
        requested=(20+30)*DT/J_PER_KCAL
        a=assess_exchange(30.,20.,10.,DT,requested_kcal=requested,unmet_kcal=requested)
        self.assertEqual(a['coverage'],0.);self.assertEqual(a['drawn_j']['m_j'],0.)
        self.assertAlmostEqual(a['native_basal_unmet_j'],20*DT,places=12)
        REPORT['law_zero_supply']={'coverage':a['coverage'],'unsupplied_m_j':a['unsupplied_j']['m_j'],'basal_unmet_j':a['native_basal_unmet_j']}
    def test_partial_supply_is_the_covered_fraction(self):
        a=assess_exchange(40.,30.,10.,DT,requested_kcal=(20+40)*DT/J_PER_KCAL,unmet_kcal=10*DT/J_PER_KCAL)
        self.assertAlmostEqual(a['coverage'],.75,places=12);self.assertEqual(a['native_basal_unmet_j'],0.)
        for t in (a['drawn_j'],a['unsupplied_j']):self.assertTrue(ledger_closes(t['m_j'],t['w_j'],t['h_j']))
    def test_decrement_is_never_charged_the_basal_deficit(self):
        a=assess_exchange(-5.,-5.,0.,DT,requested_kcal=15*DT/J_PER_KCAL,unmet_kcal=3*DT/J_PER_KCAL)
        self.assertEqual(a['coverage'],1.);self.assertEqual(a['unsupplied_j']['m_j'],0.)
        self.assertAlmostEqual(a['native_basal_unmet_j'],3*DT,places=12)
    def test_domain(self):
        for kwargs in ({'requested_kcal':1e-4,'unmet_kcal':-1e-9},{'requested_kcal':1e-4,'unmet_kcal':2e-4},
                       {'requested_kcal':float('nan'),'unmet_kcal':0.}):
            with self.subTest(**{k:str(v) for k,v in kwargs.items()}),self.assertRaises(ValueError):
                assess_exchange(1.,1.,0.,DT,**kwargs)
        with self.assertRaises(ValueError):assess_exchange(1.,.5,0.,DT,requested_kcal=1.,unmet_kcal=0.)


class RuntimeTests(unittest.TestCase):
    def test_K1_none_reproduces_pre_supply_embodied_exactly(self):
        Old=pre_supply_runtime_class()
        old=body(cls=Old);new=body()
        a=run(old,25);b=run(new,25)
        self.assertEqual(a,b);self.assertEqual(old.plant.commands,new.plant.commands)
        self.assertEqual(old.native.demands,new.native.demands)
        self.assertNotIn('metabolic_supply',b[-1]['coupling'])
        REPORT['K1_none_vs_pinned']={'frames_compared':len(a),'identical':a==b}
    def test_K1_ample_supply_delivers_request_bit_for_bit(self):
        plain=body();capped=body('energy_fraction')
        a=run(plain,25);b=run(capped,25)
        self.assertEqual(plain.plant.commands,capped.plant.commands)
        self.assertEqual(plain.native.demands,capped.native.demands)
        for x,y in zip(a,b):
            y=deepcopy(y);supply=y['coupling'].pop('metabolic_supply')
            self.assertEqual(supply['excitation_scale_applied'],1.)
            y['coupling']['metabolic_law']=x['coupling']['metabolic_law']
            self.assertEqual(x,y)
        self.assertEqual(b[-1]['coupling']['metabolic_supply']['cumulative']['unsupplied_j']['m_j'],0.)
        REPORT['K1_ample']={'exchanges':len(b),'commands_identical':plain.plant.commands==capped.plant.commands,
            'max_scale':max(f['coupling']['metabolic_supply']['excitation_scale_applied'] for f in b),
            'min_coverage':min(f['coupling']['metabolic_supply']['coverage_next'] for f in b)}
    def test_K2_zero_supply_zero_excitation_zero_active_work(self):
        runtime=body('energy_fraction',capacity_w=lambda step:0.)
        frames=run(runtime,6)
        first=frames[0]['coupling']['metabolic_supply']
        self.assertEqual(first['excitation_scale_applied'],1.);self.assertEqual(first['coverage_next'],0.)
        for frame,command in zip(frames[1:],runtime.plant.commands[1:]):
            s=frame['coupling']['metabolic_supply']
            self.assertEqual(s['excitation_scale_applied'],0.)
            self.assertEqual(command,{'a':0.,'b':0.,'c':0.})
            self.assertEqual(frame['coupling']['signed_work_increment_w'],0.)
            # Chemical increment is round-off, not 0: the fixture differences cumulative
            # float energies (first run measured 4.26e-14 W). Active work above is exact.
            self.assertLessEqual(abs(frame['coupling']['native_extra_metabolic_demand_w']),1e-12)
        self.assertEqual(runtime.plant.held,{'a':0.,'b':0.,'c':0.})
        last=frames[-1]['coupling']['metabolic_supply']
        self.assertGreater(last['cumulative']['unsupplied_j']['m_j'],0.)
        self.assertFalse(last['physically_conserved'])
        REPORT['K2_zero_supply']={'commands_after_first':runtime.plant.commands[1:],
            'active_work_increment_w_after_first':[f['coupling']['signed_work_increment_w'] for f in frames[1:]],
            'chemical_increment_w_after_first':[f['coupling']['native_extra_metabolic_demand_w'] for f in frames[1:]],
            'cumulative_unsupplied_j':last['cumulative']['unsupplied_j'],'cumulative_drawn_j':last['cumulative']['drawn_j'],
            'physically_conserved':last['physically_conserved']}
    def test_K3_ledger_closes_every_exchange(self):
        schedule=lambda step:(1e9,60.,20.,0.,35.,1e9)[(step//7)%6]
        runtime=body('energy_fraction',capacity_w=schedule)
        frames=run(runtime,84)
        raw={'m_j':0.,'h_j':0.,'w_j':0.};worst=0.;scales=set()
        for frame in frames:
            c=frame['coupling'];s=c['metabolic_supply'];scales.add(s['excitation_scale_applied'])
            raw['m_j']+=c['native_extra_metabolic_demand_w']*DT
            raw['h_j']+=c['muscle_heat_increment_w']*DT;raw['w_j']+=c['signed_work_increment_w']*DT
            cum=s['cumulative'];pending=c['metabolic_pending_energy_j']
            for key in raw:
                residual=raw[key]-cum['drawn_j'][key]-cum['unsupplied_j'][key]-pending[key]
                worst=max(worst,abs(residual),abs(s['closure_residual_j'][key]))
                self.assertLessEqual(abs(residual),1e-12*(1+cum['magnitude_j']))
            for channel in ('drawn_j','unsupplied_j'):
                t=s['interval'][channel];self.assertTrue(ledger_closes(t['m_j'],t['w_j'],t['h_j']))
        self.assertGreater(len(scales),2)
        REPORT['K3_ledger']={'exchanges':len(frames),'worst_closure_residual_j':worst,
            'distinct_scales_applied':len(scales),'capped_exchanges':frames[-1]['coupling']['metabolic_supply']['cumulative']['capped_exchanges'],
            'cumulative':frames[-1]['coupling']['metabolic_supply']['cumulative']}
    def test_K4_idempotence(self):
        schedule=lambda step:(1e9,60.,0.,35.)[(step//5)%4]
        a=run(body('energy_fraction',capacity_w=schedule),40)
        b=run(body('energy_fraction',capacity_w=schedule),40)
        self.assertEqual(a,b)
        requested=(20+30)*DT/J_PER_KCAL
        self.assertEqual(assess_exchange(30.,20.,10.,DT,requested_kcal=requested,unmet_kcal=requested/3),
                         assess_exchange(30.,20.,10.,DT,requested_kcal=requested,unmet_kcal=requested/3))
        limiter=SupplyLimiter('energy_fraction',{'a':{'excitation':.2},'b':{'excitation':.2}})
        limiter.coverage=.5;before=limiter.snapshot()
        cap1=limiter.interval();cap2=limiter.interval()
        self.assertEqual(cap1({'a':.8}),cap2({'a':.8}))
        held1=dict(before['held']);held1.update(cap1({'a':.8}))
        held2=dict(before['held']);held2.update(cap2({'a':.8}))
        self.assertEqual(held1,held2)
        assessment=assess_exchange(30.,20.,10.,DT,requested_kcal=requested,unmet_kcal=requested/3)
        kwargs=dict(raw_j={'m_j':.6,'h_j':.4,'w_j':.2},pending_j={'m_j':0.,'h_j':0.,'w_j':0.})
        kwargs['pending_j']={k:kwargs['raw_j'][k]-assessment['sent_j'][k] for k in kwargs['raw_j']}
        self.assertEqual(limiter.close(cap1,assessment,**kwargs),limiter.close(cap1,assessment,**kwargs))
        self.assertEqual(limiter.snapshot(),before)
        REPORT['K4_idempotence']={'two_runtimes_40_exchanges_identical':a==b,'close_is_pure':limiter.snapshot()==before}
    def test_K5_controls_can_fail(self):
        # (a) the ledger check raises on a 1 microjoule leak.
        limiter=SupplyLimiter('energy_fraction',{'a':{'excitation':.2}})
        assessment=assess_exchange(30.,20.,10.,DT,requested_kcal=1.,unmet_kcal=0.)
        raw={'m_j':.6,'h_j':.4,'w_j':.2};pending={k:raw[k]-assessment['sent_j'][k] for k in raw}
        limiter.close(limiter.interval(),assessment,raw_j=raw,pending_j=pending)
        with self.assertRaisesRegex(RuntimeError,'does not close'):
            limiter.close(limiter.interval(),assessment,raw_j=raw,pending_j={**pending,'m_j':pending['m_j']+1e-6})
        # (b) dropping the unsupplied channel is caught.
        short=assess_exchange(30.,20.,10.,DT,requested_kcal=1.,unmet_kcal=.3*DT/J_PER_KCAL)
        dropped=deepcopy(short);dropped['unsupplied_j']={'m_j':0.,'h_j':0.,'w_j':0.}
        with self.assertRaisesRegex(RuntimeError,'does not close'):
            limiter.close(limiter.interval(),dropped,raw_j=raw,pending_j=pending)
        # (c) K1's bit-identity comparison detects a one-part-in-a-million cap.
        plain=body();capped=body('energy_fraction');run(plain,5)
        capped.metabolic_supply.coverage=1-1e-6
        run(capped,5)
        self.assertNotEqual(plain.plant.commands,capped.plant.commands)
        # (d) a cap that scales only COMMANDED muscles fails K2: c stays excited.
        plant=Plant();plant.advance(DT,actuation={k:0.*v for k,v in {'a':.5,'b':.5}.items()})   # naive
        self.assertNotEqual(plant.held,{'a':0.,'b':0.,'c':0.})
        full=IntervalCap(0.,{'a':.2,'b':.2,'c':.2},())
        plant=Plant();plant.advance(DT,actuation=full({'a':.5,'b':.5}))
        self.assertEqual(plant.held,{'a':0.,'b':0.,'c':0.})
        REPORT['K5_controls']={'ledger_leak_1uJ':'raised','dropped_unsupplied':'raised',
            'cap_1e-6_detected':plain.plant.commands!=capped.plant.commands,'commanded_only_cap_leaves_c':.2}
    def test_K6_none_still_aborts_with_original_message(self):
        runtime=body(capacity_w=lambda step:0.)
        with self.assertRaisesRegex(RuntimeError,'^Native muscle energy demand is unmet; mechanical supply feedback is not yet supported$'):
            runtime.step({})
        self.assertTrue(runtime.failed)
    def test_capped_plant_proxy_caps_feedback_path_calls(self):
        cap=IntervalCap(.5,{'a':.2,'b':.4},())
        plant=CappedActuationPlant(Plant({'a':.2,'b':.4}),cap)
        plant.advance(DT,forces=(),actuation={'a':1.})
        self.assertEqual(plant.held,{'a':.5,'b':.2})
        with self.assertRaises(TypeError):plant.advance(DT,(),{'a':1.})
        self.assertFalse(hasattr(plant,'advance_observation'))


class FactoryTests(unittest.TestCase):
    """from_workspace routing with the regional factory's fake owners."""
    def build(self,**options):
        import tempfile
        from unittest.mock import patch
        from scripts.verify_regional_embodied_factory import FactoryTests as Fixture
        class MusclePlant:
            muscle_catalog={}
            def __init__(self,root,output,**kwargs):
                (output/'native').mkdir(parents=True);(output/'native/execution.json').write_text('{}')
            def snapshot(self):return {'time_s':0,'entities':{},'muscles':{'a':{'excitation':.01}},'metabolic_reference':{'M0_w':1,'H0_w':1,'W0_w':0}}
            def close(self):pass
        fixture=Fixture();calls=[]
        with tempfile.TemporaryDirectory() as temporary:
            root,manifests=fixture.fixture(temporary)
            with fixture.patches(root,manifests,calls),patch('ihm.assembly.selective_projection.SelectiveProjectionPlant',MusclePlant):
                body=EmbodiedRuntime.from_workspace(root,root/'output',source_pin=None,**options)
            receipt=json.loads((root/'output/manifest.json').read_text())
            body.close()
        return body,receipt
    def test_factory_records_mode_and_source(self):
        body,receipt=self.build(metabolic_supply='energy_fraction')
        self.assertEqual(receipt['metabolic_supply']['mode'],'energy_fraction')
        self.assertIn('ihm/assembly/metabolic_supply.py',receipt['sources'])
        self.assertEqual(body.metabolic_supply.snapshot()['held'],{'a':.01})
    def test_factory_none_adds_nothing(self):
        body,receipt=self.build()
        self.assertNotIn('metabolic_supply',receipt);self.assertIsNone(body.metabolic_supply)
        self.assertNotIn('ihm/assembly/metabolic_supply.py',receipt['sources'])
    def test_factory_rejects_unknown_mode_before_output(self):
        import tempfile
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ValueError,'Unknown metabolic supply mode'):
                EmbodiedRuntime.from_workspace(Path(temporary),Path(temporary)/'output',metabolic_supply='glycogen')
            self.assertFalse((Path(temporary)/'output').exists())


def native_check(out):
    """Real native mechanical plant: does the cap's command reach the engine?"""
    from ihm.native.mechanical_stream import NativeMechanicalStream
    from ihm.body_constants import MECHANICAL_TARGET_MASS_KG
    stream=NativeMechanicalStream(ROOT,out,environment='supine',target_mass_kg=MECHANICAL_TARGET_MASS_KG)
    try:
        start=stream.snapshot();muscles=sorted(start['muscles'])
        held={m:start['muscles'][m]['excitation'] for m in muscles}
        commanded={m:.3 for m in muscles[:len(muscles)//2]}   # half the muscles; the rest are held
        checkpoint=stream.checkpoint()
        def advance(actuation,steps=5):
            stream.restore(checkpoint)
            for _ in range(steps):state=stream.advance(DT,actuation=actuation)
            return state
        def fields(state):
            return {'time_s':state['time_s'],'muscle_metabolic_energy_j':state['muscle_metabolic_energy_j'],
                    'signed_active_fiber_work_j':state['signed_active_fiber_work_j'],'muscle_heat_energy_j':state['muscle_heat_energy_j'],
                    'excitation':{m:state['muscles'][m]['excitation'] for m in muscles},
                    'bodies':{k:v['transform_ground'] for k,v in state['bodies'].items()}}
        direct=fields(advance(commanded))
        unit_cap=IntervalCap(1.,held,())
        unit=fields(advance(unit_cap(commanded)))
        zero_direct=fields(advance({m:0. for m in muscles}))
        zero_cap=IntervalCap(0.,held,())
        zero=fields(advance(zero_cap(commanded)))
        half_cap=IntervalCap(.5,held,())
        half=fields(advance(half_cap(commanded)))
        expected_half={m:.5*(commanded.get(m,held[m])) for m in muscles}
        result={'muscles':len(muscles),'commanded':len(commanded),
            'scale1_bit_identical_to_uncapped':unit==direct,
            'scale0_bit_identical_to_all_zero_command':zero==zero_direct,
            'scale0_every_engine_excitation_zero':all(v==0. for v in zero['excitation'].values()),
            'scale0_5_engine_excitation_matches':all(math.isclose(half['excitation'][m],expected_half[m],rel_tol=0,abs_tol=1e-15) for m in muscles),
            'held_default_excitation_range':[min(held.values()),max(held.values())],
            'uncommanded_held_under_naive_scale0':sorted({held[m] for m in muscles if m not in commanded})[:3]}
        return result
    finally:
        stream.close()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--native',action='store_true')
    parser.add_argument('--output',type=Path,default=ROOT/'data/derived/metabolic-supply-verification')
    args=parser.parse_args()
    suite=unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__])
    outcome=unittest.TextTestRunner(verbosity=2).run(suite)
    report={'passed':outcome.wasSuccessful(),'tests':outcome.testsRun,
            'failures':[str(t[0]) for t in outcome.failures+outcome.errors],'known_answers':REPORT}
    if args.native:
        import tempfile
        args.output.mkdir(parents=True,exist_ok=True)
        report['native']=native_check(Path(tempfile.mkdtemp(prefix='plant-',dir=args.output))/'plant')
        report['passed']=report['passed'] and all(report['native'][k] for k in report['native'] if k.startswith('scale'))
    text=json.dumps(report,indent=1,default=lambda o:repr(o))
    print(text)
    if args.native:(args.output/'report.json').write_text(text+'\n')
    return 0 if report['passed'] else 1


if __name__=='__main__':
    sys.exit(main())
