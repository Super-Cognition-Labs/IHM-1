"""Metabolic supply as a limit on delivered excitation — an accounting constraint.

What "unmet" is (read from the native source, not inferred):
  ``coupling.muscle_unmet_kcal`` is written by the signed-muscle hook in
  ``data/runtime/physiology/variants/whole_body_integrity_signed_muscle_v2/Tissue.cpp``
  (the ``ihm_signed::active()->unmet_muscle_kcal = ...`` line after the muscle's
  anaerobic glycogen branch). It is the MUSCLE tissue's ``tissueNeededEnergy_kcal``
  left over, in that one 20 ms native interval, after BioGears'
  ``Tissue::CalculateMetabolicConsumptionAndProduction`` has drawn, in order:
  amino acids (hormone-scaled obligatory rate), intracellular TAG (hormone-scaled
  rate-limited), aerobic intracellular glucose, aerobic muscle glycogen (all three
  O2-limited), anaerobic intracellular glucose (2 ATP), anaerobic muscle glycogen
  (3 ATP). The mandatory anaerobic fraction (2.8% of the muscle's basal share) is
  added back into the residual when glycogen cannot cover it. The request it is a
  residual of is ``coupling.muscle_requested_kcal`` = native basal muscle share +
  the signed mechanical increment. Because anaerobic glycogen covers any residual
  while glycogen lasts, a positive value in practice means the intracellular
  glucose AND the muscle glycogen scalar ran out inside that interval.

The native request is one scalar for one whole-body muscle compartment. There is
no per-muscle supply, no ATP/PCr/Pi store, and no intracellular H+ — so no
published fatigue law found maps onto what the engine exposes (see
``docs/METABOLIC_SUPPLY.md`` for the candidates and why each fails). What is
implemented is the fallback declared in that doc: deliver the fraction of the
demand that supply covered. It is NOT a physiological fatigue model.

The law, exactly:
  After exchange k the native reports unmet U_k (kcal). The part of it attributable
  to the mechanical increment is ``u_k = min(U_k * J_PER_KCAL, max(sent_k, 0))``
  where ``sent_k`` is the chemical increment sent in that exchange (J). The
  residual is linear in the request once every branch is exhausted, so this is the
  marginal attribution; the rest, ``U_k*J_PER_KCAL - u_k``, is the native's own
  basal muscle deficit and is reported, not charged to mechanics. The covered
  fraction is ``c_k = (sent_k - u_k) / sent_k`` (1 when u_k is 0 or U_k is within
  the old gate's 1e-12 kcal). Interval k+1 delivers ``c_k`` times the requested
  excitation to EVERY muscle; ``c_k == 1`` delivers the request untouched.

What it does not do, stated so nobody reads it as more:
  * It is one interval late. The native engine has no restorable preview
    (``docs/research/MUSCLE_SUPPLY_FEASIBILITY.md``), so the interval in which
    supply first falls short has already been integrated when that is learned.
  * It scales EXCITATION, and Umberger demand is not linear in excitation and
    activation has memory; the native also sees demand through the 2 s exchange
    lag. So the energy of the next interval is not guaranteed to fit supply.
  * Whatever does not fit is not dropped: it is carried in an explicit
    ``unsupplied`` ledger channel. Every exchange checks
    raw = drawn + unsupplied + pending (each of chemical, heat, signed work) and
    chemical = heat + work on both the drawn and unsupplied triples. The ledger
    closing is bookkeeping; the body conserves energy physically only while the
    cumulative ``unsupplied`` channel is zero, and the frame says which.
"""
from copy import deepcopy
import math

# BioGears' own unit table: data/raw/physiology/biogears/share/etc/UCEDefs.conf:191
# ("UNIT Energy calorie cal 4.184 J AllPrefixes"), the file the native engine
# converts kcal with; kcal = 1000 cal = 4184 J.
J_PER_KCAL=4184.

# The threshold embodied.py applied before this module existed (its abort gate
# ``if unmet>1e-12``). Below it the demand is treated as met, as it always was.
UNMET_DEADBAND_KCAL=1e-12

MODES=('energy_fraction',)

LAW=('Accounting constraint, not a physiological fatigue law: after each exchange the fraction '
     'of the sent muscle chemical increment that native substrate covered scales the NEXT '
     'interval\'s excitation of every muscle; full coverage delivers the request untouched. '
     'Energy that was not covered is carried in an explicit unsupplied ledger channel.')


def _finite(value,label):
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):
        raise ValueError('Invalid '+label)
    return float(value)


def ledger_closes(m,w,h):
    """chemical = work + heat, to the tolerance measure_resting_metabolic_reference uses."""
    return abs(m-w-h)<=1e-9*(1+abs(m)+abs(w)+abs(h))


def assess_exchange(sent_m_w,sent_h_w,sent_w_w,dt,*,requested_kcal,unmet_kcal):
    """Pure: split one exchange's sent increment into drawn and unsupplied energy."""
    dt=_finite(dt,'exchange interval')
    if dt<=0:raise ValueError('Positive exchange interval required')
    sent_m=_finite(sent_m_w,'sent chemical increment')*dt
    sent_h=_finite(sent_h_w,'sent heat increment')*dt
    sent_w=_finite(sent_w_w,'sent work increment')*dt
    if not ledger_closes(sent_m,sent_w,sent_h):raise ValueError('Sent chemical/heat/work ledger mismatch')
    requested=_finite(requested_kcal,'native muscle requested energy')
    unmet=_finite(unmet_kcal,'native muscle unmet energy')
    if unmet<0:raise ValueError('Negative native unmet muscle energy')
    # The native residual cannot exceed the native request it is a residual of.
    if unmet>requested+1e-9*(1+abs(requested)):raise ValueError('Native unmet exceeds native request')
    unmet_j=unmet*J_PER_KCAL
    unsupplied=min(unmet_j,max(sent_m,0.))
    fraction=unsupplied/sent_m if unsupplied>0 else 0.
    drawn={'m_j':sent_m-unsupplied,'h_j':sent_h*(1.-fraction),'w_j':sent_w*(1.-fraction)}
    short={'m_j':unsupplied,'h_j':sent_h*fraction,'w_j':sent_w*fraction}
    for label,triple in (('drawn',drawn),('unsupplied',short)):
        if not ledger_closes(triple['m_j'],triple['w_j'],triple['h_j']):
            raise ValueError(label+' chemical/heat/work ledger mismatch')
    coverage=1. if unmet<=UNMET_DEADBAND_KCAL or unsupplied<=0 else (sent_m-unsupplied)/sent_m
    if not 0<=coverage<=1:raise ValueError('Supply coverage outside [0,1]')
    return {'interval_s':dt,'sent_j':{'m_j':sent_m,'h_j':sent_h,'w_j':sent_w},
            'drawn_j':drawn,'unsupplied_j':short,
            'native_requested_kcal':requested,'native_unmet_kcal':unmet,
            'native_basal_unmet_j':unmet_j-unsupplied,'coverage':coverage}


class IntervalCap:
    """One interval's excitation cap. Tracks the requested (pre-cap) held excitation.

    The plant holds an uncommanded muscle's last command, so a cap that scaled
    only the muscles named in a command would leave the rest at full excitation.
    When the cap binds, or is recovering from having bound, every muscle is
    commanded: ``scale * requested`` while binding, ``requested`` on recovery.
    """
    def __init__(self,scale,held,dirty):
        self.scale=scale;self.held=dict(held);self.dirty=set(dirty);self.calls=0
    def __call__(self,actuation):
        self.calls+=1
        if actuation is None:
            if self.scale==1. and not self.dirty:return None
            actuation={}
        unknown=set(actuation)-set(self.held)
        if unknown:raise ValueError('Supply cap received excitation for unknown muscles: '+', '.join(sorted(unknown)[:5]))
        self.held.update(actuation)
        if self.scale==1. and not self.dirty:return dict(actuation)
        if self.scale==1.:
            delivered=dict(self.held)
        else:
            delivered={name:self.scale*value for name,value in self.held.items()}
        self.dirty={name for name,value in delivered.items() if value!=self.held[name]}
        return delivered


class CappedActuationPlant:
    """Forward a plant, capping ``actuation`` on its two integrating calls only."""
    def __init__(self,plant,cap):
        object.__setattr__(self,'_plant',plant);object.__setattr__(self,'_cap',cap)
    def __getattr__(self,name):
        value=getattr(self._plant,name)
        if name not in ('advance','advance_observation'):return value
        cap=self._cap
        def capped(dt_s,*args,actuation=None,**kwargs):
            if len(args)>1:raise TypeError('Positional actuation cannot pass the supply cap')
            return value(dt_s,*args,actuation=cap(actuation),**kwargs)
        return capped


class SupplyLimiter:
    def __init__(self,mode,muscles):
        if mode not in MODES:raise ValueError('Unknown metabolic supply mode')
        if not isinstance(muscles,dict) or not muscles:raise ValueError('Plant muscle states required for a supply cap')
        held={}
        for name,row in muscles.items():
            if not isinstance(row,dict) or 'excitation' not in row:
                raise ValueError('Supply cap requires the plant\'s held excitation for every muscle')
            held[name]=_finite(row['excitation'],'held excitation')
        self.mode=mode;self.coverage=1.;self.held=held;self.dirty=frozenset()
        zero={'m_j':0.,'h_j':0.,'w_j':0.}
        self.ledger={'raw_j':dict(zero),'sent_j':dict(zero),'drawn_j':dict(zero),'unsupplied_j':dict(zero),
                     'native_basal_unmet_j':0.,'magnitude_j':0.,'exchanges':0,'capped_exchanges':0}

    def interval(self):
        """The cap for the next interval, built from committed state only."""
        return IntervalCap(self.coverage,self.held,self.dirty)

    def close(self,cap,assessment,*,raw_j,pending_j):
        """Pure: the next committed state, or raise if the ledger does not close."""
        ledger=deepcopy(self.ledger)
        for key in ('m_j','h_j','w_j'):
            ledger['raw_j'][key]+=_finite(raw_j[key],'raw increment energy')
            for channel in ('sent_j','drawn_j','unsupplied_j'):ledger[channel][key]+=assessment[channel][key]
        ledger['native_basal_unmet_j']+=assessment['native_basal_unmet_j']
        ledger['magnitude_j']+=sum(abs(raw_j[k]) for k in ('m_j','h_j','w_j'))
        ledger['exchanges']+=1;ledger['capped_exchanges']+=int(cap.scale!=1.)
        residual={}
        for key in ('m_j','h_j','w_j'):
            pending=_finite(pending_j[key],'pending energy')
            residual[key]=ledger['raw_j'][key]-ledger['drawn_j'][key]-ledger['unsupplied_j'][key]-pending
            if abs(residual[key])>1e-9*(1+ledger['magnitude_j']+abs(pending)):
                raise RuntimeError('Metabolic supply ledger does not close on '+key)
        for channel in ('drawn_j','unsupplied_j'):
            t=ledger[channel]
            if not ledger_closes(t['m_j'],t['w_j'],t['h_j']):raise RuntimeError('Cumulative '+channel+' chemical/heat/work mismatch')
        state={'coverage':assessment['coverage'],'held':dict(cap.held),'dirty':frozenset(cap.dirty),'ledger':ledger}
        record={'mode':self.mode,'law':LAW,'excitation_scale_applied':cap.scale,
                'coverage_next':assessment['coverage'],'cap_calls':cap.calls,
                'muscles_held_below_request':len(cap.dirty),
                'interval':deepcopy(assessment),'cumulative':deepcopy(ledger),'closure_residual_j':residual,
                'physically_conserved':ledger['unsupplied_j']['m_j']==0.,
                'basis':'One interval late (no native preview); scales excitation, so next-interval energy is '
                        'not guaranteed to fit supply; unsupplied energy is carried, never dropped. '
                        'See docs/METABOLIC_SUPPLY.md.'}
        return state,record

    def commit(self,state):
        self.coverage=state['coverage'];self.held=state['held'];self.dirty=state['dirty'];self.ledger=state['ledger']

    def snapshot(self):
        return {'mode':self.mode,'coverage':self.coverage,'held':dict(self.held),'dirty':sorted(self.dirty),'ledger':deepcopy(self.ledger)}
