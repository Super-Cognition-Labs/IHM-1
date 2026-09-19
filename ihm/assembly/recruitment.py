"""Size-principle motor-unit pool: net pool input -> excitation for one muscle.

WHAT THIS IS. A deterministic, rested motor-unit pool after Fuglevand, Winter &
Patla (1993), in the exact parameterisation restated by Potvin & Fuglevand (2017):

    Potvin JR, Fuglevand AJ (2017) A motor unit-based model of muscle fatigue.
    PLoS Comput Biol 13(6): e1005581. doi:10.1371/journal.pcbi.1005581 (CC-BY)
    Methods, "Motor unit pool" and "Force-frequency relation", pp. 19-20,
    Eqs 1-7. Retained at data/raw/sensorimotor/potvin_fuglevand_2017.pdf and
    pinned below by sha256. The fatigue terms (Eqs 8 onward) are NOT used.

    Fuglevand AJ, Winter DA, Patla AE (1993) Models of recruitment and rate
    coding organization in motor-unit pools. J Neurophysiol 70:2470-2488. The
    original; not retained here (publisher returns 403). Every number below is
    read from the 2017 restatement, not from memory of the 1993 paper.

The size principle (Henneman 1957) is carried by the ordering: unit i has twitch
force P(i) and recruitment threshold RTE(i), both exponential in i, so weak units
are recruited first and strong units last (Potvin & Fuglevand, Eqs 1-2).

    P(i)   = exp(ln(RP) (i-1)/(n-1))            RP = 100, n = 120        Eq 1
    RTE(i) = exp(ln(RR) (i-1)/(n-1))            RR = 50                  Eq 2
    R(i,E) = g (E - RTE(i)) + minR, <= maxR(i)  g = 1 imp/s/unit, minR = 8  Eq 3
    maxR   = 35 imp/s (unit 1) -> 25 imp/s (unit 120), "decreased uniformly"
    Emax   = RTE(n) + (maxR(n) - minR)/g = 67 excitation units
    NR     = R x CT(i)                                                    Eq 4
    NF     = 0.3 NR                (NR <= 0.4)                            Eq 5
           = 1 - exp(-2 NR^3)      (NR >  0.4)                            Eq 6
    F(i)   = NF x P(i);  muscle force = sum_i F(i)                        Eq 7

ONE PARAMETER IS GIVEN ONLY AS A RANGE. Contraction time is "an inverse function
of twitch amplitude (see [19])", 90 ms for unit 1 to 30 ms for unit 120. The law
used here is the power law CT(i) = 90 ms x P(i)^(-ln 3 / ln 100). It is checked,
not assumed: it reproduces the five labelled points of Potvin & Fuglevand Fig 1B
(MU20/40/60/80/100 read at ~76/63/52/43.5/36 ms) to <1 ms, and a law linear in
P misses MU60 by 33 ms (scripts/verify_recruitment.py, known-answer test).

WHAT THIS PRODUCES. The pool's rate-averaged isometric force as a fraction of its
own maximum, F(E)/F(Emax), in [0, 1]. That fraction is handed to native mechanics
AS AN EXCITATION. Native mechanics owns activation dynamics, fibre dynamics and
force; this module never produces activation or force for the plant. Twitch time
courses are not simulated (Eq 7 is a steady-state force-frequency mapping), and
nothing here is stochastic: there is no generator to draw, by construction.

WHAT IS NOT MUSCLE-SPECIFIC. The same pool shape (n, RP, RR, rates, contraction
times) is applied to every muscle, because the source gives one representative
pool, "not definitive characterizations of any specific skeletal muscle" (p. 19).
`max_isometric_force_n` sets only the absolute force of each unit (the pool at Emax
sums to the muscle's Fmax); it does not change the normalised curve, and ordering
units ACROSS muscles by Fmax is not a size-principle statement any source here
makes, so it is not done. Per-muscle unit counts and recruitment ranges are
unmeasured here; motor-unit counts and upper recruitment limits per muscle would
measure them.
"""
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import numpy as np

POTVIN_FUGLEVAND_URL = 'https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1005581&type=printable'
POTVIN_FUGLEVAND_DOI = '10.1371/journal.pcbi.1005581'
POTVIN_FUGLEVAND_PATH = 'data/raw/sensorimotor/potvin_fuglevand_2017.pdf'
POTVIN_FUGLEVAND_SHA256 = 'baed9b0b0ff5deea9e463839b1672e722cda7c427074eed7334f5cbf6c2f2ca6'
FUGLEVAND_1993_DOI = '10.1152/jn.1993.70.6.2470'


@dataclass(frozen=True)
class PoolParameters:
    """Potvin & Fuglevand 2017, Methods pp. 19-20. Defaults ARE the source values."""
    units: int = 120                       # n
    twitch_force_range: float = 100.       # RP, Eq 1
    recruitment_range: float = 50.         # RR, Eq 2
    minimum_rate_hz: float = 8.            # minR
    rate_gain_hz_per_unit: float = 1.      # g
    peak_rate_first_hz: float = 35.        # maxR(1)
    peak_rate_last_hz: float = 25.         # maxR(n)
    contraction_time_first_s: float = .090 # CT(1)
    contraction_time_last_s: float = .030  # CT(n)
    inversion_grid: int = 20001            # numerical inverse only; not a source value

    def __post_init__(self):
        if type(self.units) is not int or not 2 <= self.units <= 10000:
            raise ValueError('units must be an integer in [2, 10000]')
        if type(self.inversion_grid) is not int or not 1001 <= self.inversion_grid <= 10**6:
            raise ValueError('inversion_grid must be an integer in [1001, 1e6]')
        for key, value in asdict(self).items():
            if key in ('units', 'inversion_grid'):
                continue
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError(f'{key} must be finite and positive')
        if self.twitch_force_range <= 1 or self.recruitment_range <= 1:
            raise ValueError('Ranges must exceed 1')
        if not self.peak_rate_last_hz <= self.peak_rate_first_hz or self.peak_rate_last_hz <= self.minimum_rate_hz:
            raise ValueError('Peak rates must be onion-skin ordered and above the minimum rate')
        if not self.contraction_time_last_s < self.contraction_time_first_s:
            raise ValueError('Stronger units must be faster')


def normalized_force(normalized_rate):
    """Potvin & Fuglevand Eqs 5-6. Continuous at NR = 0.4 (0.12 vs 0.1201)."""
    nr = np.asarray(normalized_rate, float)
    return np.where(nr <= .4, .3 * nr, 1. - np.exp(-2. * nr ** 3))


class MotorUnitPool:
    """One representative pool. Pure: every method is a function of its arguments."""

    def __init__(self, parameters=None):
        self.parameters = parameters or PoolParameters()
        if not isinstance(self.parameters, PoolParameters):
            raise ValueError('Expected PoolParameters')
        p = self.parameters
        fraction = np.arange(p.units, dtype=float) / (p.units - 1)
        self.twitch_force = np.exp(math.log(p.twitch_force_range) * fraction)          # Eq 1
        self.threshold = np.exp(math.log(p.recruitment_range) * fraction)             # Eq 2
        self.peak_rate_hz = p.peak_rate_first_hz + (p.peak_rate_last_hz - p.peak_rate_first_hz) * fraction
        exponent = -math.log(p.contraction_time_first_s / p.contraction_time_last_s) / math.log(p.twitch_force_range)
        self.contraction_time_s = p.contraction_time_first_s * self.twitch_force ** exponent
        self.max_excitation = float(self.threshold[-1] + (self.peak_rate_hz[-1] - p.minimum_rate_hz) / p.rate_gain_hz_per_unit)
        self.max_force = float(self._force(self.max_excitation))
        # Monotone grid for the inverse. F is nondecreasing with small upward
        # jumps at each recruitment (a unit starts at minR, not at 0 imp/s).
        self._grid_excitation = np.linspace(0., self.max_excitation, p.inversion_grid)
        self._grid_output = self.output(self._grid_excitation / self.max_excitation)
        if np.any(np.diff(self._grid_output) < -1e-15):
            raise AssertionError('Pool output must be nondecreasing in excitation')
        identity = {'parameters': asdict(p), 'source_doi': POTVIN_FUGLEVAND_DOI,
                    'source_sha256': POTVIN_FUGLEVAND_SHA256, 'ct_law': 'power law in P, checked against Fig 1B'}
        self.model_sha256 = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()

    def rates_hz(self, excitation):
        """Eq 3 with saturation; 0 below threshold. `excitation` in source units."""
        e = np.asarray(excitation, float)[..., None]
        p = self.parameters
        rate = np.minimum(p.rate_gain_hz_per_unit * (e - self.threshold) + p.minimum_rate_hz, self.peak_rate_hz)
        return np.where(e >= self.threshold, rate, 0.)

    def unit_forces(self, excitation):
        """Eq 7 per unit, in units of P(1)."""
        return normalized_force(self.rates_hz(excitation) * self.contraction_time_s) * self.twitch_force

    def _force(self, excitation):
        return self.unit_forces(excitation).sum(axis=-1)

    def recruited(self, input_fraction):
        e = np.clip(np.asarray(input_fraction, float), 0., 1.) * self.max_excitation
        return (e[..., None] >= self.threshold).sum(axis=-1)

    def output(self, input_fraction):
        """Pool input as a fraction of Emax -> force fraction of F(Emax), in [0, 1]."""
        x = np.asarray(input_fraction, float)
        if not np.isfinite(x).all():
            raise ValueError('Pool input must be finite')
        return self._force(np.clip(x, 0., 1.) * self.max_excitation) / self.max_force

    def input_for(self, output_fraction):
        """Smallest grid input whose output reaches `output_fraction` (so output(input_for(y)) >= y)."""
        y = np.asarray(output_fraction, float)
        if not np.isfinite(y).all():
            raise ValueError('Pool target must be finite')
        index = np.searchsorted(self._grid_output, np.clip(y, 0., 1.), side='left')
        index = np.minimum(index, len(self._grid_excitation) - 1)
        return self._grid_excitation[index] / self.max_excitation

    def unit_forces_n(self, input_fraction, max_isometric_force_n):
        """Absolute unit forces: the pool at Emax sums to the muscle's Fmax."""
        f = float(max_isometric_force_n)
        if not math.isfinite(f) or f <= 0:
            raise ValueError('max_isometric_force_n must be finite and positive')
        e = np.clip(np.asarray(input_fraction, float), 0., 1.) * self.max_excitation
        return self.unit_forces(e) * (f / self.max_force)

    def combine(self, segmental_excitation, descending_input_fraction):
        """Sum segmental and descending input AT THE POOL, then recruit.

        `segmental_excitation` is a source reflex stimulation, which Geyer & Herr
        define at the level of muscle stimulation, i.e. already a pool OUTPUT. It is
        mapped back to the pool input that would produce it on its own; descending
        input is added there; the excitation is the source stimulation plus what
        the pool adds when the descending drive is summed on top. At zero descending
        input this returns the source stimulation exactly -- the reflex law is not
        re-shaped by the pool -- and descending drive is recruited by size.
        """
        s = np.clip(np.asarray(segmental_excitation, float), 0., 1.)
        d = np.asarray(descending_input_fraction, float)
        if not np.isfinite(d).all() or np.any(d < 0) or np.any(d > 1):
            raise ValueError('Descending pool input must be a fraction in [0, 1]')
        base = self.input_for(s)
        total = np.clip(base + d, 0., 1.)
        added = self.output(total) - self.output(base)
        return np.clip(s + added, 0., 1.), base, total
