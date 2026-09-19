"""CIE colorimetry to sRGB, written once and tested against published known answers.

Used by scripts/build_tissue_colour_palettes.py to turn published tissue colorimetry into palette
colours. Every conversion a palette entry depends on lives here, not as arithmetic in a comment,
and scripts/test_colorimetry.py checks it against a published reference whose answers were
computed independently of this code (Pascale 2006, BabelColor, "RGB coordinates of the Macbeth
ColorChecker", Tables 2 and 3: ColorChecker 2005 L*a*b* under D50 -> sRGB via Bradford).

The pipeline, and what each step assumes

  L*C*h  -> L*a*b*   a* = C* cos h, b* = C* sin h, h in degrees. Exact.
  L*a*b* -> XYZ      CIE 15:2004 inverse, with the CIE-exact constants eps = 216/24389 and
                     kappa = 24389/27, relative to the reference white of the STATED illuminant and
                     observer (WHITE below).
  XYZ    -> XYZ(D65) von Kries-type chromatic adaptation in the Bradford cone space, applied only
                     when the source white is not D65/2 degree. The simplified (linear) Bradford
                     transform, as used in the development of sRGB and in Pascale 2006, eq. 4-7.
  XYZ    -> linear   IEC 61966-2-1 sRGB matrix (D65 white, ITU-R BT.709 primaries).
  linear -> sRGB'    IEC 61966-2-1 piecewise transfer: 12.92 x below 0.0031308, else
                     1.055 x^(1/2.4) - 0.055. Out-of-gamut linear values are clipped to [0, 1] and the
                     clip is REPORTED, never silent.

An observer change is not a chromatic adaptation. A measurement reported under D65 with the 1964
10-degree observer is normalised by the 10-degree D65 white and then carried to the 2-degree D65
white with the same Bradford step. That is an approximation -- converting between observers
exactly needs the reflectance spectrum -- and the receipt says `bradford_to_D65_2` so it is
visible on every entry that used it. Its size is measured in scripts/test_colorimetry.py.
"""
from __future__ import annotations

import math

import numpy as np

# Reference whites, tristimulus x100, Y = 100.
#   2-degree A, C, D50, D55, D65: ASTM E308-01, as tabulated by B. Lindbloom
#     (brucelindbloom.com/Eqn_ChromAdapt.html, read 2026-09-18); C, D50 and D65 also appear in
#     Pascale 2006 Table 7. scripts/test_colorimetry.py checks all five against that table and checks
#     the Bradford matrices built from them against Lindbloom's published matrices.
#   10-degree D65, D50, C: carried from the first palette build (ASTM E308 10-degree values). NOT
#     re-verified against a table read in this build; only D65_10 is used (lip vermilion), and the
#     test measures how much the 10 -> 2 degree treatment moves a tissue colour.
WHITE = {'D65_2': (95.047, 100.0, 108.883), 'D50_2': (96.422, 100.0, 82.521),
         'D55_2': (95.682, 100.0, 92.149),
         'C_2': (98.074, 100.0, 118.232), 'A_2': (109.850, 100.0, 35.585),
         'D65_10': (94.811, 100.0, 107.304), 'D50_10': (96.720, 100.0, 81.427),
         'C_10': (97.285, 100.0, 116.145)}

# Bradford cone-response matrix (Lam 1985; Pascale 2006 eq. 4).
BRADFORD = np.array([[0.8951, 0.2664, -0.1614],
                     [-0.7502, 1.7135, 0.0367],
                     [0.0389, -0.0685, 1.0296]])

# IEC 61966-2-1 sRGB, linear XYZ (Y of white = 1) -> linear RGB. Seven-figure form.
M_XYZ_TO_SRGB = np.array([[3.2404542, -1.5371385, -0.4985314],
                          [-0.9692660, 1.8760108, 0.0415560],
                          [0.0556434, -0.2040259, 1.0572252]])

EPSILON = 216/24389
KAPPA = 24389/27

# The IEC transfer function. Pascale 2006 Table 6 prints a rounded form (gamma 0.42, transition
# 0.003); the test reproduces his table with his parameters as well, so the two can be told apart.
IEC_TRANSFER = {'offset': 0.055, 'gamma': 1/2.4, 'transition': 0.0031308, 'slope': 12.92}


def lch_to_lab(L, C, h_degrees):
    """CIE L*C*h(ab) -> L*a*b*. Exact."""
    h = math.radians(h_degrees)
    return (float(L), float(C*math.cos(h)), float(C*math.sin(h)))


def lab_to_lch(lab):
    L, a, b = lab
    return (float(L), float(math.hypot(a, b)), float(math.degrees(math.atan2(b, a)) % 360.0))


def _white(illuminant):
    if illuminant not in WHITE:
        raise ValueError('Unknown illuminant/observer: '+str(illuminant)+'; known: '+', '.join(sorted(WHITE)))
    return np.array(WHITE[illuminant], dtype=float)/100.0


def lab_to_xyz(lab, illuminant='D65_2'):
    """CIE L*a*b* relative to the named white -> XYZ with Y of the white = 1."""
    L, a, b = (float(v) for v in lab)
    fy = (L+16)/116
    fx, fz = fy+a/500, fy-b/200
    def finv(t):
        return t**3 if t**3 > EPSILON else (116*t-16)/KAPPA
    yr = ((L+16)/116)**3 if L > KAPPA*EPSILON else L/KAPPA
    return np.array([finv(fx), yr, finv(fz)])*_white(illuminant)


def xyz_to_lab(xyz, illuminant='D65_2'):
    r = np.asarray(xyz, dtype=float)/_white(illuminant)
    f = np.where(r > EPSILON, np.cbrt(r), (KAPPA*r+16)/116)
    return (float(116*f[1]-16), float(500*(f[0]-f[1])), float(200*(f[1]-f[2])))


def bradford_matrix(source, destination):
    """The 3x3 Bradford adaptation matrix from one named white to another (Pascale 2006 eq. 6)."""
    s = BRADFORD@_white(source)
    d = BRADFORD@_white(destination)
    return np.linalg.inv(BRADFORD)@np.diag(d/s)@BRADFORD


def encode(linear, transfer=IEC_TRANSFER):
    """Linear RGB in [0, 1] -> non-linear R'G'B' in [0, 1]."""
    x = np.asarray(linear, dtype=float)
    return np.where(x < transfer['transition'], transfer['slope']*x,
                    (1+transfer['offset'])*np.power(np.maximum(x, 0.0), transfer['gamma'])-transfer['offset'])


def decode(encoded):
    """IEC 61966-2-1 inverse transfer: R'G'B' in [0, 1] -> linear."""
    v = np.asarray(encoded, dtype=float)
    return np.where(v <= 0.04045, v/12.92, ((v+0.055)/1.055)**2.4)


def lab_to_srgb(lab, illuminant='D65_2', matrix=M_XYZ_TO_SRGB, transfer=IEC_TRANSFER):
    """CIE L*a*b* under the stated illuminant/observer -> sRGB, with a receipt.

    Returns a dict: `encoded` (float R'G'B' in [0, 1] after clipping), `rgb8` (ints), `hex`,
    `linear_unclipped`, and the receipt fields the palette records.
    """
    xyz = lab_to_xyz(lab, illuminant)
    adapted = illuminant != 'D65_2'
    if adapted:
        xyz = bradford_matrix(illuminant, 'D65_2')@xyz
    linear = np.asarray(matrix)@xyz
    out_of_gamut = bool(np.any(linear < -1e-6) or np.any(linear > 1+1e-6))
    excursion = float(np.max(np.abs(np.clip(linear, 0, 1)-linear)))
    encoded = encode(np.clip(linear, 0, 1), transfer)
    rgb8 = np.clip(np.round(encoded*255), 0, 255).astype(int)
    return {'encoded': encoded, 'rgb8': tuple(int(v) for v in rgb8),
            'hex': '#%02x%02x%02x' % tuple(int(v) for v in rgb8),
            'linear_unclipped': linear,
            'receipt': {'illuminant_observer': illuminant,
                        'chromatic_adaptation': 'bradford_to_D65_2' if adapted else 'none',
                        'transfer': 'IEC 61966-2-1 sRGB',
                        'out_of_srgb_gamut': out_of_gamut,
                        'gamut_clip_linear_excursion': round(max(0.0, excursion), 6)}}


def srgb_to_lab(rgb, illuminant='D65_2'):
    """sRGB as a '#rrggbb' hex string, or as float R'G'B' in [0, 1] -> L*a*b* under `illuminant`.

    The inverse of lab_to_srgb for in-gamut colours. For a non-D65/2 white the Bradford step is
    inverted, so srgb_to_lab(lab_to_srgb(lab, 'D50_2')['encoded'], 'D50_2') returns lab.
    """
    if isinstance(rgb, str):
        v = np.array([int(rgb[i:i+2], 16)/255 for i in (1, 3, 5)])
    else:
        v = np.asarray(rgb, dtype=float)
        if v.min() < 0.0 or v.max() > 1.0:
            raise ValueError("float R'G'B' must lie in [0, 1]; pass 8-bit values as a hex string")
    xyz = np.linalg.inv(M_XYZ_TO_SRGB)@decode(v)
    if illuminant != 'D65_2':
        xyz = bradford_matrix('D65_2', illuminant)@xyz
    return xyz_to_lab(xyz, illuminant)


def delta_e76(lab1, lab2):
    return float(np.linalg.norm(np.asarray(lab1, dtype=float)-np.asarray(lab2, dtype=float)))

