#!/usr/bin/env python3
"""Known answers for ihm/colorimetry.py, the L*a*b* -> sRGB conversion behind the tissue palette.

  .venv/bin/python scripts/test_colorimetry.py

WHY THIS EXISTS. Every `measured` colour in the realistic palette is a published CIE L*a*b* (or
L*C*h) pushed through this conversion. A conversion that has only ever been checked against
itself -- or against the sRGB primaries, whose L*a*b* is definitional and so cannot catch a
chromatic-adaptation error -- would carry every palette entry wrong in the same direction without
anything turning red. So the reference here is computed by someone else, from published inputs:

  Pascale D. RGB coordinates of the Macbeth ColorChecker. The BabelColor Company, complete update
  1 June 2006. https://www.babelcolor.com/index_htm_files/RGB%20Coordinates%20of%20the%20Macbeth
  %20ColorChecker.pdf  -- read in full on 2026-09-18.

  Table 3 ("ColorChecker 2005"): GretagMacbeth's official L*a*b* (CIE D50, 2 degree) for all 24
  patches to three decimals, and the 16-bit sRGB Pascale derived from them by L*a*b* -> XYZ (D50)
  -> Bradford to D65 -> XYZ-to-sRGB matrix -> transfer. Table 2 gives the same in 8 bits.

This exercises the whole non-D65 path, Bradford included, which the builder's own self-test does
not. Two checks against it:

  A  the 16-bit Table 3, with the IEC 61966-2-1 parameters the palette actually uses. This is the
     test of the PIPELINE -- same inputs, same steps, same answer to within a few counts of 65535.
  B  the 8-bit Table 2: every in-gamut channel within 1 count, and the colour difference reported
     in dE*ab. This is the error bar to attach to any palette hex from conversion alone.

A NEGATIVE, recorded because it is how the test learned what it tests. Pascale's Table 6 prints
the sRGB transfer as gamma 0.42 / transition 0.003 and the matrix to four decimals. The first
version of check A used those printed parameters and FAILED at 204 counts of 65535, with a
systematic sign (ours darker at mid-tones, ~0 at white and black) -- the signature of an exponent
error, not of a matrix or adaptation error. Re-run with the exact IEC exponent 1/2.4 it passes at
6 counts, so his table was computed with 1/2.4 and the 0.42 in Table 6 is a rounded display. The
tolerance (40 counts) was not moved. The printed-parameter variant is kept below as a control that
MUST fail: it shows check A can see a 0.3 % transfer error, so its pass means something.

The cyan patch (18) is out of the sRGB gamut and Pascale clips R to 0; the module must report the
clip rather than hide it.

Plus a null case (Bradford from a white to itself is the identity), a round trip, idempotence
(called twice at the same input), the L*C*h -> a*b* step, and the size of the D65/10-degree ->
D65/2-degree approximation the lip entry depends on, so that it is measured rather than assumed.
"""
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ihm import colorimetry as cm  # noqa: E402

# Pascale 2006 Table 3, "ColorChecker 2005": L*a*b* (CIE D50) and 16-bit sRGB. Patch 0 is the
# illuminant itself. Transcribed from the PDF text layer (pdftotext -layout), commas as decimals.
PASCALE_T3 = [
    (0, 'illuminant', (100.0, 0.0, 0.0), (65535, 65535, 65535)),
    (1, 'dark skin', (37.986, 13.555, 14.059), (29684, 20794, 17311)),
    (2, 'light skin', (65.711, 18.130, 17.810), (51033, 37831, 33071)),
    (3, 'blue sky', (49.927, -4.880, -21.925), (23285, 31447, 40035)),
    (4, 'foliage', (43.139, -13.095, 21.905), (23061, 27664, 16548)),
    (5, 'blue flower', (55.112, 8.844, -25.399), (33299, 32893, 45254)),
    (6, 'bluish green', (70.719, -33.397, -0.199), (23760, 48805, 44209)),
    (7, 'orange', (62.661, 36.067, 57.096), (57637, 31797, 12000)),
    (8, 'purplish blue', (40.020, 10.410, -45.964), (17444, 23445, 43738)),
    (9, 'moderate red', (51.124, 48.239, 16.248), (50970, 21055, 24945)),
    (10, 'purple', (30.325, 22.976, -21.587), (24062, 14904, 27134)),
    (11, 'yellow green', (72.532, -23.709, 57.255), (40800, 48564, 16148)),
    (12, 'orange yellow', (71.941, 19.363, 67.857), (59221, 41533, 10089)),
    (13, 'blue', (28.778, 14.179, -50.297), (9090, 16275, 37805)),
    (14, 'green', (55.261, -38.342, 31.370), (17200, 38272, 19051)),
    (15, 'red', (42.101, 53.378, 28.190), (46236, 12506, 14638)),
    (16, 'yellow', (81.733, 4.039, 79.819), (61244, 50998, 5069)),
    (17, 'magenta', (51.935, 49.986, -14.574), (49611, 21580, 38695)),
    (18, 'cyan', (51.038, -28.631, -28.638), (0, 35002, 43613)),
    (19, 'white 9.5', (96.539, -0.425, 1.186), (62954, 63018, 62371)),
    (20, 'neutral 8', (81.257, -0.638, -0.335), (51492, 51965, 52019)),
    (21, 'neutral 6.5', (66.766, -0.734, -0.504), (41301, 41847, 41958)),
    (22, 'neutral 5', (50.867, -0.153, -0.270), (31014, 31145, 31239)),
    (23, 'neutral 3.5', (35.656, -0.421, -1.231), (21187, 21613, 22046)),
    (24, 'black 2', (20.461, -0.079, -0.973), (12507, 12685, 13032)),
]

# Pascale 2006 Table 2, same chart, the "sRGB" column (derived from L*a*b*, 8-bit). Not the
# "sRGB (GMB)" column, which GretagMacbeth computed from spectra by an unstated route.
PASCALE_T2_SRGB8 = {
    0: (255, 255, 255), 1: (116, 81, 67), 2: (199, 147, 129), 3: (91, 122, 156), 4: (90, 108, 64),
    5: (130, 128, 176), 6: (92, 190, 172), 7: (224, 124, 47), 8: (68, 91, 170), 9: (198, 82, 97),
    10: (94, 58, 106), 11: (159, 189, 63), 12: (230, 162, 39), 13: (35, 63, 147), 14: (67, 149, 74),
    15: (180, 49, 57), 16: (238, 198, 20), 17: (193, 84, 151), 18: (0, 136, 170),
    19: (245, 245, 243), 20: (200, 202, 202), 21: (161, 163, 163), 22: (121, 121, 122),
    23: (82, 84, 86), 24: (49, 49, 51)}

# Pascale 2006 Table 6, sRGB row: his printed XYZ->RGB matrix and transfer parameters.
PASCALE_MATRIX = np.array([[3.2405, -1.5371, -0.4985], [-0.9693, 1.8760, 0.0416], [0.0556, -0.2040, 1.0572]])
PASCALE_TRANSFER = {'offset': 0.055, 'gamma': 0.42, 'transition': 0.003, 'slope': 12.92}

# Lindbloom B. "Chromatic Adaptation", brucelindbloom.com/Eqn_ChromAdapt.html, read 2026-09-18: the
# reference whites (ASTM E308-01, 2 degree) and his precomputed Bradford matrices. Transcribed from
# the page text by a parser, not by hand. D55 is here because the Hosoki 2007 oral-mucosa data the
# palette now uses were measured under illuminant D55.
LINDBLOOM_WHITE = {'A_2': (1.09850, 1.00000, 0.35585), 'C_2': (0.98074, 1.00000, 1.18232),
                   'D50_2': (0.96422, 1.00000, 0.82521), 'D55_2': (0.95682, 1.00000, 0.92149),
                   'D65_2': (0.95047, 1.00000, 1.08883)}
LINDBLOOM_BRADFORD_TO_D65 = {
    'A_2': (0.8446965, -0.1179225, 0.3948108, -0.1366303, 1.1041226, 0.1291718, 0.0798489, -0.1348999, 3.1924009),
    'C_2': (0.9904476, -0.0071683, -0.0116156, -0.0123712, 1.0155950, -0.0029282, -0.0035635, 0.0067697, 0.9181569),
    'D50_2': (0.9555766, -0.0230393, 0.0631636, -0.0282895, 1.0099416, 0.0210077, 0.0122982, -0.0204830, 1.3299098),
    'D55_2': (0.9726856, -0.0135482, 0.0361731, -0.0167463, 1.0049102, 0.0120598, 0.0070026, -0.0116372, 1.1869548)}

# Tolerances, fixed before the comparison was run and not moved after it.
TOL_A_16BIT = 40          # counts of 65535 (0.06 %)
TOL_B_8BIT = 1            # counts of 255, IEC parameters, in-gamut channels
TOL_B_DE = 1.0            # dE*ab (CIE 1976) between our 8-bit colour and his, both read back to Lab
TOL_ROUNDTRIP = 1e-9


def main():
    failures = []

    def check(name, ok, detail):
        print(('PASS ' if ok else 'FAIL ')+name+' :: '+detail)
        if not ok:
            failures.append(name)

    # A. the pipeline, IEC parameters, against the 16-bit table
    def worst16(**kwargs):
        worst = 0
        for no, name, lab, rgb16 in PASCALE_T3:
            got = cm.lab_to_srgb(lab, 'D50_2', **kwargs)
            mine = np.round(got['encoded']*65535).astype(int)
            worst = max(worst, int(np.max(np.abs(mine-np.array(rgb16)))))
        return worst
    worst_a = worst16()
    check('A: Pascale Table 3 (16-bit) reproduced with the IEC parameters', worst_a <= TOL_A_16BIT,
          f'24 patches + white, max |d| = {worst_a} of 65535 (tolerance {TOL_A_16BIT})')
    printed = worst16(matrix=PASCALE_MATRIX, transfer=PASCALE_TRANSFER)
    check('control: the rounded parameters printed in his Table 6 are DETECTED', printed > TOL_A_16BIT,
          f'gamma 0.42 + 4-decimal matrix gives max |d| = {printed} of 65535; check A can fail')

    # B. the conversion the palette uses (IEC), against the 8-bit table
    worst_b, worst_de, clip_reported = 0, 0.0, None
    for no, name, lab, rgb16 in PASCALE_T3:
        got = cm.lab_to_srgb(lab, 'D50_2')
        ref = np.array(PASCALE_T2_SRGB8[no])
        if no == 18:
            clip_reported = got['receipt']['out_of_srgb_gamut'] and got['rgb8'][0] == 0
            channels = [1, 2]
        else:
            channels = [0, 1, 2]
        worst_b = max(worst_b, int(np.max(np.abs(np.array(got['rgb8'])[channels]-ref[channels]))))
        de = cm.delta_e76(cm.srgb_to_lab(got['hex']), cm.srgb_to_lab('#%02x%02x%02x' % tuple(ref)))
        worst_de = max(worst_de, de)
    check('B: IEC sRGB within 1 count of Pascale Table 2 on every in-gamut channel', worst_b <= TOL_B_8BIT,
          f'max |d| = {worst_b} of 255 (tolerance {TOL_B_8BIT})')
    check('B: colour difference to Pascale Table 2', worst_de <= TOL_B_DE,
          f'max dE*ab = {worst_de:.3f} under D65/2 (tolerance {TOL_B_DE})')
    check('out-of-gamut cyan is clipped AND reported', bool(clip_reported),
          'patch 18: R clipped to 0 and out_of_srgb_gamut=True')

    # Known D65 answers: the sRGB primaries and white, whose L*a*b* under D65/2 is definitional.
    primaries = {'#ff0000': (53.2408, 80.0925, 67.2032), '#00ff00': (87.7347, -86.1827, 83.1793),
                 '#0000ff': (32.2970, 79.1875, -107.8602), '#ffffff': (100.0, 0.0, 0.0)}
    got = {h: cm.lab_to_srgb(lab)['hex'] for h, lab in primaries.items()}
    worst = max(cm.delta_e76(cm.srgb_to_lab(h), lab) for h, lab in primaries.items())
    check('sRGB primaries and white under D65/2', all(got[h] == h for h in primaries),
          f'{sum(got[h] == h for h in primaries)}/4 exact; inverse max dE = {worst:.4f}')

    # Second, independent published reference: whites and Bradford matrices (Lindbloom).
    white_err = max(abs(cm.WHITE[k][i]/100.0-v[i]) for k, v in LINDBLOOM_WHITE.items() for i in range(3))
    check('reference whites equal the ASTM E308-01 table', white_err < 1e-9,
          f'{len(LINDBLOOM_WHITE)} whites, max |d| = {white_err:.1e}')
    brad_err = max(float(np.max(np.abs(cm.bradford_matrix(k, 'D65_2')-np.array(v).reshape(3, 3))))
                   for k, v in LINDBLOOM_BRADFORD_TO_D65.items())
    check('Bradford matrices to D65 equal the published ones (A, C, D50, D55)', brad_err < 1e-6,
          f'max element |d| = {brad_err:.1e} (published to 7 decimals)')

    # Null case: adapting a white to itself must do nothing.
    eye = max(float(np.max(np.abs(cm.bradford_matrix(w, w)-np.eye(3)))) for w in cm.WHITE)
    check('null case: Bradford from any white to itself is the identity', eye < 1e-12, f'max |M - I| = {eye:.1e}')
    d50_white = cm.lab_to_srgb((100.0, 0.0, 0.0), 'D50_2')
    check('D50 perfect diffuser adapts to sRGB white', d50_white['hex'] == '#ffffff'
          and d50_white['receipt']['chromatic_adaptation'] == 'bradford_to_D65_2', d50_white['hex'])

    # Round trip in floating point (no 8-bit rounding), through the adaptation, in gamut only.
    rt = 0.0
    for no, name, lab, _ in PASCALE_T3:
        got = cm.lab_to_srgb(lab, 'D50_2')
        if np.any(got['linear_unclipped'] < 0) or np.any(got['linear_unclipped'] > 1):
            continue      # a clipped value has no inverse; the D50 white clips by 3e-7 in linear R
        rt = max(rt, cm.delta_e76(cm.srgb_to_lab(got['encoded'], 'D50_2'), lab))
    check('float round trip Lab(D50) -> sRGB -> Lab(D50)', rt < TOL_ROUNDTRIP, f'max dE = {rt:.1e}')

    # Idempotence: the same input twice gives the same answer.
    a = cm.lab_to_srgb((41.9, 19.565, 13.197), 'D65_10')
    b = cm.lab_to_srgb((41.9, 19.565, 13.197), 'D65_10')
    check('called twice at the same input', a['hex'] == b['hex'] and np.array_equal(a['encoded'], b['encoded']),
          a['hex'])

    # L*C*h -> a*b*: a hue of 0 is pure +a*, 90 is pure +b*; and the lip entry's published numbers.
    ok = (np.allclose(cm.lch_to_lab(50, 20, 0), (50, 20, 0)) and np.allclose(cm.lch_to_lab(50, 20, 90), (50, 0, 20), atol=1e-12)
          and np.allclose(cm.lab_to_lch(cm.lch_to_lab(41.9, 23.6, 34.0)), (41.9, 23.6, 34.0)))
    lip = cm.lch_to_lab(41.9, 23.6, 34.0)
    check('L*C*h -> L*a*b*', ok, f'Vergnaud 2024 L*41.9 C*23.6 h34.0 -> a*{lip[1]:.3f} b*{lip[2]:.3f}')

    # The size of the observer approximation: D65/10 treated by Bradford vs treated as if 2 degree.
    spread = []
    for lab in [(41.9, 19.565, 13.197), (52.9, 23.3, 14.9), (73.5, 2.2, 11.9), (60.0, 10.0, 17.0)]:
        adapted = cm.lab_to_srgb(lab, 'D65_10')
        naive = cm.lab_to_srgb(lab, 'D65_2')
        spread.append(cm.delta_e76(cm.srgb_to_lab(adapted['encoded']), cm.srgb_to_lab(naive['encoded'])))
    check('D65/10 -> D65/2 approximation is small and measured', max(spread) < 1.0,
          f'Bradford(10->2 deg white) vs treating 10 deg as 2 deg: max dE = {max(spread):.3f} over 4 tissue colours')

    print(f'known-answer error: {worst_a}/65535 counts (Pascale 2006 Table 3); '
          f'IEC {worst_b}/255 counts, {worst_de:.3f} dE*ab (Pascale 2006 Table 2)')
    if failures:
        raise SystemExit('FAILED: '+'; '.join(failures))


if __name__ == '__main__':
    main()
