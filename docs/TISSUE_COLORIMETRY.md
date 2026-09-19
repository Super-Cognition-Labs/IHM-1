# Tissue colorimetry: what the `realistic` palette rests on

The `realistic` palette claims "fresh, in vivo tissue under neutral white light". This doc records, for
each of its 74 tissue-class roles, what published colorimetry exists, what was used, and what was
rejected and why. It covers the literature search of 18 Sep 2026. The didactic palette (the app
default, 17 of 17 roles synthesized by design) is out of scope.

- **Code:** `scripts/tissue_colour_tables.py` (data), `scripts/build_tissue_colour_palettes.py`
  (builder), `ihm/colorimetry.py` (the conversion), `scripts/test_colorimetry.py` (its known answers).
- **Output:** `data/derived/tissue-colour-palette-candidate-v1/` (gitignored). Regenerate it with
  `.venv/bin/python scripts/build_tissue_colour_palettes.py --self-test`.
  `ihm.app.experiments.read_palette` fails closed on any hash mismatch.

## Result

| | measured | transferred | synthesized |
|---|---|---|---|
| before (7 Sep build) | 5 | 4 | 65 |
| **after (18 Sep)** | **5** | **5** | **64** |

The headline moved by one role. The rest of the change is in what the numbers stand on:

- **One promotion.** `mucosa_pharyngeal` went from synthesized to transferred. Its source is
  instrument colorimetry of buccal lining mucosa under a stated illuminant.
- **Three re-sourced roles, tier unchanged.**
  - `skin_palmoplantar` now comes from a source that states both illuminant and observer.
  - `mucosa_palatal` now borrows from the right mucosa class (lining, not masticatory).
  - `mucosa_lingual` no longer depends on a composite that combined a camera L*/a* with a b* taken
    from the gingiva.
- **Two measured roles fall below the bar this campaign applied.** `mucosa_gingival` (Ho 2015) and
  `nail_plate` (Horibata 2025) come from sources that state neither illuminant nor observer. Both full
  texts were re-read on 18 Sep to check. The previous build accepted them with D65/2° assumed, and
  that assumption is carried on the record. They were left as `measured` rather than silently demoted
  or silently kept. Every converted entry now carries `illuminant_observer_basis`, and the build
  manifest lists the measured roles by basis:
  - **stated:** `mucosa_lip_vermilion` (D65/10°), `skin_palmoplantar` (D65/10°), `tooth` (D65/2°)
  - **assumed:** `mucosa_gingival`, `nail_plate`

  **Whether an assumed-illuminant measurement may keep the `measured` tier is a decision for the
  owner, not for this campaign.** If the answer is no, measured drops to 3.

**No published colorimetry exists for any visceral, glandular, neural, vascular, connective or ocular
tissue that meets the bar.** Those are 56 of the 64 synthesized roles. The other eight are nasal,
gastric, anal-margin and tonsil mucosa, the two skin roles, hair, and the field overlay, which is not a
tissue. The near misses are listed per role below, because the next campaign should start from them.

## The bar

- **measured.** A measurement of this tissue, in vivo or fresh, in humans. The source must state both
  illuminant and observer, and the value is converted by the code below. Only numbers read in the
  source during this campaign count. An abstract's number that was not seen in the paper does not.
- **transferred.** A measurement of a different but closely related tissue, species or site, with the
  reason it transfers. An unstated observer was accepted here with 2° assumed. That assumption is
  written on the entry, and its size is measured: 0.125 ΔE\*ab at D65 over four tissue colours
  (`scripts/test_colorimetry.py`). An unstated illuminant was not accepted.
- **Rejected everywhere.**
  - Camera-derived L\*a\*b\* whose white point is not stated. It has no defined conversion to sRGB.
  - Hyperspectral data starting at 500 nm (e.g. TIVITA). It cannot produce a colour.
  - Spectra that are only plotted.
  - Anything measured on a print or a processed carcass the authors say was altered.

## The conversion, and its known-answer error

`ihm/colorimetry.py` runs these steps in order:

1. L\*C\*h → a\*b\*, where the source publishes L\*C\*h.
2. L\*a\*b\* → XYZ, relative to the white of the stated illuminant and observer. The constants are the
   CIE-exact ε and κ.
3. A Bradford adaptation to D65/2°, when the source white is not D65/2°.
4. The IEC 61966-2-1 matrix and transfer function.
5. A gamut clip. The clip is reported in the palette record, never applied silently.

No arithmetic is done in the tables. Entries hold the published numbers verbatim (`lab` or `lch`) and
the published illuminant.

`scripts/test_colorimetry.py` checks the conversion against two independent published references.

**Pascale 2006 (BabelColor, "RGB coordinates of the Macbeth ColorChecker"), read in full.** The test
converts the GretagMacbeth 2005 L\*a\*b\* (D50) of all 24 patches, plus white, to sRGB through Bradford:

- **Table 3 (16-bit):** maximum error **6 of 65535 counts**. Tolerance 40, fixed before the run.
- **Table 2 (8-bit):** maximum **1 count of 255** on every in-gamut channel, **0.50 ΔE\*ab**. Both
  sides are quantised to 8 bits.
- The cyan patch is out of gamut. The test checks that it is clipped to R=0 and that the clip is
  reported.
- **Negative, recorded:** the first run used the transfer parameters Pascale prints in his Table 6
  (γ 0.42, 4-decimal matrix) and failed at 204 counts. Its errors had one sign and vanished at black and
  white, which is the signature of an exponent error. With the IEC exponent 1/2.4 the error is 6
  counts, so his table was computed with 1/2.4 and the 0.42 is a rounded display. The tolerance was not
  moved. The 0.42 variant stays in the test as a control that must fail, and it does (204 > 40).

**Lindbloom, "Chromatic Adaptation" (brucelindbloom.com).**

- The reference whites for A, C, D50, D55 and D65 (2°, ASTM E308-01) match exactly.
- The Bradford matrices built from those whites match his published A, C, D50 and D55 → D65 matrices to
  **4.8e-8** per element.
- D55 was added for this campaign, because Hosoki 2007 measured under D55.

**Also tested:**

- The sRGB primaries round-trip exactly.
- Bradford from any white to itself is the identity. This is the null case.
- The float round trip through D50 is exact to 3e-13.
- Two calls with the same input give identical output.

**Not tested against a published table:** the 10° whites (D65_10, D50_10, C_10). They are carried over
from the first build. Only D65_10 is used, by the lip and palm entries.

Treating D65/10° data by Bradford rather than as if it were 2° moves a tissue colour by at most 0.125
ΔE\*ab. Adapting D55 data to D65 rather than ignoring the white moves it by 0.7–0.9 ΔE\*ab. The
adaptation is therefore not cosmetic, and a source's stated illuminant does matter.

## Cross-source consistency

Both sides of each pair below were carried to L\*a\*b\* under D65/2° by the module.

| comparison | ΔE\*ab | reading |
|---|---|---|
| gingiva: Ho 2015 (assumed D65/2°) vs Hosoki 2007 (D55) | 3.35 | Ho's assumed-illuminant value is corroborated to about 3 units by an independent instrument under a stated illuminant |
| palm: Horibata 2025 (assumed) vs Wang/Luo 2017 (D65/10°) | 4.13 | consistent; different cohorts |
| lip: Vergnaud 2024 multi-ethnic mean vs Hosoki 2007 Japanese lower lip | 7.95 | mostly L\*; the Vergnaud Caucasian mean (L\* 46.0) sits between them |
| tongue side: Cho 2025 (camera, no white point) vs Hosoki 2007 margin | 10.58 | device, not tissue: the reason a camera value is never mixed with an instrument value |

## Per role

Tiers are as of this build. **Used** is the source behind the colour. **Rejected / near miss** means a
number that was read and not used. **Searched, nothing** means the queries returned nothing
convertible.

### Mucosa and transitional zones

**`mucosa_lip_vermilion` — measured (unchanged).**
- Used: Vergnaud 2024, *Skin Res Technol* 30:e13583, n=410, hyperspectral, D65/10°. The published
  L\*C\*h is now converted in code (`lch` → `lch_to_lab`) instead of as precomputed a\*b\*. The hex is
  unchanged (#88564e).
- Corroboration: Hosoki 2007 lower lip, L\*49.6 a\*22.5 b\*14.3 (D55).
- Rejected: the Horibata 2025 lip value. It was measured **on an inkjet print of a photograph**, not on
  tissue. The earlier source record said only the conjunctiva was indirect, which was wrong, and it has
  been corrected.

**`mucosa_gingival` — measured, illuminant ASSUMED.**
- Used: Ho 2015, *Sci Rep* 5:18498, n=238, PR-670. The re-read found no illuminant or observer anywhere
  in the paper.
- Corroborated by Hosoki 2007, attached gingiva, L\*49.9 a\*24.9 b\*14.8, D55, n=62, to 3.35 ΔE.
- Not replaced by Hosoki, which states its illuminant but not its observer. Ho is the larger cohort,
  and neither source meets the full bar.

**`mucosa_palatal` — transferred (re-sourced).**
- Used: Hosoki M 2007, "Analysis of color changes of oral mucosa by smoking", *Kokubyo Gakkai Zasshi*
  74(2):108–118, doi:10.5357/koubyou.74.108. Read in full; Table 6 was read from the page image.
- Cohort: n=62 healthy nonsmokers, 30–83 y. Instrument: Konica Minolta CS-100 non-contact chroma meter,
  45/0 under an artificial-sunlight lamp. Illuminant stated as "illuminant D55"; observer not stated.
- Values, buccal mucosa: **L\*58.9 (3.16) a\*28.7 (3.03) b\*19.4 (4.60)** → #c4796d.
- Why it transfers: the soft palate and uvula carry lining mucosa (non-keratinized, loose vascular
  lamina propria), the same class as the cheek.
- Previously the palate borrowed keratinized, bound-down gingiva, which the old note itself called the
  wrong class. The new value is 8.4 ΔE away.
- Searched, nothing: hard or soft palate and palatal mucosa × CIELAB, SpectroShade, Easyshade,
  spectrophotometer (PubMed, Europe PMC full text). Deferm 2018 is HSV from an intraoral scanner, and
  only operator differences are reported.

**`mucosa_lingual` — transferred (re-sourced).**
- Used: Hosoki 2007, tongue margin, **L\*42.2 (4.21) a\*24.6 (3.00) b\*14.2 (4.18)**, D55 → #8e534e.
- Why it transfers: same organ and same healthy cohort, but the site is the lateral border rather than
  the papillated dorsum.
- Near misses: every healthy tongue-dorsum value found is camera-derived with no stated white point.
  - Cho 2025, *Oncol Rev* 19:1697252. b\* does exist, in Supplementary File 1: body L\*53.6 a\*25.3
    b\*11.1. Oncology cohort.
  - Kim 2024, *J Clin Med* 13:3549, TAS-4000, n=14: L\*50.8 a\*23.0 b\*12.7.
  - Noguchi 2023, *Sci Rep* 13:1334.
  - Tian 2024, n=1448: L\*69.5, b\* 4. Already rejected as a white-balance artefact.
- Could not read: Zeng 2011, *Zhong Xi Yi Jie He Xue Bao* 9:948 (spectrometer, 10°, n=516; the
  abstract gives x₁₀y₁₀Y but no illuminant); Takamoto 2013 (paywalled).

**`mucosa_pharyngeal` — synthesized → transferred.**
- Used: Hosoki 2007 buccal mucosa, as for the palate.
- Why it transfers: the oro- and laryngopharynx are lined by the same non-keratinized lining mucosa,
  continuous with the cheek. It does **not** transfer well to the nasopharynx (respiratory epithelium),
  which shares the role and is stated on the entry.
- The old synthesized "darker and redder than oral mucosa" offset had no measurement behind it and was
  dropped rather than applied.
- Searched, nothing: pharynx, oropharynx, posterior pharyngeal wall and tonsil × colorimetry, L\*a\*b\*,
  chromaticity, reflectance. The TXI pharynx study reports only ΔE between lesion and background.

**`mucosa_nasal` — synthesized.**
- Near miss: Joko 2002, *Am J Rhinol* 16:11 (abstract only). Inferior turbinate in 60 normal subjects,
  chromaticity x 0.4264 y 0.3204. There is no luminance and no illuminant, so it cannot give a colour.
- Respiratory epithelium over erectile turbinate is not the same class as oral lining mucosa, so the
  buccal value was **not** transferred here.

**`mucosa_gastric` — synthesized.**
- Every source is camera-derived with no stated white point or RGB space:
  - Oki 2026, *DEN Open* e70297: Photoshop Lab from EVIS X1 images; the background was eradication
    mucosa with map-like redness.
  - Mizukami 2017; Ishikawa 2021; Kanzaki 2023: ΔE only.

**`mucosa_glans` — transferred (unchanged).**
- Sommers 2013 labial mucosa, illuminant assumed. Nothing new found.

**`mucosa_anal_margin` — synthesized.**
- Searched, nothing: perianal, anoderm, anal margin/verge, anogenital, perineal × colorimetry,
  chromameter, L\*a\*b\*, Mexameter, erythema index. 21 PubMed hits, none relevant.

**`lymphoid_tonsil` — synthesized.**
- Searched with the pharynx, nothing found. The buccal value was not transferred: the tonsil's visible
  colour includes crypts and lymphoid tissue under a thin epithelium.

### Integument

**`skin`, `skin_region` — synthesized in the default build.**
- These take the selected skin-tone option. The measured Xiao 2017 options now carry "stated (D65, 2°)".
  That was confirmed from the accepted manuscript: "CIE XYZ … assuming a 2° standard observer … the
  illuminant is set to the CIE standard D65".
- The Sommers 2019 options carry "ASSUMED".

**`skin_palmoplantar` — measured (re-sourced).**
- Used: Wang Y, Luo MR, Wang M, Xiao K, Pointer M 2017, "Spectrophotometric measurement of human skin
  colour", *Color Res Appl* 42(6):764–774, doi:10.1002/col.22143. Read in full (accepted manuscript).
- Cohort: n=10 Chinese women. Instrument: Datacolor 600, de:8, 8 mm aperture.
- Illuminant and observer stated: "under CIE D65 illuminant and the CIE 1964 standard colorimetric
  observer".
- Value, Table 3: palm **L\*66.69 C\*16.92 h69.78°** → #b79e86. No SD is published. The male palm is
  darker by ΔL\* 3.60.
- Replaces Horibata 2025 thenar (n=67, no illuminant stated), which is 4.1 ΔE away.
- Searched, nothing, for the sole: Takiwaki 1994 is abstract only, and the ENCoDE plantar-toe spectra are
  behind PhysioNet credentialing. Plantar skin still carries the palmar value.

**`nail_plate` — measured, illuminant ASSUMED.**
- Horibata 2025 thumbnail; no illuminant stated.
- Near miss: Leeb 2024, *EBioMedicine* 102:105051, fingernail at D65/10°, but ITA only, with no
  L\*a\*b\*.

**`hair` — synthesized.**
- Nearest miss: Itou 2019, *Int J Mol Sci* 20:3739, and Itou 2022, 23:14459. CR-400 "with illuminant
  D65", observer not stated, per-subject tables, black Japanese hair only (e.g. L\*17.5 a\*4.2 b\*5.3).
- This fails only on the observer (about 0.1 ΔE) and on covering one hair colour. It is recorded as
  `rejected_measurement` on the entry, as the first option for the hair parameter the entry says should
  exist.
- Could not read: Mengel-From 2009 (ranges only), Lozano 2017, Vaughn 2008, Norton 2016.

### Musculoskeletal and connective

**`skeletal_muscle` — transferred (unchanged).**
- Bovine, illuminant C/2° stated.
- Considered: Knecht 2021, *Animals* 11:1282. Porcine, unbloomed, D65, observer unstated. Six muscles
  span L\* 46–60, so picking one would be picking a colour. Recorded on the entry.
- No human muscle colorimetry found, again.

**`adipose` — synthesized.**
- Rejected: Parkinson 2024, *Transl Anim Sci* 8:txae071. Bovine subcutaneous fat, CR-410, D65,
  L\*73.4 a\*10.2 b\*17–20. The authors report that the slaughter-plant steam cabinet lightens fat
  colour, and that they found this only after analysis. Recorded on the entry.
- No human fat colorimetry found.

**`cortical_bone` — synthesized.**
- Rubio 2020 (D65, "8°" observer, which is probably the geometry) shows its unheated controls only as
  bar charts. Schafer 2001 (defleshed skull) was already rejected.

**`cartilage_hyaline` — synthesized.**
- Could not read: Ishimoto 2009, *Osteoarthritis Cartilage* 17:1204. Human, spectrocolorimetric, the
  strongest lead; every route returned 403.
- The rabbit study from the same group (Hattori 2008) gives its values in figures only.

**`cartilage_fibro`, `tendon`, `ligament`, `fascia`, `serous_membrane` — synthesized.**
- Searched, nothing. Pilin 2007 (post-mortem tendon and disc, RGB/IHS) is abstract only and
  device-relative.

**`tooth` — measured (unchanged).**
- Wee 2023, D65/2° stated.

### Heart, vessels, blood, lymph

**`cardiac_muscle`, `valve_leaflet` — synthesized.**
- Searched, nothing: porcine and bovine heart and offal L\*a\*b\*; pericardium bioprosthesis colour.

**`artery_wall`, `vein_wall` — synthesized.**
- The only vascular colorimetry is coronary angioscopy of the inner plaque surface (Inami 2008). It is
  the wrong surface and video-derived.

**`blood_cavity` (and the four chamber overrides) — synthesized.**
- Near miss: the Basson 2021 impala dataset (figshare 12403112, the pre-review version). It has
  per-sample L\*a\*b\* against SaO₂ (e.g. predicted at 95%: L\*20.2 a\*44.2 b\*30.5), but no device,
  illuminant, observer or path length, and the species is impala.
- Could not read: Shibata 2012 (human haemodialysis blood, handheld colorimeter).

**`lymph_node`, `lymph_vessel` — synthesized.**
- A thin search, nothing found.

### Nervous

**`brain_grey`, `brain_white`, `peripheral_nerve`, `choroid_plexus`, `csf_space` — synthesized.**
- Fabelo 2018 has in-vivo human cortex hyperspectral data from 400 nm, but the spectra are only
  plotted. Its public database could give a colour if someone integrated normal-tissue pixels; that
  would be our computation, not a published colour.
- Stelzle 2011 (porcine nerve) stops at 650 nm.

### Viscera and glands

**`liver` — synthesized.**
- Near misses:
  - De 2007, MMVR 15: in-vivo porcine liver, camera plus chart, illuminant not stated, one animal.
    Its sRGB contradicts its own expert-matched range.
  - Kanamori 2021, *Lab Invest* 101:1098: murine, blood-flushed, ex vivo.
- Could not read: Nunes 2017 (in-vivo liver BRDF), McLaughlin 2010.

**`kidney` — synthesized.**
- Baran 2012 (fresh human nephrectomy) covers 630–800 nm only.

**`gallbladder`, `bile_duct` — synthesized.**
- Could not read: Maitland 1993, *Appl Opt* 32:586. Fresh human gallbladder and bile, 350–2450 nm. The
  best lead in the group; the abstract lists figures only.

**`spleen`, `pancreas`, `stomach_serosa`, `oesophagus`, `small_intestine`, `large_intestine`, `mesentery_omentum`, `lung`, `airway_wall`, `urothelium_wall`, `thymus` — synthesized.**
- Searched, nothing.
- Bowel: Clancy 2021 reports SO₂ and haemoglobin only.
- Peritoneum: van de Weerd 2026 uses TIVITA at 500–1000 nm.

**`gland_thyroid`, `gland_parathyroid` — synthesized.**
- Tseregorodtseva 2025, *Sci Rep* 15:22097, has calibrated in-vivo reflectance from 300 nm, but plotted
  only, with data "on reasonable request". This is the strongest human data request in the whole
  campaign.

**`gland_salivary`, `gland_serous` — synthesized.**
- Wisotzky 2020 has plots only. Could not read: Wisotzky 2018, eight human tissues, 400–700 nm.

**`gland_adrenal`, `gland_pituitary`, `gland_pineal` — synthesized.**
- Searched, nothing.

**`testis`, `epididymis`, `duct_muscular`, `seminal_vesicle`, `prostate`, `erectile_tissue` — synthesized.**
- Searched, nothing.

### Eye

**`sclera` — synthesized.**
- Every ocular-prosthesis source measures prosthesis materials or ocularist-made patches (e.g. Reinhard
  2024, *Nat Commun*, D50/2°), never the patient's eye.
- The jaundice studies are camera-derived and in sick cohorts.
- Could not read: Bashkatov 2010, ex-vivo optical coefficients.

**`iris` — synthesized.**
- Near miss: Edwards 2016, *Pigment Cell Melanoma Res* 29:141. n=1448, RAW photographs converted
  "with D55 / 2°" as a software setting, with no chart. Two identical camera bodies differ by 5–7 L\*,
  which is direct evidence the values are device-relative.
- Could not read: Melgosa 2000, *Ophthalmic Physiol Opt* 20:252 (spectroradiometric, n=40). The best
  lead for the iris.

**`retina` — synthesized.**
- Delori & Pflibsen 1989, *Appl Opt* 28:1061, tabulates in-vivo fundus reflectance, but from 445 nm
  and as "equivalent reflectance" through the ocular media. It was read only through a table preview.
  Not convertible without inventing the band below 445 nm.

**`lens`, `cornea`, `choroid`, `ciliary_body`, `ocular_humour` — synthesized.**
- Sparrow 1988 grades lens colour against Munsell chips by eye. Sasaki 1985 and Artigas 2012 could not
  be read.

**`field_overlay`** is not a tissue.

## Library requests that would move the count

Ranked by likely yield:

1. Maitland 1993 (gallbladder and bile spectra).
2. Tseregorodtseva 2025 raw spectra (thyroid and parathyroid; ask the authors).
3. Ishimoto 2009 (articular cartilage).
4. Melgosa 2000 (iris).
5. Zeng 2011 (tongue and lip, 10°).
6. Wisotzky 2018 (eight tissues, 400–700 nm).
7. Takiwaki 1994 (sole).
8. PhysioNet ENCoDE (palm, toe and earlobe spectra, credentialed).
9. Shibata 2012 (human blood colour).

## Record corrections made in passing

- **`horibata2025sites`:** the lip was measured on a print, like the conjunctiva. The conjunctiva
  entries in `UNMAPPED` and on `sclera` no longer say a tissue measurement exists.
- **`scripts/verify_tissue_colour_palettes.py`** asserted 7390 structures. The display manifest has
  8979, so the verifier failed every run. It now takes the count from the manifest it verifies.
- **`docs/WORKBENCH_AUTHENTICITY.md` Tier 7** still reads "5 measured, 4 transferred, 65 synthesized".
  It is outside this change's territory and is left for its owner. The current figure is 5 / 5 / 64,
  with 2 of the 5 measured resting on an assumed illuminant.
