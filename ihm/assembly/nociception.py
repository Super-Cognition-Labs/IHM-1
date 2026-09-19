"""Cutaneous nociceptor transduction at the 1,326 innervated skin patches.

Stimulus in -- normal contact pressure on the skin, Pa -- and afferent firing out --
A-delta and C, Hz -- delayed along each patch's own route by fibre class.  That is
the whole contract.  Read the three sections below before using the output for
anything, because each one is a place this module could be mistaken for more than
it is.

### what this is not

**Not a pain scale, not a reward, not a claim about experience.**  Nociceptor
firing and pain are different quantities and the human recordings this module is
built on say so directly: Van Hees and Gybels 1981 recorded mechanically evoked
C-nociceptor discharge "up to more than 10 spikes/s" that was not necessarily
painful, and Adriaensen et al. 1984 titled their paper "a paradox" because pain
rose through a 120 s squeeze while the polymodal C-fibre discharge adapted away.
Anything that turns these rates into a scalar to be maximised or minimised is a
separate decision, made elsewhere, with its own evidence.

### where every number comes from

Every constant below is in `EVIDENCE` with its source, the exact words or figure
it was read from, and how it is used.  The shape of the argument:

* **Threshold and half-saturation come from ONE stimulus geometry.**  Both human
  anchors are forceps squeezes on 30 mm2 faces, so the conversion from force to
  pressure is the same division for both and no constant crosses a change of
  units or of probe area:
    - pain onset: "forces greater than 4 N exerted on forceps faces of 30 mm2
      elicited pain" (Adriaensen et al. 1984) -> 4 N / 30 mm2 = 133.3 kPa;
    - strongly painful tonic pressure: "14N at 30 mm(2); 120 s", which evoked
      augmenting pain and recruited mechano-insensitive C units (Schmidt et al.
      2000) -> 14 N / 30 mm2 = 466.7 kPa.
* **The ceiling is an OBSERVED maximum, not a fitted asymptote**: mechanically
  evoked C discharge "up to more than 10 spikes/s" (Van Hees and Gybels 1981).
  "More than" makes 10 Hz a floor on what was seen.
* **The curve between them is a modelling choice, stated as one.**  The sources
  give a threshold, a strongly painful level and a ceiling; none gives a curve.
  The law is a rectified hyperbola, silent at and below threshold, half of the
  ceiling at the strong tonic level, approaching the ceiling above it:

      x = max(0, P - P_threshold) / (P_half - P_threshold)
      rate = rate_ceiling * x / (1 + x)

  Placing half-saturation at 466.7 kPa is the anchoring decision; it is not a
  measurement of where human nociceptors half-saturate.

* **An independent human check, not used in the law.**  Algometer pressure pain
  thresholds over a 1 cm2 probe (DFNS protocol, Rolke et al. 2006) in 130 healthy
  adults (Pan et al. 2024): thenar 290.3 +- 91.3 kPa (men), 256.7 +- 94.8 (women);
  masseter 178.5 +- 56.7 and 156.6 +- 58.4.  Pooled by n, 270.9 kPa thenar and
  165.9 kPa masseter -- 2.03x and 1.24x the forceps threshold.  Spatial summation
  (Lautenbacher et al. 2005) predicts the LARGER probe gives the LOWER threshold;
  the difference runs the other way, so site and loading mode dominate it, not
  area.  The honest reading is that human mechanical pain thresholds span
  roughly 130-290 kPa across these two instruments and sites, and the model sits
  at the low end of that span.

### the three things this does not model, each named

1. **Sub-pain nociceptor firing.**  The threshold is where humans REPORT pain.
   Nociceptor units fire below it (the von Frey thresholds of human polymodal C
   units are 2.3-13.1 g, Van Hees and Gybels 1981, and < 160 mN, Schmidt et al.
   2000).  Those thresholds are forces on hair tips whose area the abstracts do
   not give, so they cannot be converted to a pressure without inventing one,
   and they are recorded in `EVIDENCE` as NOT USED.  This channel is therefore
   silent over a range where real C-nociceptors are not: it understates sub-pain
   activity by construction.
2. **Adaptation and sensitisation.**  Both are measured -- polymodal C units give
   "an initial high frequency dynamic discharge followed by adaptation"
   (Adriaensen et al. 1984); mechano-insensitive units were activated only
   "after more than 20s" and responded more strongly to a second identical
   stimulus (Schmidt et al. 2000) -- and neither source gives a time constant.
   The transduction here is static.  Inventing the time constants would be the
   ledger's "constant typed from memory".
3. **A separate A-delta parameterisation.**  No human A-delta mechanical
   stimulus-response function was read for this module (Adriaensen et al. 1983,
   J Neurophysiol 49:111-122, is the human A-delta source; its content was not
   available to read).  The A-delta channel therefore uses the SAME law as C.
   That is an unmeasured assumption, recorded in `EVIDENCE['adelta_law']`, and
   what would measure it is that paper's high-threshold-mechanoreceptor data.
   What does differ between the classes, and is the point of routing them
   separately, is conduction velocity.

### stimulus sources

* **Supine surface quadrature** -- per-sample normal force and area, and the
  index of the skin triangle each sample lies on.  `patch_rates_from_quadrature`
  binds a sample to a patch by that triangle, through the exact membership
  recovered by `scripts/build_nociceptor_patch_binding.py`, never by nearest
  centroid.  It transduces at each sample and reports the patch's SKIN-AREA-
  WEIGHTED MEAN rate -- the mean firing of endings spread evenly over the patch,
  with unloaded skin contributing zero -- and the peak sample rate beside it.
  The quadrature's own areas are PROJECTED; `supine_sample_surface_area`
  converts them to skin area, and mixing the two is a bug this module shipped
  for one gate run (docs/NOCICEPTION.md).
* **A per-patch normal force** -- `patch_rates_from_patch_forces` divides by the
  patch area.  That is the patch-mean pressure and a LOWER BOUND on the peak the
  endings feel whenever the load is concentrated on part of the patch.
* **The crude body's contact elements are not converted.**  They are spheres
  carrying a force and no contact area; a pressure needs an area they do not
  have, and deepening dependence on the scaffold is what ACTUATION_STAGES.md
  says not to do.

### routing

`NociceptorDelayLine` delays each patch's rate by `path_length_m / v` for A-delta
and for C, using the fibre velocities IHM already carries as IBM's snapshot
(`peripheral.json['fibre_velocity_m_s']`, typical 15 and 1.0 m/s).  The human C
conduction velocity measured by Van Hees and Gybels, 0.86 +- 0.17 m/s, sits 0.8
SD below the 1.0 m/s used.  **The route lengths are schematic, not measured**:
each is an authored polyline from the patch centroid through a waypoint to a
relay (`length_method` on every patch), not a dissected nerve.  IBM-1's
`ibm/topologies/ihm_bridge.py` carries the standing caveat -- found when the
spinal relays sat 89-337 mm below their cord segments -- that every spinal-route
delay is a schematic LOWER BOUND.  The relays have since moved onto the dura
centreline (docs/BODY_PERIPHERAL.md: patch lengths x1.29 median) and the caveat
has not been lifted, so it is carried here unchanged: quote these delays as
schematic, never as measured conduction times.  Central and synaptic delays are
excluded (`route_contract` in peripheral.json).  Which trunk carries
which class is IBM's to say; the IBM join
(`ibm.topologies.ihm_bridge.cutaneous_nociceptive_routes`) asserts every patch
trunk carries both A-delta and C.
"""
from __future__ import annotations

import copy
import hashlib
import heapq
import json
import math
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# the published anchors.  every number used below is derived from these lines.
# ---------------------------------------------------------------------------

#: Adriaensen H, Gybels J, Handwerker HO, Van Hees J (1984) Nociceptor discharges
#: and sensations due to prolonged noxious mechanical stimulation -- a paradox.
#: Hum Neurobiol 3(1):53-58.  PMID 6330012.
ADRIAENSEN_1984 = dict(
    citation="Adriaensen H, Gybels J, Handwerker HO, Van Hees J (1984). Nociceptor "
             "discharges and sensations due to prolonged noxious mechanical "
             "stimulation -- a paradox. Hum Neurobiol 3(1):53-58.",
    pmid="6330012", species="human", method="psychophysics + microneurography, "
    "forceps squeeze of a dorsal-hand skin fold, constant force for 120 s",
    quote="forces greater than 4 N exerted on forceps faces of 30 mm2 elicited pain")

#: Schmidt R, Schmelz M, Torebjork HE, Handwerker HO (2000) Mechano-insensitive
#: nociceptors encode pain evoked by tonic pressure to human skin.
#: Neuroscience 98(4):793-800.  doi:10.1016/s0306-4522(00)00189-5, PMID 10891622.
SCHMIDT_2000 = dict(
    citation="Schmidt R, Schmelz M, Torebjork HE, Handwerker HO (2000). "
             "Mechano-insensitive nociceptors encode pain evoked by tonic pressure "
             "to human skin. Neuroscience 98(4):793-800.",
    doi="10.1016/s0306-4522(00)00189-5", pmid="10891622", species="human",
    method="microneurography, superficial peroneal nerve; tonic pressure",
    quote="tonic pressure stimulation (14N at 30 mm(2); 120 s) ... induced "
          "augmenting pain responses which were matched by the discharges of "
          "initially mechano-insensitive ('silent') C-units")

#: Van Hees J, Gybels J (1981) C nociceptor activity in human nerve during
#: painful and non painful skin stimulation.  J Neurol Neurosurg Psychiatry
#: 44(7):600-607.  doi:10.1136/jnnp.44.7.600, PMID 7288447.
VAN_HEES_GYBELS_1981 = dict(
    citation="Van Hees J, Gybels J (1981). C nociceptor activity in human nerve "
             "during painful and non painful skin stimulation. J Neurol Neurosurg "
             "Psychiatry 44(7):600-607.",
    doi="10.1136/jnnp.44.7.600", pmid="7288447", species="human",
    method="microneurography, >100 single C fibres, radial nerve",
    quotes=["mechanical threshold, measured with von Frey hairs, varied between "
            "2.3 and 13.1 g",
            "conduction velocities ... mean value of 0.86 m/s (SD: 0.17)",
            "discharge even up to more than 10 spikes/s ... not necessarily "
            "accompanied by pain sensation"])

#: Pan LH et al. (2024) The normative values of pain thresholds in healthy
#: Taiwanese.  Brain Behav 14(4):e3485.  doi:10.1002/brb3.3485, PMID 38648375.
PAN_2024 = dict(
    citation="Pan LH, Ling YH, Lai KL, Wang YF, Hsiao FJ, Chen SP, Liu HY, Chen WT, "
             "Wang SJ (2024). The normative values of pain thresholds in healthy "
             "Taiwanese. Brain Behav 14(4):e3485.",
    doi="10.1002/brb3.3485", pmid="38648375", species="human",
    method="electronic algometer, 1 cm2 probe, ~50 kPa/s, DFNS protocol",
    n_male=55, n_female=75,
    ppt_kpa={"thenar": {"male": (290.3, 91.3), "female": (256.7, 94.8)},
             "masseter": {"male": (178.5, 56.7), "female": (156.6, 58.4)}})

#: Rolke R et al. (2006) Quantitative sensory testing in the German Research
#: Network on Neuropathic Pain (DFNS): standardized protocol and reference
#: values.  Pain 123:231-243.  The protocol Pan et al. follow.
ROLKE_2006 = dict(
    citation="Rolke R, Baron R, Maier C, et al. (2006). Quantitative sensory testing "
             "in the German Research Network on Neuropathic Pain (DFNS): "
             "standardized protocol and reference values. Pain 123:231-243.",
    quote="pressure gauge device ... with a probe area of 1 cm2 (probe diameter of "
          "1.1 cm) ... slowly increasing ramp of 50 kPa/s")

#: Lautenbacher S, Kunz M, Strate P, Nielsen J, Arendt-Nielsen L (2005) Age
#: effects on pain thresholds, temporal summation and spatial summation of heat
#: and pressure pain.  Pain 115(3):410-418.  doi:10.1016/j.pain.2005.03.025.
LAUTENBACHER_2005 = dict(
    citation="Lautenbacher S, Kunz M, Strate P, Nielsen J, Arendt-Nielsen L (2005). "
             "Age effects on pain thresholds, temporal summation and spatial "
             "summation of heat and pressure pain. Pain 115(3):410-418.",
    doi="10.1016/j.pain.2005.03.025", pmid="15876494")

FORCEPS_FACE_AREA_M2 = 30e-6       # "30 mm2" in both Adriaensen 1984 and Schmidt 2000
PAIN_ONSET_FORCE_N = 4.0           # Adriaensen 1984
STRONG_TONIC_FORCE_N = 14.0        # Schmidt 2000
OBSERVED_C_MAX_RATE_HZ = 10.0      # Van Hees and Gybels 1981, "more than 10 spikes/s"

THRESHOLD_PA = PAIN_ONSET_FORCE_N / FORCEPS_FACE_AREA_M2          # 133,333 Pa
HALF_SATURATION_PA = STRONG_TONIC_FORCE_N / FORCEPS_FACE_AREA_M2  # 466,667 Pa
RATE_CEILING_HZ = OBSERVED_C_MAX_RATE_HZ

FIBRE_CLASSES = ("adelta", "c")


def pooled_ppt_kpa(site: str) -> float:
    """Pan et al. 2024's sex-specific means pooled by group size -- derived, so a
    reader can re-derive it rather than trust a typed number."""
    d = PAN_2024["ppt_kpa"][site]
    nm, nf = PAN_2024["n_male"], PAN_2024["n_female"]
    return (nm * d["male"][0] + nf * d["female"][0]) / (nm + nf)


EVIDENCE = {
    "threshold_pa": dict(
        value=THRESHOLD_PA, units="Pa",
        derivation="4 N / 30 mm2", source=ADRIAENSEN_1984,
        kind="human psychophysical pain threshold, used as the receptor channel's "
             "threshold",
        transfer="a pain-REPORT threshold, not a unit's firing threshold; human "
                 "nociceptors fire below it (see sub_pain_von_frey).  applied to "
                 "every patch region, though it was measured on dorsal hand skin."),
    "half_saturation_pa": dict(
        value=HALF_SATURATION_PA, units="Pa",
        derivation="14 N / 30 mm2", source=SCHMIDT_2000,
        kind="the one human suprathreshold pressure, in the same 30 mm2 geometry, "
             "reported strongly and increasingly painful",
        transfer="half-saturation is PLACED here.  that placement is the modelling "
                 "choice; the source does not report a half-saturation."),
    "rate_ceiling_hz": dict(
        value=RATE_CEILING_HZ, units="Hz",
        source=VAN_HEES_GYBELS_1981,
        kind="observed maximum of mechanically evoked human C-nociceptor discharge",
        transfer="'more than 10 spikes/s' makes 10 Hz a FLOOR on the observed "
                 "maximum; used as the asymptote, so the channel may understate "
                 "peak discharge."),
    "law": dict(
        form="rate = ceiling * x / (1 + x), x = max(0, P - threshold) / "
             "(half_saturation - threshold)",
        kind="modelling choice",
        why="the sources give a threshold, a strongly painful level and a ceiling, "
            "not a curve.  a rectified hyperbola is silent at and below threshold, "
            "monotone, and saturates without a second free constant."),
    "adelta_law": dict(
        kind="UNMEASURED ASSUMPTION",
        statement="the A-delta channel uses the C-fibre law above unchanged",
        what_would_measure_it="human A-delta high-threshold mechanoreceptor "
            "thresholds and discharge: Adriaensen H, Gybels J, Handwerker HO, "
            "Van Hees J (1983) J Neurophysiol 49(1):111-122, "
            "doi:10.1152/jn.1983.49.1.111 -- not read for this module."),
    "sub_pain_von_frey": dict(
        kind="NOT USED", source=[VAN_HEES_GYBELS_1981, SCHMIDT_2000],
        values="polymodal C units: 2.3-13.1 g (22.6-128.5 mN); < 160 mN",
        why="forces on von Frey hair tips of undeclared area; converting to a "
            "pressure would need a tip area neither abstract gives."),
    "algometry_crosscheck": dict(
        kind="independent human check, NOT used in the law",
        source=[PAN_2024, ROLKE_2006, LAUTENBACHER_2005],
        pooled_ppt_kpa={s: pooled_ppt_kpa(s) for s in ("thenar", "masseter")},
        ratio_to_threshold={s: pooled_ppt_kpa(s) * 1e3 / THRESHOLD_PA
                            for s in ("thenar", "masseter")},
        reading="1 cm2 algometry reads 1.24-2.03x the 30 mm2 forceps threshold; "
                "spatial summation predicts the opposite sign, so site and loading "
                "mode dominate the difference."),
    "unmodelled_dynamics": dict(
        kind="measured, NOT modelled",
        source=[ADRIAENSEN_1984, SCHMIDT_2000],
        what="polymodal C adaptation; mechano-insensitive recruitment after >20 s "
             "and sensitisation on repeat",
        why="no time constant is given by either source"),
    "c_velocity_crosscheck": dict(
        kind="human measurement beside the IBM velocity used",
        source=VAN_HEES_GYBELS_1981, measured_m_s=(0.86, 0.17), used_m_s=1.0),
}

ROUTE_SCOPE = ("schematic lower bound (IBM-1 ibm/topologies/ihm_bridge.py caveat): "
               "route length is an authored endpoint-to-relay polyline "
               "(dermatomes.json path_length_scope), not a dissected nerve; "
               "central and synaptic delays excluded")
SCOPE = ("cutaneous nociceptor transduction: contact pressure -> A-delta and C "
         "afferent rate.  NOT a pain scale and NOT a reward.")


# ---------------------------------------------------------------------------
# transduction: a pure function of pressure.
# ---------------------------------------------------------------------------

def firing_rate_hz(pressure_pa, *, threshold_pa: float = THRESHOLD_PA,
                   half_saturation_pa: float = HALF_SATURATION_PA,
                   ceiling_hz: float = RATE_CEILING_HZ) -> np.ndarray:
    """High-threshold, saturating receptor law.  Zero at and below threshold.

    Accepts a scalar or an array of pressures in Pa.  Compressive pressure is
    positive; a negative value (tension, or no contact) transduces to zero.
    """
    if not (0.0 < threshold_pa < half_saturation_pa) or not ceiling_hz > 0:
        raise ValueError("need 0 < threshold < half_saturation and ceiling > 0")
    p = np.asarray(pressure_pa, dtype=float)
    if not np.isfinite(p).all():
        raise ValueError("non-finite pressure")
    x = np.maximum(0.0, p - threshold_pa) / (half_saturation_pa - threshold_pa)
    return ceiling_hz * x / (1.0 + x)


# ---------------------------------------------------------------------------
# the patches and their exact triangle membership
# ---------------------------------------------------------------------------

def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_patches(root) -> dict:
    root = Path(root)
    path = root / "data/derived/canonical/dermatomes.json"
    data = json.loads(path.read_text())
    data["_sha256"] = _sha256(path)
    return data


def load_face_binding(root, dermatomes: dict | None = None) -> np.ndarray:
    """face index -> patch index, refusing a binding built for other inputs."""
    root = Path(root)
    dermatomes = load_patches(root) if dermatomes is None else dermatomes
    path = root / "data/derived/nociception/face_patch.npz"
    if not path.exists():
        raise FileNotFoundError(f"{path} absent; run "
                                "scripts/build_nociceptor_patch_binding.py")
    z = np.load(path)
    if str(z["dermatomes_sha256"]) != dermatomes["_sha256"]:
        raise ValueError("face binding was built for a different dermatomes.json")
    if str(z["geometry_sha256"]) != dermatomes["skin_geometry"]["sha256"]:
        raise ValueError("face binding was built for a different skin geometry")
    return z["face_patch"]


def supine_sample_surface_area(root, quadrature_dir) -> tuple[np.ndarray, dict]:
    """The SKIN area each supine quadrature sample stands for, in m2.

    The quadrature's `area_m2` is PROJECTED area: a 5 mm raster cell on the plane
    normal to the bed (`supine_contact.posterior_envelope`).  The skin under that
    cell is larger by 1 / |cos| of the angle between the face normal and the ray.
    Pooling a rate over a patch by projected area against the patch's TRUE surface
    area mixes two measures, and did: the first gate run pooled a uniform load on
    one patch to 7.0985 Hz when every sample on it read 7.0588 Hz -- a mean above
    its own maximum (docs/NOCICEPTION.md).

    The ray direction is the source x axis, carried to canonical coordinates by
    the manifest's `registration.source_to_canonical_ground`.  That map is checked,
    not trusted: every sample, carried back, must lie on the plane of the face it
    is recorded on.

    Two things this does NOT do, both tried or considered and both wrong.  It
    does not cap a sample at the area of the face it hit: a 5 mm cell covers 25 mm2
    and 81% of exterior skin faces are smaller (median 7.7 mm2), so its footprint
    spans many faces and only its centre ray hits the recorded one, and a per-face cap cut
    the posterior skin from 0.868 to 0.489 m2 -- BELOW the 0.535 m2 projected area
    it came from, which a projection can never exceed.  And it does not floor
    |cos|, which would need a constant nobody measured.  Grazing cells (minimum
    |cos| 0.006) and raster excess are handled where they matter, in the pooling:
    a patch whose samples sum past its own area is normalised by that sum, and
    counted.
    """
    import gzip
    root, qdir = Path(root), Path(quadrature_dir)
    man = json.loads((qdir / "manifest.json").read_text())
    q = np.load(qdir / "quadrature.npz")
    geo_rel = "data/derived/canonical/geometry/body-bp3d-FJ2810.json.gz"
    if man["source_files"].get(geo_rel) != _sha256(root / geo_rel):
        raise ValueError("quadrature was built on a different skin geometry")
    g = json.loads(gzip.decompress((root / geo_rel).read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3)
    F = np.asarray(g["indices"], int).reshape(-1, 3)
    T = np.asarray(man["registration"]["source_to_canonical_ground"], float)
    R, t = T[:3, :3], T[:3, 3]
    if np.abs(R.T @ R - np.eye(3)).max() > 1e-9:
        raise ValueError("source_to_canonical_ground is not a rotation")
    f = q["face_indices"]
    tri = V[F[f]]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    face_area = 0.5 * np.linalg.norm(n, axis=1)
    n = n / (2.0 * face_area[:, None])
    canon = q["reference_points_source_m"] @ R.T + t
    off_plane = np.abs(((canon - tri[:, 0]) * n).sum(1))
    if off_plane.max() > 1e-9:
        raise ValueError(f"samples do not lie on their faces after the recorded "
                         f"map (max {off_plane.max():.3e} m)")
    cos = np.abs(n @ R[:, 0])
    true = q["area_m2"] / np.maximum(cos, 1e-300)
    return true, dict(max_off_plane_m=float(off_plane.max()), min_cos=float(cos.min()),
                      samples_cos_below_0_1=int((cos < 0.1).sum()),
                      projected_area_m2=float(q["area_m2"].sum()),
                      surface_area_m2=float(true.sum()))


def patch_rates_from_quadrature(face_indices, normal_force_n, area_m2,
                                surface_area_m2, face_patch: np.ndarray,
                                patch_area_m2) -> dict:
    """Transduce each contact sample, then pool per patch by skin area.

    Per sample, as the supine surface quadrature carries them:
      `normal_force_n`  the column force along the bed normal (the +x component
                        of `supine_contact.foundation(...)['point_forces_n']`);
      `area_m2`         the PROJECTED cell area that force was computed over, so
                        `normal_force_n / area_m2` is the foundation's own column
                        pressure;
      `surface_area_m2` the SKIN area the sample stands for
                        (`supine_sample_surface_area`), the pooling weight.
    The patch rate is sum(rate_i x skin_area_i) / max(patch_area, sum skin_area_i):
    the mean over the patch's endings, spread evenly, with unloaded skin counting
    zero.  The max() is there because raster excess and grazing cells can make a
    patch's samples sum past the patch itself; such patches are normalised by
    their samples and COUNTED (`patches_overcovered`).  Either way the weights sum
    to at most 1, so a patch rate can never exceed the highest sample rate on it
    -- asserted, not assumed.  A sample on a face in no patch is COUNTED and
    reported, never silently dropped.
    """
    f = np.asarray(face_indices, dtype=int)
    fn = np.asarray(normal_force_n, dtype=float)
    a = np.asarray(area_m2, dtype=float)
    w = np.asarray(surface_area_m2, dtype=float)
    if not (f.shape == fn.shape == a.shape == w.shape) or f.ndim != 1:
        raise ValueError("per-sample arrays must be 1-D and the same length")
    if not (np.isfinite(fn).all() and np.isfinite(a).all() and np.isfinite(w).all()) \
            or (a <= 0).any() or (w <= 0).any():
        raise ValueError("finite forces and positive areas required")
    if f.min(initial=0) < 0 or f.max(initial=0) >= len(face_patch):
        raise ValueError("face index outside the skin mesh")
    patch_area = np.asarray(patch_area_m2, dtype=float)
    pressure = np.maximum(0.0, fn) / a
    rate = firing_rate_hz(pressure)
    patch = face_patch[f]
    bound = patch >= 0
    n_patches = len(patch_area)
    covered = np.zeros(n_patches)
    np.add.at(covered, patch[bound], w[bound])
    over = covered > patch_area
    mean = np.zeros(n_patches)
    np.add.at(mean, patch[bound], rate[bound] * w[bound])
    mean /= np.maximum(patch_area, covered)
    peak = np.zeros(n_patches)
    np.maximum.at(peak, patch[bound], rate[bound])
    if (mean > peak * (1 + 1e-12) + 1e-15).any():
        raise AssertionError("a patch mean exceeds its own maximum sample rate")
    peak_pressure = np.zeros(n_patches)
    np.maximum.at(peak_pressure, patch[bound], pressure[bound])
    return dict(rate_hz=mean, peak_sample_rate_hz=peak,
                peak_pressure_pa=peak_pressure, covered_fraction=covered / patch_area,
                patches_overcovered=int(over.sum()),
                worst_overcover_ratio=float((covered / patch_area).max(initial=0.0)),
                samples=int(len(f)), samples_unbound=int((~bound).sum()),
                unbound_force_n=float(np.maximum(0.0, fn[~bound]).sum()),
                samples_above_threshold=int((pressure > THRESHOLD_PA).sum()),
                pooling="skin-area-weighted mean over the patch; unloaded skin counts 0")


def patch_rates_from_patch_forces(normal_force_n: dict, patches: list[dict]) -> dict:
    """Per-patch normal force -> patch-mean pressure -> rate.  A LOWER BOUND on the
    peak pressure whenever the load covers only part of the patch."""
    index = {p["id"]: i for i, p in enumerate(patches)}
    unknown = set(normal_force_n) - set(index)
    if unknown:
        raise ValueError(f"unknown patch ids: {sorted(unknown)[:5]}")
    rate = np.zeros(len(patches))
    pressure = np.zeros(len(patches))
    for pid, force in normal_force_n.items():
        force = float(force)
        if not math.isfinite(force):
            raise ValueError(f"{pid}: non-finite force")
        i = index[pid]
        pressure[i] = max(0.0, force) / patches[i]["area_m2"]
    rate[:] = firing_rate_hz(pressure)
    return dict(rate_hz=rate, pressure_pa=pressure,
                pooling="patch-mean pressure: a lower bound on the peak")


# ---------------------------------------------------------------------------
# routing: one delay per (patch, fibre class)
# ---------------------------------------------------------------------------

class NociceptorDelayLine:
    """Delay each patch's rate along its own route, separately per fibre class.

    Deterministic: events are ordered by (arrival time, serial), and the serial
    is assigned in a fixed patch order.  `step` returns what has ARRIVED at the
    relay end of the route by the end of the step.
    """

    def __init__(self, patches: list[dict], fibre_velocity_m_s: dict):
        self.patches = patches
        self.ids = [p["id"] for p in patches]
        self.velocity = {}
        for cls in FIBRE_CLASSES:
            v = fibre_velocity_m_s[cls]
            self.velocity[cls] = float(v["typical"] if isinstance(v, dict) else v)
            if not self.velocity[cls] > 0:
                raise ValueError(f"non-positive {cls} velocity")
        self.delay_s = {cls: np.array([p["path_length_m"] / self.velocity[cls]
                                       for p in patches]) for cls in FIBRE_CLASSES}
        self.time_s = 0.0
        self.sent = {cls: np.zeros(len(patches)) for cls in FIBRE_CLASSES}
        self.arrived = {cls: np.zeros(len(patches)) for cls in FIBRE_CLASSES}
        self.events: list = []
        self.serial = 0

    @classmethod
    def from_root(cls, root, dermatomes: dict | None = None):
        root = Path(root)
        dermatomes = load_patches(root) if dermatomes is None else dermatomes
        per = json.loads((root / "data/derived/canonical/peripheral.json").read_text())
        return cls(dermatomes["patches"], per["fibre_velocity_m_s"])

    def step(self, dt_s: float, rate_hz) -> dict:
        dt = float(dt_s)
        if not (math.isfinite(dt) and dt > 0):
            raise ValueError("dt_s must be positive")
        rate = np.asarray(rate_hz, dtype=float)
        if rate.shape != (len(self.ids),) or not np.isfinite(rate).all() or (rate < 0).any():
            raise ValueError("one finite non-negative rate per patch required")
        start, end = self.time_s, self.time_s + dt
        for cls in FIBRE_CLASSES:
            changed = np.flatnonzero(np.abs(rate - self.sent[cls]) > 1e-12)
            for i in changed:
                self.serial += 1
                heapq.heappush(self.events, (start + self.delay_s[cls][i],
                                             self.serial, cls, int(i), float(rate[i])))
                self.sent[cls][i] = rate[i]
        while self.events and self.events[0][0] <= end + 1e-12:
            _, _, cls, i, value = heapq.heappop(self.events)
            self.arrived[cls][i] = value
        self.time_s = end
        nerve = {cls: {} for cls in FIBRE_CLASSES}
        for cls in FIBRE_CLASSES:
            for i in np.flatnonzero(self.arrived[cls]):
                n = self.patches[i]["nerve_id"]
                nerve[cls][n] = nerve[cls].get(n, 0.0) + float(self.arrived[cls][i])
        return dict(time_s=end,
                    arrived_hz={cls: self.arrived[cls].copy() for cls in FIBRE_CLASSES},
                    nerve_hz=nerve, pending_events=len(self.events),
                    route_scope=ROUTE_SCOPE, scope=SCOPE, biological_validation=False)

    def snapshot(self) -> "NociceptorDelayLine":
        return copy.deepcopy(self)
