"""Gate the cutaneous nociceptor: known answers, idempotence, exact binding, routing.

Every check prints PASS or FAIL and the verdict is written as measured to
`out/gate_nociception.json`.  A failed check is recorded as FAILED; the bars below
were written before the first run and are not to be moved after one.

  K1  sub-threshold is silent      P in {0, 0.5 P0, P0(1-1e-9), P0}  -> exactly 0 Hz
  K2  half-saturation is placed    P = P_half                         -> ceiling / 2
  K3  far above saturates          P = 1e4 P0 -> within 0.1% of the ceiling, never above
  K4  monotone                     non-decreasing over 2,000 log-spaced pressures
  I1  idempotent law               firing_rate_hz twice at one input -> array_equal
  I2  idempotent pooling           patch_rates_from_quadrature twice -> equal
  I2s idempotent skin area         supine_sample_surface_area twice -> equal
  I3  idempotent routing step      same state, same input, stepped twice -> equal
  B1  binding is the patching      per-patch face count == dermatomes triangle_count
  B2  supine samples all bind      every supine quadrature sample lands in a patch
  B3  pooling known answer         uniform P on one patch's real samples -> the
                                   patch rate equals rate(P) x covered skin/patch area
  B3m a mean cannot exceed its max uniform P on EVERY supine sample -> no patch
                                   rate above rate(P)
                                   (added after run 1: see docs/NOCICEPTION.md)
  B4  null case on real geometry   0.999 P0 on every supine sample -> every patch 0
  S1  the projected->skin map      every sample, carried back by the manifest's
                                   rigid map, lies on its own face (< 1e-9 m)
  S2  a projection never grows     total skin area >= total projected area
                                   (added after run 2: see docs/NOCICEPTION.md)
  R1  A-delta before C             step stimulus on all 1,326 patches: first A-delta
                                   arrival strictly earlier than first C, every patch
  R2  delay is length / velocity   each first arrival within one dt of L / v
  R3  nothing arrives from nothing zero stimulus -> zero events, zero arrivals
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly import nociception as N  # noqa: E402

RESULTS: list[dict] = []


def record(key: str, ok: bool, detail: str, **values) -> bool:
    RESULTS.append(dict(check=key, verdict="PASS" if ok else "FAIL", detail=detail,
                        **values))
    print(f"  {key:3s} {'PASS' if ok else 'FAIL'}  {detail}")
    return ok


def _supine_quadrature():
    """the newest supine surface quadrature built from the CURRENT mechanics.json."""
    import hashlib
    mech = hashlib.sha256((ROOT / "data/derived/canonical/mechanics.json")
                          .read_bytes()).hexdigest()
    found = []
    for d in sorted((ROOT / "data/derived").glob("supine-surface-contact-*")):
        m = json.loads((d / "manifest.json").read_text())
        if m["source_files"].get("data/derived/canonical/mechanics.json") == mech:
            found.append((d.stat().st_mtime, d))
    if not found:
        raise SystemExit("no supine surface quadrature built from the current "
                         "mechanics.json; run scripts/build_supine_surface_contact.py")
    d = max(found)[1]
    return d, np.load(d / "quadrature.npz")


def main() -> int:
    P0, PH, R = N.THRESHOLD_PA, N.HALF_SATURATION_PA, N.RATE_CEILING_HZ
    print(f"threshold {P0:.1f} Pa, half-saturation {PH:.1f} Pa, ceiling {R} Hz")

    # -- the law ----------------------------------------------------------------
    sub = np.array([0.0, 0.5 * P0, P0 * (1 - 1e-9), P0])
    r = N.firing_rate_hz(sub)
    record("K1", bool((r == 0).all()), f"rates {r.tolist()} Hz below/at threshold",
           pressures_pa=sub.tolist(), rates_hz=r.tolist())
    rh = float(N.firing_rate_hz(PH))
    record("K2", abs(rh - R / 2) < 1e-12, f"rate at P_half = {rh!r} Hz (ceiling/2 = {R/2})",
           rate_hz=rh)
    far = float(N.firing_rate_hz(1e4 * P0))
    record("K3", R * 0.999 <= far <= R, f"rate at 1e4 P0 = {far:.6f} Hz of {R}",
           rate_hz=far)
    grid = np.logspace(3, 9, 2000)
    rg = N.firing_rate_hz(grid)
    record("K4", bool((np.diff(rg) >= 0).all()), "non-decreasing over 1 kPa..1 GPa")

    # -- idempotence of the law -------------------------------------------------
    probe = np.random.default_rng(0).uniform(0, 2e6, 4096)
    record("I1", np.array_equal(N.firing_rate_hz(probe), N.firing_rate_hz(probe)),
           "firing_rate_hz called twice on 4,096 pressures")

    # -- the exact binding on the real skin -------------------------------------
    derm = N.load_patches(ROOT)
    patches = derm["patches"]
    face_patch = N.load_face_binding(ROOT, derm)
    counts = np.bincount(face_patch[face_patch >= 0], minlength=len(patches))
    want = np.array([p["triangle_count"] for p in patches])
    record("B1", bool(np.array_equal(counts, want)),
           f"{int((counts == want).sum())} of {len(patches)} patches match their "
           f"recorded triangle count", faces_bound=int((face_patch >= 0).sum()))

    qdir, q = _supine_quadrature()
    faces, area = q["face_indices"], q["area_m2"]
    patch_area = np.array([p["area_m2"] for p in patches])
    bound = face_patch[faces] >= 0
    record("B2", bool(bound.all()),
           f"{int(bound.sum())} of {len(faces)} supine samples bound "
           f"({qdir.name})", quadrature=str(qdir.relative_to(ROOT)))

    skin, sdiag = N.supine_sample_surface_area(ROOT, qdir)
    record("S1", sdiag["max_off_plane_m"] < 1e-9,
           f"max sample off its face {sdiag['max_off_plane_m']:.2e} m; projected "
           f"{sdiag['projected_area_m2']:.4f} m2 -> skin {sdiag['surface_area_m2']:.4f} m2 "
           f"(min |cos| {sdiag['min_cos']:.4f}, {sdiag['samples_cos_below_0_1']} samples "
           f"below 0.1)", **sdiag)
    record("S2", sdiag["surface_area_m2"] >= sdiag["projected_area_m2"],
           "skin area >= the projected area it came from (a projection never grows)")

    # B3: a known uniform pressure on ONE patch's real samples
    target = int(np.bincount(face_patch[faces][bound]).argmax())
    on = face_patch[faces] == target
    P = 2.0 * PH
    force = np.where(on, P * area, 0.0)
    out = N.patch_rates_from_quadrature(faces, force, area, skin, face_patch, patch_area)
    expect = float(N.firing_rate_hz(P)) * skin[on].sum() / max(patch_area[target],
                                                               skin[on].sum())
    got = float(out["rate_hz"][target])
    others = float(np.abs(np.delete(out["rate_hz"], target)).max())
    record("B3", abs(got - expect) < 1e-12 and others == 0.0 and got <= float(N.firing_rate_hz(P)),
           f"patch {patches[target]['id']}: {got:.9f} Hz, expected {expect:.9f}, "
           f"sample rate {float(N.firing_rate_hz(P)):.9f}; covered {skin[on].sum()*1e4:.2f} "
           f"of {patch_area[target]*1e4:.2f} cm2 skin; every other patch {others}",
           patch=patches[target]["id"], rate_hz=got, expected_hz=expect,
           peak_sample_rate_hz=float(out["peak_sample_rate_hz"][target]))
    allP = N.patch_rates_from_quadrature(faces, P * area, area, skin, face_patch, patch_area)
    rP = float(N.firing_rate_hz(P))
    cf = allP["covered_fraction"]
    loaded = cf > 0
    record("B3m", bool((allP["rate_hz"] <= rP * (1 + 1e-12)).all()),
           f"uniform 2 P_half everywhere: max patch rate {allP['rate_hz'].max():.6f} Hz "
           f"<= sample rate {rP:.6f}; {int(loaded.sum())} patches loaded, covered "
           f"fraction median {np.median(cf[loaded]):.3f}, max {cf.max():.3f}; "
           f"{allP['patches_overcovered']} patches over-covered (normalised by samples)",
           patches_loaded=int(loaded.sum()), patches_overcovered=allP["patches_overcovered"],
           worst_overcover_ratio=allP["worst_overcover_ratio"],
           covered_fraction_median=float(np.median(cf[loaded])))
    out2 = N.patch_rates_from_quadrature(faces, force, area, skin, face_patch, patch_area)
    record("I2", all(np.array_equal(out[k], out2[k]) for k in
                     ("rate_hz", "peak_sample_rate_hz", "peak_pressure_pa",
                      "covered_fraction")),
           "patch_rates_from_quadrature called twice on the same field")
    skin2, _ = N.supine_sample_surface_area(ROOT, qdir)
    record("I2s", bool(np.array_equal(skin, skin2)),
           "supine_sample_surface_area called twice on the same quadrature")

    null = N.patch_rates_from_quadrature(faces, 0.999 * P0 * area, area, skin,
                                         face_patch, patch_area)
    record("B4", bool((null["rate_hz"] == 0).all()) and null["samples_above_threshold"] == 0,
           f"0.999 P0 on all {len(faces)} samples -> max patch rate "
           f"{null['rate_hz'].max()} Hz")

    # -- routing ----------------------------------------------------------------
    line = N.NociceptorDelayLine.from_root(ROOT, derm)
    dt = 1e-3
    zero = np.zeros(len(patches))
    quiet = line.snapshot()
    total_events = 0
    for _ in range(50):
        o = quiet.step(dt, zero)
        total_events += o["pending_events"]
    arrived_any = any(o["arrived_hz"][c].any() for c in N.FIBRE_CLASSES)
    record("R3", total_events == 0 and not arrived_any and quiet.serial == 0,
           f"zero stimulus: {quiet.serial} events sent, arrivals {arrived_any}")

    step_rate = np.full(len(patches), float(N.firing_rate_hz(2.0 * PH)))
    a = line.snapshot()
    b = line.snapshot()
    oa, ob = a.step(dt, step_rate), b.step(dt, step_rate)
    same = (oa["pending_events"] == ob["pending_events"] and
            all(np.array_equal(oa["arrived_hz"][c], ob["arrived_hz"][c])
                for c in N.FIBRE_CLASSES) and oa["nerve_hz"] == ob["nerve_hz"])
    record("I3", same, "one delay-line state stepped twice with one input")

    first = {c: np.full(len(patches), np.nan) for c in N.FIBRE_CLASSES}
    t_end = max(line.delay_s["c"]) + 10 * dt
    steps = int(np.ceil(t_end / dt))
    for k in range(steps):
        o = line.step(dt, step_rate)
        for c in N.FIBRE_CLASSES:
            new = np.isnan(first[c]) & (o["arrived_hz"][c] > 0)
            first[c][new] = o["time_s"]
    arrived = ~np.isnan(first["adelta"]) & ~np.isnan(first["c"])
    lead = first["c"] - first["adelta"]
    ok_r1 = bool(arrived.all() and (lead > 0).all())
    record("R1", ok_r1,
           f"{int((lead > 0).sum())} of {len(patches)} patches: A-delta first; "
           f"lead min {np.nanmin(lead)*1e3:.1f} ms, median {np.nanmedian(lead)*1e3:.1f} ms, "
           f"max {np.nanmax(lead)*1e3:.1f} ms",
           lead_ms={"min": float(np.nanmin(lead) * 1e3),
                    "median": float(np.nanmedian(lead) * 1e3),
                    "max": float(np.nanmax(lead) * 1e3)})
    err = {c: np.abs(first[c] - line.delay_s[c]) for c in N.FIBRE_CLASSES}
    worst = max(float(np.nanmax(err[c])) for c in N.FIBRE_CLASSES)
    record("R2", worst <= dt + 1e-9, f"worst |first arrival - L/v| = {worst*1e3:.3f} ms "
           f"(dt {dt*1e3:.0f} ms)")

    delays = {c: {"min_ms": float(line.delay_s[c].min() * 1e3),
                  "median_ms": float(np.median(line.delay_s[c]) * 1e3),
                  "max_ms": float(line.delay_s[c].max() * 1e3),
                  "velocity_m_s": line.velocity[c]} for c in N.FIBRE_CLASSES}
    for c, d in delays.items():
        print(f"  {c:6s} v={d['velocity_m_s']:5.1f} m/s  delay min {d['min_ms']:7.1f}  "
              f"median {d['median_ms']:7.1f}  max {d['max_ms']:7.1f} ms  "
              f"(schematic lower bound)")

    declared = [ln.split()[0] for ln in __doc__.splitlines()
                if ln.startswith("  ") and ln.split() and ln.split()[0][:1] in "KIBSR"
                and any(ch.isdigit() for ch in ln.split()[0])]
    ran = [r["check"] for r in RESULTS]
    record("ALL", sorted(declared) == sorted(ran),
           f"every declared check ran: declared {len(declared)}, ran {len(ran)}; "
           f"missing {sorted(set(declared) - set(ran)) or '-'}")
    failed = [r["check"] for r in RESULTS if r["verdict"] == "FAIL"]
    report = dict(gate="cutaneous nociceptor transduction", scope=N.SCOPE,
                  route_scope=N.ROUTE_SCOPE, verdict="FAILED" if failed else "PASSED",
                  failed=failed, checks=RESULTS, delays=delays,
                  law=dict(threshold_pa=P0, half_saturation_pa=PH, ceiling_hz=R),
                  evidence=N.EVIDENCE, dermatomes_sha256=derm["_sha256"],
                  patches=len(patches), supine_samples=int(len(faces)))
    path = ROOT / "out/gate_nociception.json"
    path.write_text(json.dumps(report, indent=1, default=lambda o: getattr(
        o, "tolist", lambda: str(o))()) + "\n")
    print(f"\n{report['verdict']}: {len(RESULTS) - len(failed)} of {len(RESULTS)} "
          f"checks pass; written to {path.relative_to(ROOT)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
