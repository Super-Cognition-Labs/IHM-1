"""Exact triangle -> skin-patch membership for the 1,326 dermatome patches.

WHY THIS EXISTS.  A nociceptor transduces the pressure at the skin where it sits.
The supine surface quadrature (`scripts/build_supine_surface_contact.py`) already
carries, per sample, the index of the skin triangle it lies on; the 1,326 patches
in `data/derived/canonical/dermatomes.json` carry a centroid, an area and a
triangle COUNT -- but not the triangles.  Assigning a contact sample to the patch
with the nearest centroid is the coarse proxy IHM-1's CLAUDE.md records failing
three times on a mesh with a long facet tail ("the nearest centroid is routinely
not the face the association sits on").  So this recovers the real membership.

HOW, AND WHY IT CAN BE TRUSTED.  `scripts/build_dermatome_patches.py` is
deterministic and computes the membership (`cl`, per patch) but does not save it.
This re-runs that builder READ-ONLY -- its output is redirected to a temporary
directory, the canonical file is never written -- while recording each cluster
`pca_split` returns.  Then it asserts the rebuilt `dermatomes.json` is BYTE-
IDENTICAL to the canonical one.  Byte identity is the known answer: if the rebuild
differs in any way the membership it produced is not the membership that made the
canonical patches, and the script refuses to write a binding.  Every patch's
triangle count and area are then re-checked against the canonical record.

OUTPUT (gitignored, regenerable): `data/derived/nociception/face_patch.npz`
    face_patch   int32 [n_faces]  patch index into dermatomes.json['patches'],
                                  -1 for faces outside the exterior selection
    keyed by the sha256 of dermatomes.json and of the skin geometry, so a loader
    can refuse a binding built for different inputs.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "data/derived/canonical"
OUT_DIR = ROOT / "data/derived/nociception"
BINDING = OUT_DIR / "face_patch.npz"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _load_builder():
    spec = importlib.util.spec_from_file_location(
        "_dermatome_builder", ROOT / "scripts/build_dermatome_patches.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def recover_membership() -> dict:
    mod = _load_builder()
    clusters: list[np.ndarray] = []
    original_split = mod.pca_split

    def recording_split(idx, C, area, target):
        out = original_split(idx, C, area, target)
        for cl in out:
            clusters.append((np.asarray(cl).copy(), float(area[cl].sum())))
        return out

    mod.pca_split = recording_split
    with tempfile.TemporaryDirectory(prefix="nocicept-binding-") as tmp:
        tmp = Path(tmp)
        # the builder reads anatomy.json and peripheral.json from OUT and writes
        # dermatomes.json there.  symlink the two inputs; the write lands in tmp.
        for name in ("anatomy.json", "peripheral.json"):
            os.symlink(CANON / name, tmp / name)
        mod.OUT = tmp
        mod.build(verbose=False)
        rebuilt = (tmp / "dermatomes.json").read_bytes()
    canonical = (CANON / "dermatomes.json").read_bytes()
    if rebuilt != canonical:
        raise SystemExit(
            "REFUSED: the rebuilt dermatomes.json is not byte-identical to the "
            "canonical one, so the recorded membership is not the membership that "
            f"produced it (rebuilt sha {hashlib.sha256(rebuilt).hexdigest()[:16]}, "
            f"canonical {hashlib.sha256(canonical).hexdigest()[:16]}).")
    # the builder skips clusters with zero area; the patches are the rest, in order
    kept = [cl for cl, a in clusters if a > 0]
    return dict(mod=mod, clusters=kept, canonical_sha=hashlib.sha256(canonical).hexdigest())


def build(verbose: bool = True) -> dict:
    derm = json.loads((CANON / "dermatomes.json").read_text())
    patches = derm["patches"]
    geo_path = ROOT / derm["skin_geometry"]["path"]
    if sha256(geo_path) != derm["skin_geometry"]["sha256"]:
        raise SystemExit("skin geometry changed under dermatomes.json")
    rec = recover_membership()
    kept = rec["clusters"]
    if len(kept) != len(patches):
        raise SystemExit(f"REFUSED: {len(kept)} clusters recovered for "
                         f"{len(patches)} patches")
    g = json.loads(gzip.decompress(geo_path.read_bytes()))
    V = np.asarray(g["positions"], float).reshape(-1, 3)
    F = np.asarray(g["indices"], int).reshape(-1, 3)
    terr = json.loads((ROOT / "data/research/engineered_skin_territories/"
                               "materialization.json").read_text())
    keep = np.asarray(terr["contact_eligible_triangle_ids"], int)
    P = V[F[keep]]
    tri_area = 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1)

    face_patch = np.full(len(F), -1, dtype=np.int32)
    worst_area = 0.0
    for i, (cl, p) in enumerate(zip(kept, patches)):
        if len(cl) != p["triangle_count"]:
            raise SystemExit(f"REFUSED: patch {p['id']} triangle count "
                             f"{len(cl)} != recorded {p['triangle_count']}")
        a = float(tri_area[cl].sum())
        worst_area = max(worst_area, abs(a - p["area_m2"]) / p["area_m2"])
        faces = keep[cl]
        if (face_patch[faces] != -1).any():
            raise SystemExit(f"REFUSED: patch {p['id']} overlaps an earlier patch")
        face_patch[faces] = i
    if worst_area > 1e-9:
        raise SystemExit(f"REFUSED: worst relative area mismatch {worst_area:.3e}")
    assigned = int((face_patch >= 0).sum())
    if assigned != len(keep):
        raise SystemExit(f"REFUSED: {assigned} faces assigned of {len(keep)} exterior")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(BINDING, face_patch=face_patch,
                        dermatomes_sha256=np.array(rec["canonical_sha"]),
                        geometry_sha256=np.array(derm["skin_geometry"]["sha256"]),
                        n_patches=np.array(len(patches)))
    report = dict(binding=str(BINDING.relative_to(ROOT)),
                  dermatomes_sha256=rec["canonical_sha"],
                  geometry_sha256=derm["skin_geometry"]["sha256"],
                  faces_total=int(len(F)), faces_bound=assigned,
                  patches=len(patches),
                  rebuilt_dermatomes_byte_identical=True,
                  worst_relative_area_mismatch=worst_area,
                  triangle_counts_match=True, patches_disjoint=True)
    (OUT_DIR / "face_patch_report.json").write_text(json.dumps(report, indent=2) + "\n")
    if verbose:
        print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    build(verbose=not ap.parse_args().quiet)
