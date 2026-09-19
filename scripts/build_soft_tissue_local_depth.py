#!/usr/bin/env python3
"""This body's soft-tissue depth, per VERTEX of each skin contact mesh, in the segment's own frame.

WHY.  `soft_tissue_layer` made its rigid core from ONE number per segment: the median of the
depth map over that segment's skin (`skin-layer-map-v1`, `layer.thickness_m`).  On calcn_l that
is 18.61 mm.  Measured on the same depth map, the plantar heel is 6-13.5 mm to the calcaneus and
the plantar forefoot 4-9 mm; the median is set by the sides and back of the heel.  Applied to the
forefoot, which is thinner than 2 x 18.6 mm everywhere, the one-number rule leaves it with NO core:
the metatarsals are deleted and the tissue under them is carried by a midfoot core 30 mm away.
The depth map already holds the local value; this puts it where the layer can read it.

HOW.  The depth points are vertices of the canonical skin (soft-tissue-depth-v1 asserts it).  The
skin bundle's meshes are those same vertices through `build_skin_contact_meshes.binding_registration`
and the segment's reference-pose frame.  The same map is applied here, by importing the bundle
builder's own function rather than re-deriving it.

KNOWN ANSWERS (gates; the run fails if either fails)
  A1  Every mapped depth point lands on a vertex of its segment's bundle mesh.  Bar 1e-6 m (a
      micron; the transform chain is float64 and the bundle meshes are its output, so anything
      above round-off is a wrong map).  Measured before this file was written: <= 5.0e-10 m.
  A2  The median of the mapped depths per segment is EXACTLY the bundle's declared thickness
      (the same numbers selected the same way; float equality, not a tolerance).

FIRST RUN (2026-09-18): A1 FAILED for radius_l and passed for the other 19.  One of radius_l's 306
depth points is 10.8 mm from every vertex of its mesh: canonical vertex 89206 carries binding weight
0.401 radius_l against 0.393 ulna_l, so the layer map (vertex argmax) counted it for the radius,
while the bundle (triangle argmax of the three vertices' mean) gave all three of its triangles to the
ulna.  The known answer found a real mismatch between the two partitions, not a bad map.  The bar is
not moved: the artefact is written with each segment's verdict, the failing segment is marked, and
`soft_tissue_layer.local_depth` refuses the local rule for it.  The exit status stays non-zero while
any segment fails, so nobody reads this run as clean.

Unsampled skin vertices take the depth of the nearest sampled vertex (Euclidean, parameter-free).
Cap vertices (the centroid fans `cap_boundaries` adds at each joint cut) are not skin and get NaN;
the faces that use them are marked `cap`.

    PYTHONPATH=. prlimit --as=4294967296 .venv/bin/python scripts/build_soft_tissue_local_depth.py
"""
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'scripts'))

from ihm.assembly import soft_tissue_layer as stl        # noqa: E402

DEPTH = 'data/derived/soft-tissue-depth-v1/depth.npz'
OUT = ROOT / stl.LOCAL_DEPTH


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    import build_skin_contact_meshes as builder
    bundle = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_text())
    if bundle['registration']['choice'] != 'binding' or bundle['registration'].get('warp'):
        raise SystemExit('this bundle was not built on the unwarped binding map; the map below would be wrong')
    inverse_similarity, frames, _ = builder.binding_registration()
    mechanics = json.loads((ROOT / 'data/derived/canonical/mechanics.json').read_text())
    skin = next(e for e in mechanics['entities'] if e['role'] == 'skin')
    geometry = json.loads(gzip.decompress((ROOT / skin['reference_geometry']['path']).read_bytes()))
    canonical = np.asarray(geometry['positions'], float).reshape(-1, 3)
    binding = json.loads(gzip.decompress((ROOT / builder.BINDING).read_bytes()))
    segments = [s['id'] for s in binding['segments']]
    owner = np.asarray(binding['weights'], np.float32).argmax(1)
    z = np.load(ROOT / DEPTH)
    gap, vertex = cKDTree(canonical).query(z['skin_points_m'])
    if gap.max() != 0.0:
        raise SystemExit('depth points are not vertices of the canonical skin')
    source = z['skin_points_m'] @ inverse_similarity[:3, :3].T + inverse_similarity[:3, 3]
    canonical_source = canonical @ inverse_similarity[:3, :3].T + inverse_similarity[:3, 3]
    out, failures = {}, []
    for record in bundle['records']:
        body = record['body']
        v, f = stl.load_obj(ROOT / stl.BUNDLE / 'meshes' / record['mesh_file'])
        world = np.linalg.inv(frames[body])
        chosen = owner[vertex] == segments.index(body)
        local = source[chosen] @ world[:3, :3].T + world[:3, 3]
        depth = z['depth_m'][chosen]
        tree = cKDTree(v)
        miss, on = tree.query(local)
        a1 = float(miss.max())
        a2 = float(np.median(depth)) == record['layer']['thickness_m']
        # which mesh vertices are canonical skin (the rest are cap centroids)
        mine = canonical_source[owner == segments.index(body)] @ world[:3, :3].T + world[:3, 3]
        skin_gap, _ = cKDTree(mine).query(v)
        # a vertex off the segment's own canonical skin is a cap centre OR a neighbour's vertex the
        # partition gave this piece; cap centres are the ones no canonical vertex of ANY segment hits
        any_gap, _ = cKDTree(canonical_source @ world[:3, :3].T + world[:3, 3]).query(v)
        cap_vertex = any_gap > 1e-6
        sampled = np.zeros(len(v), bool)
        value = np.full(len(v), np.nan)
        # several depth points can share a vertex only if the depth map sampled it twice; keep the
        # first in the depth map's own order, so the choice is fixed by the data
        order = np.argsort(on, kind='stable')
        first = order[np.r_[True, on[order][1:] != on[order][:-1]]]
        value[on[first]] = depth[first]
        sampled[on[first]] = True
        fill = ~sampled & ~cap_vertex
        if fill.any():
            _, nearest = cKDTree(v[sampled]).query(v[fill])
            value[fill] = value[sampled][nearest]
        cap_face = cap_vertex[f].any(axis=1)
        ok = a1 <= 1e-6 and a2
        off_mesh = int((miss > 1e-6).sum())
        if not ok:
            failures.append(body)
        print(f'  {"PASS" if ok else "FAIL"} {body:10s} points {chosen.sum():5d} -> vertices {sampled.sum():5d} of {len(v):5d} '
              f'(cap {cap_vertex.sum()}, filled {fill.sum()}); A1 max gap {a1:.2e} m; A2 median {np.median(depth)*1e3:.4f} '
              f'== declared {record["layer"]["thickness_m"]*1e3:.4f} mm: {a2}; own-skin gap max {skin_gap[~cap_vertex].max():.1e}; '
              f'skin faces {int((~cap_face).sum())}, cap faces {int(cap_face.sum())}', flush=True)
        out[body] = {'mesh_file': record['mesh_file'], 'mesh_sha256': record['written_sha256'],
                     'vertex_depth_m': [None if not np.isfinite(x) else float(x) for x in value],
                     'sampled_vertex': sampled.astype(int).tolist(), 'cap_face': cap_face.astype(int).tolist(),
                     'depth_points': int(chosen.sum()), 'A1_max_gap_m': a1, 'A2_median_equals_declared': bool(a2),
                     'depth_points_off_this_mesh': off_mesh, 'known_answers_pass': bool(ok),
                     'declared_thickness_m': record['layer']['thickness_m']}
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = {'schema': 'ihm.soft-tissue-local-depth.v1', 'bundle': stl.BUNDLE,
                'bundle_manifest_sha256': sha(ROOT / stl.BUNDLE / 'manifest.json'),
                'depth_map': DEPTH, 'depth_map_sha256': sha(ROOT / DEPTH),
                'map': 'build_skin_contact_meshes.binding_registration, then the segment frame at the binding reference pose',
                'fill': 'unsampled skin vertex: depth of the nearest sampled vertex; cap vertex: none (not skin)',
                'known_answers': 'A1 every depth point on a bundle vertex (<= 1e-6 m); A2 per-segment median == declared thickness exactly',
                'failed_known_answers': failures,
                'bodies': out}
    (OUT / 'manifest.json').write_text(json.dumps(manifest) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}/manifest.json')
    if failures:
        raise SystemExit('FAILED known answers (recorded in the artefact; the local rule is refused for them): '
                         + ', '.join(failures))
    print('ALL KNOWN ANSWERS PASS')


if __name__ == '__main__':
    main()
