#!/usr/bin/env python3
"""Can the fitted, local-depth soft tissue layer be built for EVERY segment?  Census and gates.

The battery (scripts/verify_soft_tissue.py) builds the fitted layers on the heel only.  An identity
offered to the plant (`plant_options.SOFT_TISSUE_LAYERS['layer_fitted_local_*']`) must build on every
segment it names, so this builds each one, in its OWN child process (so each peak RSS is that
segment's), one after another, under the caller's address-space cap.

GATES, written into this header before this file was first run
  F1  Every segment the identity names (all admitted skin records, less `refused_segments`) builds at
      5 mm, has a core (at least one base node), has no tet of non-positive volume, and puts every
      skin node and every interface node within 1e-3 x spacing of its surface (G2's bar, per segment).
  F2  `build_selected_layers` refuses the identity's refused segment (radius_l) whether it is asked for
      by name or implied by `bodies=None`, and builds calcn_l with surface 'fitted' and depth 'local'
      in its receipt.
MEASURED, not gated: DOF, tetrahedra, layer volume, build seconds, peak MB, minimum and 1st-percentile
element quality, per segment and for the whole body.

    PYTHONPATH=. prlimit --as=4294967296 .venv/bin/python -u scripts/census_soft_tissue_fitted.py
"""
import json
import resource
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from ihm.assembly import soft_tissue_layer as stl                  # noqa: E402

IDENTITY = 'layer_fitted_local_confined'
OUT = ROOT / 'data/derived/soft-tissue-layer-v1/census-fitted-local.json'


def child(body):
    from ihm.assembly.plant_options import SOFT_TISSUE_LAYERS
    spec = SOFT_TISSUE_LAYERS[IDENTITY]
    began = time.perf_counter()
    out = {'body': body}
    try:
        layer = stl.segment_layer(ROOT, body, spacing_m=spec['spacing_m'], mapping=spec['mapping'],
                                  surface=spec['surface'], depth=spec['depth'])
    except (ValueError, RuntimeError, MemoryError) as error:
        out['error'] = f'{type(error).__name__}: {str(error)[:300]}'
    else:
        fit = layer.meta['fit']
        out.update({'dof': layer.dof, 'tetrahedra': int(len(layer.tets)), 'nodes': int(len(layer.local)),
                    'base_nodes': int(layer.base.sum()),
                    'layer_volume_l': float(layer.volumes.sum() * 1e3),
                    'minimum_volume_m3': float(layer.volumes.min()), **fit})
    out['build_seconds'] = time.perf_counter() - began
    out['peak_mb'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print('RESULT ' + json.dumps(out, default=float))


def main():
    if len(sys.argv) > 2 and sys.argv[1] == '--child':
        child(sys.argv[2])
        return 0
    from ihm.assembly.plant_options import SOFT_TISSUE_LAYERS, SOFT_TISSUE_UNANCHORED
    spec = SOFT_TISSUE_LAYERS[IDENTITY]
    manifest = json.loads((ROOT / stl.BUNDLE / 'manifest.json').read_text())
    bodies = [r['body'] for r in manifest['records'] if r['body'] not in SOFT_TISSUE_UNANCHORED
              and r['body'] not in spec['refused_segments']]
    rows, bad = [], []
    for body in bodies:
        p = subprocess.run([sys.executable, '-u', __file__, '--child', body], capture_output=True, text=True, cwd=ROOT)
        line = [ln for ln in p.stdout.splitlines() if ln.startswith('RESULT ')]
        row = json.loads(line[0][7:]) if line else {'body': body, 'error': 'child died: ' + p.stderr[-300:]}
        rows.append(row)
        if 'error' in row:
            bad.append(body)
            print(f'  FAIL {body:10s} {row["error"]}', flush=True)
            continue
        bar = 1e-3 * spec['spacing_m']
        ok = (row['base_nodes'] > 0 and row['minimum_volume_m3'] > 0 and row['skin_residual_max_m'] <= bar
              and row['interface_residual_max_m'] <= bar)
        if not ok:
            bad.append(body)
        print(f'  {"ok  " if ok else "FAIL"} {body:10s} DOF {row["dof"]:7d} tets {row["tetrahedra"]:7d} base {row["base_nodes"]:6d} '
              f'volume {row["layer_volume_l"]:.3f} L  skin {row["skin_residual_max_m"]*1e6:8.3f} um  interface '
              f'{row["interface_residual_max_m"]*1e6:8.3f} um  quality min {row["minimum_quality"]:.1e} p1 '
              f'{row["quality_percentiles_1_5_50"][0]:.3f}  build {row["build_seconds"]:.1f} s  peak {row["peak_mb"]:.0f} MB',
              flush=True)
    f1 = not bad
    print(f'  {"PASS" if f1 else "FAIL"} F1  every named segment builds with a core, no flat tet, nodes on their surfaces'
          + (f'  failing {bad}' if bad else ''), flush=True)
    selection = {'soft_tissue': {'id': IDENTITY, 'bundle': spec['bundle'], 'mapping': spec['mapping'],
                                 'segments': [r['body'] for r in manifest['records'] if r['body'] not in SOFT_TISSUE_UNANCHORED]}}
    refused = 0
    for bodies_arg in (None, ['radius_l']):
        try:
            stl.build_selected_layers(ROOT, selection, bodies_arg)
        except ValueError as error:
            refused += 'refuses radius_l' in str(error)
    built = stl.build_selected_layers(ROOT, selection, ['calcn_l'])['calcn_l']
    f2 = refused == 2 and built.meta['surface'] == 'fitted' and built.meta['depth_rule'] == 'local'
    print(f'  {"PASS" if f2 else "FAIL"} F2  radius_l refused by name and by bodies=None ({refused}/2); calcn_l built fitted/local',
          flush=True)
    ok_rows = [r for r in rows if 'error' not in r]
    total = {'dof': sum(r['dof'] for r in ok_rows), 'tetrahedra': sum(r['tetrahedra'] for r in ok_rows),
             'layer_volume_l': sum(r['layer_volume_l'] for r in ok_rows),
             'build_seconds': sum(r['build_seconds'] for r in rows), 'peak_mb_max': max(r['peak_mb'] for r in rows)}
    print('  whole body: ' + json.dumps(total, default=float), flush=True)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({'identity': IDENTITY, 'rows': rows, 'total': total, 'F1': f1, 'F2': f2}, indent=1, default=float) + '\n')
    return 0 if f1 and f2 else 1


if __name__ == '__main__':
    sys.exit(main())
