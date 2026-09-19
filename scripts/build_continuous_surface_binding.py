"""Freeze the shared geometric skin embedding from retained registration.

One blend per anatomy plant (`ihm.assembly.articulated.SURFACE_BINDINGS`), every one by the
SAME `materialize` over the segment supports of that plant's own registration.json:

    # base, 22 segments -- the asset every result before 18 Sep 2026 used
    PYTHONPATH=. .venv/bin/python scripts/build_continuous_surface_binding.py \\
      --registration <retained-runtime>/mechanics/registration.json
    # articulated_spine_v1, 25 segments. A variant ArticulatedBodyPlant writes registration.json
    # (and canonical_mechanics.json beside it) before it loads the blend, so a construction that
    # refuses for want of this asset leaves behind exactly the input this needs.
    PYTHONPATH=. .venv/bin/python scripts/build_continuous_surface_binding.py \\
      --plant articulated_spine_v1 --registration <variant-plant-output>/registration.json

`--plant` selects only the default OUTPUT path and a check that the registration's bodies are
that plant's; the weights depend on nothing but the registration's groups and the skin.
"""
import argparse,gzip,json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET
from ihm.assembly.continuous_surface_binding import materialize
from ihm.assembly.articulated import SURFACE_BINDINGS
from ihm.assembly.anatomy_pose import PLANTS

if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    parser=argparse.ArgumentParser();parser.add_argument('--registration',type=Path,required=True)
    parser.add_argument('--plant',default='engineering_stance_v1',choices=sorted(SURFACE_BINDINGS))
    parser.add_argument('--output',type=Path,default=None);args=parser.parse_args()
    output=args.output or root/SURFACE_BINDINGS[args.plant]
    manifest=json.loads(args.registration.read_text());payload=json.loads((args.registration.parent/'canonical_mechanics.json').read_text())
    bodies={b.get('name') for b in ET.parse(root/PLANTS[args.plant][0]).getroot().iter('Body')}
    if set(manifest['groups'])!=bodies:
        raise SystemExit(f'registration registers {len(manifest["groups"])} bodies, plant {args.plant} declares {len(bodies)}: '
                         +', '.join(sorted(set(manifest['groups'])^bodies)))
    registration=SimpleNamespace(groups=manifest['groups'],specs={e['id']:e for e in payload['entities']},manifest=lambda:manifest)
    result=materialize(root,registration);raw=json.dumps(result,separators=(',',':'),allow_nan=False).encode();output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(gzip.compress(raw,mtime=0))
    print(json.dumps(dict(plant=args.plant,path=str(output),segments=len(result['segments']),binding_identity=result['binding_identity'],algorithm=result['algorithm'],bytes=len(raw),compressed_bytes=output.stat().st_size),indent=2))
