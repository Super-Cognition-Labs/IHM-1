#!/usr/bin/env python3
"""Build a matched CONTROL kernel from a retained IHM-1 kernel bundle.

`build_ibm_kernel_bundle.py --shuffle-seed` needs the 1 GB donor checkpoint in IBM-1.
The retained bundle already holds the exact `dyn.embed` the donor contributed, so a
control derived from it is the same control without the donor: the retained kernel's
own digest is recorded as the parent, and the donor provenance is carried over.

Two controls, each destroying a different thing:

* ``permute``  -- site rows permuted under a seed. Every per-site vector, its norm and
  the whole marginal distribution are kept; which site sits where (the learned
  association between sites) is destroyed. Identical to the donor builder's control.
* ``random``   -- every row replaced by an i.i.d. Gaussian direction rescaled to that
  row's original norm. Norms are kept per site; every direction, and so every learned
  relation between site vectors, is destroyed. A kernel that was never trained.

docs/STANCE_KERNEL_CONTROLS.md is the experiment these exist for.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    import torch
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--parent', default='data/models/ibm_curriculum16_kernel_v1')
    ap.add_argument('--mode', choices=('permute', 'random'), required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    parent = ROOT / args.parent
    parent_manifest = json.loads((parent / 'manifest.json').read_text())
    raw = (parent / 'kernel.pt').read_bytes()
    if sha(raw) != parent_manifest['kernel_sha256']:
        raise ValueError('Parent kernel does not match its own manifest')
    if parent_manifest.get('is_ablation_control'):
        raise ValueError('Parent is itself a control; derive controls from the real kernel')
    checkpoint = torch.load(parent / 'kernel.pt', map_location='cpu', weights_only=True)
    embed = checkpoint['dyn.embed']
    generator = torch.Generator().manual_seed(args.seed)
    if args.mode == 'permute':
        permutation = torch.randperm(embed.shape[0], generator=generator)
        if torch.equal(permutation, torch.arange(embed.shape[0])):
            raise ValueError('Degenerate control permutation')
        retained = embed[permutation]
        basis = ('Site rows permuted. Per-site embedding vectors, their norms and the whole marginal '
                 'distribution are identical to the parent; only which site sits where, and therefore '
                 'the learned association between sites, is destroyed.')
    else:
        noise = torch.randn(embed.shape, generator=generator, dtype=torch.float64).to(embed.dtype)
        norms = embed.norm(dim=1, keepdim=True)
        retained = noise / noise.norm(dim=1, keepdim=True) * norms
        if not torch.allclose(retained.norm(dim=1), embed.norm(dim=1), rtol=1e-5, atol=1e-7):
            raise ValueError('Random control does not preserve per-site norms')
        basis = ('Every site row replaced by an i.i.d. Gaussian direction rescaled to that row\'s '
                 'original norm. Per-site norms are kept; every direction, and so every learned '
                 'relation between site vectors, is destroyed: an untrained kernel.')
    kernel = dict(checkpoint)
    kernel['dyn.embed'] = retained.contiguous()
    out = ROOT / args.output
    out.mkdir(parents=True, exist_ok=False)
    torch.save(kernel, out / 'kernel.pt')
    kernel_bytes = (out / 'kernel.pt').read_bytes()
    reloaded = torch.load(out / 'kernel.pt', map_location='cpu', weights_only=True)
    if not torch.equal(reloaded['dyn.embed'], kernel['dyn.embed']):
        raise ValueError('Control kernel does not round-trip under the runtime loader')
    (out / 'pretrain_video_loop.py').write_bytes((parent / 'pretrain_video_loop.py').read_bytes())
    (out / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    manifest = {k: v for k, v in parent_manifest.items()}
    manifest.update({
        'kernel_sha256': sha(kernel_bytes),
        'kernel_bytes': len(kernel_bytes),
        'parent_bundle': args.parent,
        'parent_kernel_sha256': parent_manifest['kernel_sha256'],
        'site_permutation_seed': args.seed if args.mode == 'permute' else None,
        'control_mode': args.mode,
        'control_seed': args.seed,
        'is_ablation_control': True,
        'ablation_basis': basis,
    })
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({k: manifest[k] for k in ('kernel_sha256', 'control_mode', 'control_seed')}, indent=2))


if __name__ == '__main__':
    main()
