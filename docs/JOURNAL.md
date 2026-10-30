# IHM-1 lab notebook

Dated, append-only. What was tried, with the exact command, what happened, and the dead
ends, so the next person (or agent) does not pay for them twice. Results that matter are
also written up in their own docs; this file is the trail.

## 2026-10-30: closing runs on a second host (gb10-direct)

All runs on `gb10-direct` (machine-id `96ac3ff3…`, aarch64, GB10), `nice -n 15`, at most
10 jobs, sharing the box with emulator work. IHM-1 and IBM-1 copied to the same paths as
the original host; IHM-1 `.venv` rebuilt from the lockfile
(`uv sync --frozen --extra brain --extra plot --extra data`).

**Second-host gotcha 1: IBM-1's payloads are not in the copy.** On the original host
`IBM-1/ckpt` and 375 `IBM-1/data/sources/*/raw` directories are symlinks to a portable
drive (`/media/brandonin/Portable2TB/IBM-1/...`), so an rsync of `~/Documents/IBM-1` carries
dangling links. Fixed on the peer by moving each link aside (`*.portable-link`) and copying
the real contents: all of `ckpt/` (7.7 GB; `ibm1_implicit.pt` sha256 `e6681f04…` matches) and
every `raw/` payload under 300 MB (358 of them). The 17 larger payloads (cmu-mocap, narratives,
things-images, h01, …) were not needed by anything run here and were not copied.

**Thoracic reassignability register (the WIP at 2ce7f49): reproduced exactly.**
`python -m scripts.measure_thoracic_reassignability` on the second host rewrote
`out/thoracic_reassignability.json` identically to the committed one (every field equal
once timestamps are ignored): 44 torso-attached muscles, 10 reassignable (5 per side),
verdict UNBLOCKED. The WIP commit's content stands; this entry is its resolution.

**Atlas renderer: heights were in the wrong frame.** `scripts/render_anatomy_atlas.py`
placed the transverse cuts, the trunk dissection and the head views at absolute heights
(1.30 m, 1.15 m, 1.60 m…), but the canonical anatomy frame is centred: the skin spans
y = −0.865 … +0.855 m. So the "transverse thorax" cut was above the head (nothing removed;
the panel showed a whole body from above), and the head sheet was empty. Heights are now
given above the sole and converted (`floor = min skin y`). Two more fixes from looking at
the sheets: the "organs" layer of the peel-away included the connective system (a beige
fascia sheet that hid everything under it), and the by-level panel clipped bones with
`y >= 9.0`, i.e. removed all of them. Sections now show real slices (thorax: heart and
lungs; abdomen: bowel and aorta; pelvis), and the dissection reads skin → muscle → organs →
bone.

**A process-lifetime trap.** A render started with `setsid nohup … &` inside a non-
interactive `ssh` died silently after four sheets when the ssh session ended (no traceback,
memory was fine). Long jobs now run as `systemd-run --user --unit=…` with the log as a file.

**120 s supine session (`verify_unified_world.py --environment supine --arms full`):
not reproducible on the second host; stopped after eight attempts.** The implicit IBM
controller reads IBM-1 *source* at run time (`scripts/pretrain_video_loop.py`,
`ibm/processes/cord.py`, `ibm/anatomy/muscles.py`) from `~/Documents/IBM-1`, so the result
depends on which IBM-1 commit is checked out:
- at `267ad26` (needed for the stance-policy rebuilds) and at `c5606d7` (09-09 04:07, the
  last IHM change to `ibm_controller.py`): `ValueError: antagonist must be one valid integer
  muscle index per muscle`. IHM passes `-1` for unpaired muscles; IBM-1 only accepts it from
  `a9d0adc` (09-09 05:54).
- at `c8ada7a` (09-09 10:27, after the cord fix and before the loop script's surface/atlas
  rewrite): `RuntimeError: mat1 and mat2 shapes cannot be multiplied (1x59 and 2048x256)`.
- at `c75621f` (09-21, contemporary with IHM HEAD): `(1x84 and 2048x256)`.
Along the way: the newer IBM-1 needs `nibabel` (installed into the peer `.venv` only) and
the fsaverage/Desikan payloads (see gotcha 1). The 2048-wide input is the checkpoint's; the
59/84-wide inputs are what the controller builds from the body, so no IBM-1 revision on this
host matches the `ibm1_implicit.pt` that the original host used with IHM HEAD. The honest
record: the M0_w re-base already landed in `e6b88a6` (HANDOFF §5 is stale on that), and the
supine *duration* after it is unmeasured. Reproducing it needs the exact IBM-1 working tree
the original host had at the time (it was likely uncommitted). The tested IBM-1 states are
kept as worktrees (`~/Documents/IBM-1-sep09`, `IBM-1-sep21`); a HOME shim
(`~/ihm-run/home-sep*/Documents/IBM-1 -> worktree`) selects one without moving the pinned
checkout.
