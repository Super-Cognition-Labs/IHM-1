# Interactive scene mechanics — the reduced-kinematics engine, which is NOT the body

**Renamed 18 September 2026.** This engine is served at
`/api/reduced-kinematics/sessions`. It was served at `/api/scene/sessions`, which
a caller could reasonably read as the body; that path is retired and returns
`410` naming `/api/embodied/sessions`, with no alias. Every response this engine
returns now carries `is_body_simulation: false`, `use_instead`, and the three
facts below at the **top level** rather than only inside `scope`:
`body_rotations` (the body cannot rotate), `body_object_contact: false` (dragged
objects pass straight through the body) and `body_environment` (no body-surface
mattress or floor contact solve). The body is
`ihm/assembly/embodied.py` (OpenSim/Simbody + BioGears) at
`/api/embodied/sessions`, and the shipped UI has always driven that one
(`app/src/scene-interaction.js:63`).

The workbench provides Select, Gimbal and Force modes. The latter two apply a
bounded cursor spring to the selected mechanical owner; the gimbal does not
teleport anatomy. Studio, floor and bed environments use separate, recorded
live sessions. Historical physiological replay pauses while a live scene owns
the display, and Reset releases that ownership.

The body uses 2,408 canonical linked translations and affine tissue regions.
Reference orientations remain constrained: off-centroid body moments are
constraint reactions, not free articulated rotations. Floor and bed use named
ideal supports about an **assumed balanced gravity preload**. Gravitational
prestress, support distribution and body-surface contact are not solved. This
is an incremental model, not validated standing, lying or walking dynamics.

The free sphere has finite mass and rotational inertia, ground contact,
inelastic normal impulse and Coulomb friction. Its force point must lie inside
or on the sphere and follows material rotation across solver substeps and
cursor commands. Body-object contact, clothing contact and feedback into
physiology are not yet coupled to this scene; the monitor states these limits.

Each session retains the exact consumed canonical mechanics bytes, import-time
source bytes and normalized loaded-code fingerprints. Changed source files
require restarting the server. Immutable compressed events are atomically
published with predecessor hashes. Failed publication rolls back mechanical
state; failed close preserves the active session. Strict integer sequences
reject duplicated advances. Following an uncertain HTTP response the client
reads current state and pauses, without replaying a potentially committed
force.

Endpoints are local-workbench-only: POST `/api/reduced-kinematics/sessions`, GET
`/api/reduced-kinematics/sessions/{id}`, and POST
`/api/reduced-kinematics/sessions/{id}/step`, `/insert` or `/close`. The old
`/api/scene/sessions` spellings of all four return `410`. GET
`/api/scene/catalog` keeps its name and is unaffected: it is the
environment/object/scene tile catalogue, it claims no physics of its own, and the
live embodied body draws its environment ids from it. It is not disclosed as
"not the body" for that reason, and it rewrites the retired session path out of
the derived catalogue's `selection` blocks on the way out
(`ihm/app/scenes.py`). Events are retained under
`data/derived/interactive-scenes/{id}/`. A server restart does not resume an
in-memory session from those records.

The actual-body force check retained at
`data/derived/audits/interactive-scene-wdowc0_3/verification.json` verifies a
0.02 N s total impulse and sphere energy/friction accounting. Later lifecycle
hardening uses a one-owner fixture rather than repeating whole-body runs:

```sh
nice -n 10 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .venv/bin/python scripts/verify_interactive_scene.py --light
nice -n 10 node --test app/test/scene-interaction.test.js
```

To share the machine with IBM-1 training, rendering is capped at 30 FPS with
pixel ratio 1, geometry loads use two concurrent requests, and hidden tabs
perform no render or simulation updates. Hair strand dynamics are opt-in;
visible synthesized strands remain available without continuous rod solves.
Source meshes and full follicle populations retain their original precision.
No background scene stepping occurs without the browser's explicit Start.
