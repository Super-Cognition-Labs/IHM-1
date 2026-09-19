"""The two bodies' defining numbers, declared once, with where each came from.

IHM-1 runs **two bodies of different stature and different mass**, and it always
has.  They are different objects and this module does not pretend otherwise:

=========================  ===================  ====================
quantity                   mechanical body      anatomical body
=========================  ===================  ====================
stature                    1.7972725296074763   1.7194712 m
mass                       77.6122029 kg        70.7713 kg
=========================  ===================  ====================

The mechanical mass was recorded across this repository as having **no
derivation**.  It has one, found by this audit on 2026-09-18 and written out
under ``MECHANICAL_TARGET_MASS_KG`` below: it is the BioGears StandardMale
*physiology* patient's weight at t = 0, 170 lb plus half a litre of stomach
water, exact in IEEE-754 double arithmetic.  It is not a mechanical number at
all, and ``scripts/verify_body_constants.py`` re-derives it on every run.

Making them one body is a **modelling decision** -- it means choosing a subject,
re-fitting the OpenSim path polynomials and re-composing the mass ledger -- and
nothing in this module does it.  What this module does is make the seam
*declared in one place* so it cannot drift: every live call site in ``ihm/``
reads a name from here, and ``scripts/verify_body_constants.py`` fails if a raw
literal of any of these numbers reappears in ``ihm/`` outside this file.

**What is deliberately NOT unified.**  ``data/derived/**`` artefacts, ``out/``,
``logs/`` and the dated records in ``docs/`` carry these numbers as *recorded
history*: they say what a particular run was handed, and rewriting them would
falsify the record.  They are audited (``docs/BODY_PARAMETERS.md``, section
"The two-body constants, counted") and left alone.

**Two mirrors that cannot import this module** are declared in ``MIRRORS``
below and checked by the guard instead of being rewritten.  A mirror that is
checked cannot drift either; it just cannot be an import.

This module imports nothing.  That is on purpose: it must be readable by a
static parser and importable from anywhere, and a constant with a dependency is
a constant with a failure mode.
"""

# ---------------------------------------------------------------------------
# The mechanical body -- the OpenSim Rajagopal subject that ``ihm.native``
# integrates.  22 rigid bodies, 80 muscles, 33 coordinates.  This is the body
# that touches the world.
# ---------------------------------------------------------------------------

#: Path, relative to the repository root, of the model every mechanical number
#: below is measured on.
MECHANICAL_MODEL = ('data/raw/mechanics/opensim-core/OpenSim/Examples/Moco/'
                    'example3DWalking/subject_walk_scaled.osim')

#: Sum of the 22 ``<Body><mass>`` entries in that ``.osim``.
#:
#: DERIVED, and recomputable: ``ihm.body_parameters.assert_measured_defaults``
#: re-adds them from the file and raises if this drifts by more than 1e-9.
#: The native engine divides it into the requested ``target_mass_kg`` to get its
#: uniform ``mass_scale`` (``scripts/native_mechanical_stream.cpp`` line 51).
MECHANICAL_SOURCE_MASS_KG = 85.26984854173146

#: The mass every mechanical run actually asks for.
#:
#: **DERIVED, and the derivation is NOT mechanical.  Found 2026-09-18; this
#: constant was recorded as "no derivation found" in IHM-1 ``CLAUDE.md``, in
#: ``ihm/body_parameters.py`` and in ``docs/WORKBENCH_AUTHENTICITY.md`` until
#: then, and every one of those statements was wrong.**
#:
#: It is the **BioGears StandardMale patient's weight at t = 0**::
#:
#:     170 lb x 0.45359237 kg/lb      = 77.1107029 kg   the declared patient
#:     + 0.5   L  water               =  0.5      kg  \
#:     + 500   mg calcium             =  0.0005   kg   |  StomachContents at t=0
#:     + 1     g  sodium              =  0.001    kg  /
#:     -------------------------------------------------
#:                                    = 77.6122029 kg
#:
#: That sum equals this constant in **exact IEEE-754 double arithmetic**, not to
#: a tolerance.  The number also appears verbatim in this repository, as
#: ``<Weight unit="kg" value="77.6122029"/>`` (line 8) and
#: ``<RestingPatientMass_kg>`` (line 1235) of
#: ``data/raw/physiology/biogears/share/data/states/StandardMale@0s.xml``, and in
#: 7 further male ``*@0s`` state files -- 8 of the 19 BioGears initial states.
#:
#: **What that means, and it is a modelling fact rather than a tidy-up.**  The
#: mechanical plant is scaled to the weight of a 170 lb BioGears *physiology*
#: reference patient carrying half a litre of stomach water.  It has nothing to
#: do with the OpenSim subject, and it is the SAME patient whose 170 lb this
#: repository's own composition ledger rejected as unreachable in the acquired
#: envelope (see ``SUPERSEDED_BIOGEARS_STANDARDMALE_MASS_KG``).  The mechanics
#: therefore run at the physiology patient's fed weight while the anatomy runs at
#: a mass composed from its own geometry, and that is the mass half of the
#: two-body seam, now with a name.
#:
#: The four candidate derivations checked before this one are all still refuted,
#: and are kept because a future reader will re-check them:
#:
#: * not the source model's summed body masses (85.26984854173146 kg);
#: * not the anatomical body's composed mass (70.7713 kg);
#: * not the superseded BioGears StandardMale patient weight alone
#:   (77.1107029 kg -- it is that PLUS the stomach);
#: * not recoverable from the trial's ground reaction force -- the mean vertical
#:   GRF over ``grf_walk.mot`` gives 60.83 kg, which is a duty-cycle artefact of
#:   averaging over swing as well as stance, not a weight.
#:
#: The superseded hypothesis, written down so it is not re-proposed: "most likely
#: the AddBiomechanics subject's measured mass".  It is refuted -- the identity
#: above is exact and the AddBiomechanics subject has no such number.
#:
#: **What is established and what is inferred.**  The arithmetic identity and the
#: eight verbatim occurrences are established.  The *intent* -- that the literal
#: was lifted so the mechanics and the physiology would agree on one body weight
#: -- is inference, from the commit that introduced it (``756078c``, "Integrate
#: continuing native articulation, metabolic states and canonical feedback",
#: 5 Sep 2026, as ``target_mass_kg=77.6122029`` in an acceptance script written
#: beside the metabolic-state work).  Nothing states it.
#:
#: What it IS legitimate to use this for: the default ``target_mass_kg`` of a
#: mechanical run.  It is not the anatomical body's mass and it does not describe
#: the anatomy's composition.  It is not "the subject's mass" either: quote it as
#: the BioGears reference patient's fed weight, which is what it is.
MECHANICAL_TARGET_MASS_KG = 77.6122029

#: The pieces of the identity above, kept separate so ``scripts/
#: verify_body_constants.py`` can re-derive the constant instead of trusting the
#: prose.  Units are declared because the derivation crosses two of them: the
#: patient weight is an imperial round number and the stomach is metric.
BIOGEARS_STANDARDMALE_WEIGHT_LB = 170.0
POUND_KG = 0.45359237          # international avoirdupois pound, exact by definition
BIOGEARS_INITIAL_STOMACH_KG = 0.5015   # 0.5 L water + 500 mg calcium + 1 g sodium
BIOGEARS_STANDARDMALE_STATE = ('data/raw/physiology/biogears/share/data/states/'
                               'StandardMale@0s.xml')

#: Vertical distance, at the model's default pose, from the plane of the
#: AddBiomechanics ``*Ground`` virtual markers (R/L ``HeelGround``, ``MT5Ground``,
#: ``ToeGround``) to the ``Head`` marker.
#:
#: MEASURED, and recomputable: ``ihm.native.model_scaling.head_marker_height_m``
#: reads it out of the ``.osim`` and ``assert_measured_defaults`` gates it at
#: 1e-9.  It is a stature PROXY and not a standing height: the ``Head`` marker
#: sits on the head, not at the vertex.
#:
#: Legitimate use: the denominator of ``stature_scale`` for anything measured on
#: the MECHANICAL body -- segment lengths, muscle paths, contact geometry, the
#: displayed anatomy's placement transform.  It is **not** the denominator for
#: anything measured on the anatomical skin; see ``ANATOMICAL_STATURE_M``.
MECHANICAL_STATURE_M = 1.7972725296074763


# ---------------------------------------------------------------------------
# The anatomical body -- the BodyParts3D atlas plus BioGears.  3,816 entities,
# posed kinematically from the mechanical segments.  This body feels nothing.
# ---------------------------------------------------------------------------

#: Path, relative to the repository root, of the profile the anatomical numbers
#: below are gated against.
ANATOMICAL_PROFILE = 'data/derived/canonical/profile.json'

#: Head-to-toe extent of the canonical skin mesh (BodyParts3D ``FJ2810``), in
#: metres, as ``ihm/assembly/profile.py`` writes it into ``profile.json``.
#:
#: MEASURED, and recomputable: ``assert_measured_defaults`` reads
#: ``profile.json['height_m']`` and gates it at 1e-7.
#:
#: Legitimate use: the denominator of the scale factor for anything MEASURED ON
#: THE SKIN OR THE ATLAS -- nerve route lengths, skin-patch route lengths,
#: garment registration.  Using ``MECHANICAL_STATURE_M`` for those is the bug
#: that gave a "2.03 m" variant the nerves of a 1.942 m body
#: (``docs/BODY_PERIPHERAL.md``, 18 Sep 2026).
ANATOMICAL_STATURE_M = 1.7194712

#: Total body mass composed over this specimen's own measured interior: an 8 mm
#: voxel partition of the acquired envelope filled at sourced per-constituent
#: densities.
#:
#: DERIVED AND GATED: ``ihm/assembly/profile.py`` refuses to build a profile
#: unless ``data/derived/interstitial-composition-prior-v1/ledger.json`` still
#: composes to this within 5e-5 kg, and ``assert_measured_defaults`` then gates
#: ``profile.json`` against it at 1e-7.  Sensitivity to the two free choices in
#: the ledger is 0.024 kg (8 mm vs 10 mm partition) and 0.534 kg (residual
#: density swept over [1000, 1060] kg/m3).
#:
#: Legitimate use: statements about the anatomy's composition, and the
#: physiology patient record.  It is **not** what the mechanical plant is scaled
#: to and never has been.
ANATOMICAL_MASS_KG = 70.7713

#: The BioGears StandardMale patient's DECLARED weight: 170 lb, exactly, which
#: the patient file writes in pounds and every state file writes in kilograms.
#: This repository replaced it as the anatomical mass.
#:
#: HISTORICAL as an anatomical mass, and NOT unrelated to the mechanical one:
#: ``MECHANICAL_TARGET_MASS_KG`` is this number plus the 0.5015 kg of stomach
#: contents BioGears seeds at t = 0.  They are the same patient, fasted and fed.
#: It resembles the mechanical mass to two figures because it *is* the mechanical
#: mass, less a stomach -- which is exactly the kind of near-coincidence that
#: gets one number quoted for the other, and is why both are named here.
#:
#: It is not reachable in the acquired envelope: it would need the interstitial
#: void to weigh 1185.8 kg/m3, denser than every ICRU-44 soft tissue, which by
#: the Siri relation is a negative fat fraction.  **Never use it as a mass for
#: the anatomical body.**
SUPERSEDED_BIOGEARS_STANDARDMALE_MASS_KG = 77.1107029


# ---------------------------------------------------------------------------
# The seam between them
# ---------------------------------------------------------------------------

#: Fitted similarity scale taking OpenSim body geometry onto the BodyParts3D
#: atlas (``scripts/bind_anatomy_to_segments.py``;
#: ``docs/ANATOMY_SEGMENT_BINDING.md``).  This is the only number that relates
#: the two bodies' sizes, and it is a FIT, not an identity.
ANATOMY_REGISTRATION_SCALE = 0.963


def _rel(x, y):
    return abs(x - y) / y


#: How far apart the two bodies are, as a fraction.  Printed rather than
#: absorbed, because it is the honest error bar on any statement of the form
#: "the body is N metres tall" or "the body weighs N kilograms".
STATURE_DISAGREEMENT = _rel(MECHANICAL_STATURE_M * ANATOMY_REGISTRATION_SCALE,
                            ANATOMICAL_STATURE_M)
MASS_DISAGREEMENT = _rel(MECHANICAL_TARGET_MASS_KG, ANATOMICAL_MASS_KG)

#: The seam, in one sentence, for anything that has to disclose it.
KNOWN_SEAM = (
    'Two bodies. Nerve and skin-patch route lengths scale by '
    'stature / ANATOMICAL_STATURE_M (%.7f m, the atlas skin extent); the '
    'displayed anatomy and the mechanical plant scale by '
    'stature / MECHANICAL_STATURE_M (%.7f m, the OpenSim head-marker height). '
    'They disagree by %.2f%% in stature and %.1f%% in mass. Unifying them is a '
    'modelling decision, not a refactor.'
    % (ANATOMICAL_STATURE_M, MECHANICAL_STATURE_M,
       100 * STATURE_DISAGREEMENT, 100 * MASS_DISAGREEMENT))


# ---------------------------------------------------------------------------
# Mirrors: copies that exist for a reason and are checked instead of removed
# ---------------------------------------------------------------------------

#: ``(path, assignment name, constant name, why it is not an import)``.
#:
#: ``scripts/verify_body_constants.py`` parses each file with ``ast`` -- it does
#: NOT import them, because two of the three cannot be imported where they run
#: -- and fails if the mirrored value stops matching the constant here.  A
#: checked mirror cannot drift; it just is not an import.
MIRRORS = (
    ('ihm/assembly/garment_wardrobe.py', 'REFERENCE_HEIGHT_M',
     'ANATOMICAL_STATURE_M',
     'That module declares in its own docstring that it imports numpy and scipy '
     'ONLY, because scripts/fit_garments_to_envelope.py loads it out of process '
     'inside the vendored libigl venv (data/runtime/geometry/libigl-2.6.2). '
     'An import was checked and does work there today -- "import ihm" pulls '
     'numpy and scipy and nothing else as of 2026-09-18 -- but taking it would '
     'silently make garment fitting depend on everything ihm/__init__.py ever '
     'grows. A checked mirror keeps the declared import surface and still '
     'cannot drift.'),
    ('ihm/assembly/anatomy_pose.py', 'REGISTRATION_SCALE',
     'ANATOMY_REGISTRATION_SCALE',
     'anatomy_pose.py takes ihm.body_parameters only as a deferred import '
     'inside a function, deliberately, to keep module import free of the '
     'parameter layer; it is the module verified to 7.8e-16 against Simbody '
     'and is not edited for tidiness. Its docstring also quotes the rounded '
     '1.7195 m and 0.963, which is why the text scan skips mirror files: the '
     'value that could drift is pinned by the assignment check above.'),
)

#: Files outside ``ihm/`` that are known to carry a raw literal and are NOT
#: category (c).  Listed so the audit in ``docs/BODY_PARAMETERS.md`` can be
#: reproduced, not so the guard can skip them -- the guard only ever looks
#: inside ``ihm/``.
RECORDED_HISTORY_NOTE = (
    'data/derived/**, data/research/**, out/, logs/ and dated entries in docs/ '
    'carry these numbers as a record of what a particular run was handed. '
    'data/research/hair_*.json additionally carries 77.61218724193209 kg, '
    'which is MECHANICAL_TARGET_MASS_KG minus 1.5658e-05 kg of separately '
    'partitioned hair -- a derived residual, not a transcription error. None of '
    'it is rewritten.')
