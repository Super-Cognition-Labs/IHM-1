"""Server-owned mechanical fidelity options for a live body.

WHY THIS EXISTS.  `NativeMechanicalStream` has accepted joint stops, real segment
contact surfaces and derived tissue force elements for some time.  Nothing that
serves a live body could ask for them: `ArticulatedBodyPlant.__init__` did not
forward the keywords, so the only callers were offline `scripts/*.py`.  The live
workbench therefore ran a body with **no joint limits**, standing on **spheres
inscribed in inertia ellipsoids**, carrying **none** of the 117 tissue elements --
while the replacements sat built, measured and documented on disk.
`docs/WORKBENCH_AUTHENTICITY.md` Tier 1 is that gap; this module is its fix.

WHAT IT IS.  A resolver, in the same shape as `controller_selection.py`: the
client names an identity, never a path, and gets back the resolved parameters
plus a disclosure it cannot render a result without having been handed.  Anything
not named is off, and off is the historical behaviour, so an existing caller is
unaffected to the float.

WHAT IT IS NOT.  It does not make any of these true by default.  Each option
below is an improvement with a measured cost, and two of them are actively worse
in some configurations -- the selections carry those numbers rather than hiding
them.
"""
from pathlib import Path
import json
import xml.etree.ElementTree as ET

MODEL = 'data/models/engineering_stance_v1/model.osim'

# A coordinate declaring +-10 rad is not declaring a range, it is declaring the
# absence of one; the six shoulder coordinates do exactly that.  Stopping them at
# a number the model never gave would be inventing a limit, so they are left free
# and `unranged_coordinates` says which they are on every resolution.
_RANGE_SENTINEL_RAD = 18.0
_TRANSLATIONS = ('pelvis_tx', 'pelvis_ty', 'pelvis_tz')

# Measured, not chosen.  scripts/crawl.py swept stiffness at IDENTICAL port gains
# over 3 s of the seed pattern (its own note records that the first sweep compared
# a no-stop row against stopped rows at different gains, which is this programme's
# recurring "compared against the wrong thing"):
#
#   k (N.m/rad)   s / advance   worst excursion past the declared range
#   none              0.515     1.450 rad   (ankle_angle_r -- the folded foot)
#    300              0.937     0.139 rad
#    100              0.505     0.192 rad
#     30              0.264     0.220 rad
#
# 30 is HALF the wall clock of no stop at all and holds the plant 6.6x closer to
# its declared range, because a plant kept out of absurd configurations is a plant
# the error controller can integrate.  These are engineering constants stated by
# the caller, NOT measured ligament properties, and the resolution says so.
#
# DAMPING, and a defect that is now fixed. Until 18 Sep 2026 the engine
# (scripts/native_mechanical_stream.cpp) converted a stop's limits, stiffness and
# transition from radians to degrees for CoordinateLimitForce and passed damping
# THROUGH, while OpenSim declares that property Nm/(degree/s). A declared 1.5 was
# applied as 1.5 * 180/pi = 85.94 N.m.s/rad. The engine now converts it, and the
# constant is re-declared at 85.94 -- the value every stopped run in this repo was
# actually measured at -- so the fix moves no plant: a 50-step stopped trajectory
# reproduces bit-for-bit across the fix.
#
# The profiles are the stiffnesses crawl.py actually swept, at the damping and
# transition that sweep held fixed. An earlier version of this module gave `firm`
# and `stiff` damping 2.0/3.0 and transitions 0.25/0.20 that no measurement
# supported; they were typed, and are withdrawn.
_SWEPT_DAMPING = 1.5 * 180.0 / 3.141592653589793       # 85.94 N.m.s/rad, as measured
_SWEPT_TRANSITION = 0.35
JOINT_STOP_PROFILES = {
    'measured_soft': {'stiffness_nm_per_rad': 30.0, 'damping_nm_s_per_rad': _SWEPT_DAMPING,
                      'transition_rad': _SWEPT_TRANSITION},
    'firm': {'stiffness_nm_per_rad': 100.0, 'damping_nm_s_per_rad': _SWEPT_DAMPING,
             'transition_rad': _SWEPT_TRANSITION},
    'stiff': {'stiffness_nm_per_rad': 300.0, 'damping_nm_s_per_rad': _SWEPT_DAMPING,
              'transition_rad': _SWEPT_TRANSITION},
}
JOINT_STOP_DAMPING_NOTE = (
    'Damping 85.94 N.m.s/rad is the value every stopped run in this repo was measured at. '
    'It was historically declared as 1.5 because the engine passed damping to OpenSim '
    'CoordinateLimitForce (which reads Nm/(degree/s)) unconverted; the conversion is '
    'fixed and the constant re-declared so no plant moved. It is roughly 14x critical '
    'damping for a limb segment at k=30 N.m/rad: an engineering constant, not a '
    'measured ligament property.')
DEFAULT_JOINT_STOP_PROFILE = 'measured_soft'



# The real segment surfaces, in place of the upright environment's COM spheres.
# `skin` is the one the programme is actually about: docs/ACTUATION_STAGES.md says
# contact with the world is never bone against world, and this is the bundle where
# the skin is what meets the floor.
# `replaces_source_feet` is the whole question, not a detail.
#
# `scripts/measure_segment_contact_meshes.py` records, and EXCLUDES the arm from its
# own default run for it: **the skin never reaches the floor in the stance pose.**
# The skin surface of the foot sits above the source foot contact spheres' effective
# plane, so a bundle that REPLACES those spheres leaves a standing body with nothing
# under it, and it collapses. That is not the skin failing to hold the body; it is a
# pose in which the skin is not yet touching.
#
# So both arms are offered and neither is called the truth: `skin` replaces the feet
# (the intended end state, and the one that collapses from the stance pose), and
# `skin_carried` keeps them (the cost of CARRYING the geometry, separated from the
# cost of a plant that is collapsing -- two things one number would mix).
SEGMENT_CONTACT_BUNDLES = {
    'skin': {
        'path': 'data/derived/segment-contact-meshes/skin',
        'layer': 'skin',
        'label': 'Skin exterior, replacing the source feet',
        'replaces_source_feet': True,
        'caveat': 'The skin does not reach the floor in the stance pose. With the source '
                  'foot spheres replaced, a standing body has nothing under it and collapses. '
                  'Use skin_carried to separate the cost of the geometry from that.',
    },
    'skin_carried': {
        'path': 'data/derived/segment-contact-meshes/skin',
        'layer': 'skin',
        'label': 'Skin exterior, source feet kept',
        'replaces_source_feet': False,
        'caveat': 'The source foot spheres still carry the body; the skin is present and '
                  'only loads where it actually touches. This measures carrying the '
                  'geometry, NOT the skin holding the body up.',
    },
    'skin_layer_map': {
        'path': 'data/derived/segment-contact-meshes/skin-layer-map-v1',
        'layer': 'skin',
        'label': 'Skin exterior with per-patch measured depth, replacing the source feet',
        'replaces_source_feet': True,
        'caveat': 'Same caveat as skin: it does not reach the floor in the stance pose.',
    },
    'bone_all': {
        'path': 'data/derived/segment-contact-meshes/stance-bone-all',
        'layer': 'bone',
        'label': 'Every bone surface',
        'replaces_source_feet': False,
    },
    'bone_proxy': {
        'path': 'data/derived/segment-contact-meshes/stance-bone-proxy',
        'layer': 'bone',
        'label': 'Bone surfaces, proxy subset',
        'replaces_source_feet': False,
    },
}

# docs/TISSUE_MECHANICS.md measured both sets on three drops.  The full set makes
# the plant WORSE; the 66 elements that survive the kinematic check make it no
# worse and sometimes much better, and added to the stops improve the worst
# excursion past the declared ranges on every drop tested by 4-24%.  51 of the 117
# fail that check -- the derived ACL reads 77% strain at 90 degrees of knee flexion,
# which is an attachment in the wrong place.  So `admissible` is the only selection
# offered by default, and `all` exists to reproduce the negative result.
TISSUE_BUNDLES = {
    'admissible': {'path': 'data/derived/tissue-force-elements-v1',
                   'classes': ['ligament', 'joint_capsule'], 'admissible_only': True,
                   'label': 'Ligaments and capsules, kinematically admissible only'},
    'all': {'path': 'data/derived/tissue-force-elements-v1',
            'classes': ['ligament', 'joint_capsule'], 'admissible_only': False,
            'label': 'Ligaments and capsules, unfiltered (measured WORSE than none)'},
}

# A DEFORMABLE soft-tissue layer over each segment (ihm/assembly/soft_tissue_layer.py,
# docs/SOFT_BODY.md): 3-D compressible neo-Hookean tetrahedra between the skin and this
# body's own measured soft-tissue depth, carried by the segment, deforming under
# contact and returning the force and moment it transmits.  Unlike every other entry in
# this file it is NOT a NativeMechanicalStream keyword: the native plant does not
# integrate it.  Resolving it puts the layer's identity and disclosure on the body, and
# `soft_tissue_layer.build_selected_layers` builds it for a caller that will pose it.
#
# The two identities differ only in how the in-vivo APPARENT layer modulus is read as a
# 3-D material, which the measurement does not fix (soft_tissue_layer.LAYER_MODULUS_MAPPINGS).
SOFT_TISSUE_LAYERS = {
    'layer_map_confined': {
        'bundle': 'data/derived/segment-contact-meshes/skin-layer-map-v1',
        'mapping': 'confined', 'spacing_m': 0.005,
        'label': 'Deformable neo-Hookean layer at measured depth; modulus read as confined (softest)'},
    'layer_map_unconfined': {
        'bundle': 'data/derived/segment-contact-meshes/skin-layer-map-v1',
        'mapping': 'unconfined', 'spacing_m': 0.005,
        'label': 'Deformable neo-Hookean layer at measured depth; modulus read as Young\'s (stiffest)'},
}
# Measured, and re-derived by scripts/verify_soft_tissue.py, which fails if it moves: at
# 5 mm these segments' skin patches are thinner than their own measured depth
# everywhere, so no cell is core and nothing carries the layer.
SOFT_TISSUE_UNANCHORED = ('patella_l', 'patella_r')
# Which world axis the support pushes along, per environment, in the plant's own
# conventions: the supine foundation measures penetration as `plane_x - x`
# (supine_contact.foundation); the upright ground is OpenSim's y-up floor.
SOFT_TISSUE_SUPPORT = {
    'supine': {'axis': 0, 'sign': 1.0, 'basis': 'supine_contact.foundation: penetration = plane_x - x'},
    'upright': {'axis': 1, 'sign': 1.0, 'basis': 'OpenSim ground frame, y up'},
}

DISCLOSURE = {
    'joint_stops':
        'Joint stops at the coordinate ranges the source model already declares. '
        'Nothing else in this plant enforces them: every rotational coordinate '
        'carries <clamped>true</clamped>, the model holds zero CoordinateLimitForce, '
        'and OpenSim does not clamp during forward dynamics. The LIMIT is the '
        "model's own; the stiffness, damping and transition width are explicit "
        'engineering constants, not measured ligament properties. The declared range '
        "and the model's own passive stops DISAGREE -- a resting prone body already "
        'sits 0.240 rad outside its declared hip_rotation_l -- so these stops are '
        'soft by design and a coordinate may still be found outside its range.',
    'segment_contact_skin':
        'The body stands on its SKIN: the canonical exterior skin surface, cut per '
        'segment and capped, as OpenSim ContactMesh over SimTK TriangleMesh carried '
        'by ElasticFoundationForce -- an independent spring at every triangle centroid '
        'below the plane. No convex hull is taken. The skin is carried RIGIDLY by its '
        'segment: no in-plane stretch, no sliding, no deformable continuum. Every '
        'segment boundary is a seam the real body does not have, and both tali are '
        'refused for having fewer than 64 exterior triangles.',
    'segment_contact_bone':
        'The body stands on BONE surfaces. This is a collider against other bones and '
        'its own soft tissue, NOT what meets the floor in a real body; select the skin '
        'bundle for that. The layer travels into the native emit so no report can say '
        'the body stood on its skin about a run that stood on its femurs.',
    'tissue_ligaments':
        'Blankevoort1991Ligament elements over attachments derived from each '
        "structure's OWN surface. A CONSTRUCTION from mesh geometry and a published "
        'cadaver modulus -- not measured insertion footprints, not a subject-specific '
        'ligament property. These forces are internal: they can change how the plant '
        'moves and cannot change its momentum balance. They do NOT replace the joint '
        'stops.',
    'soft_tissue_layer':
        'A DEFORMABLE soft-tissue layer per segment: compressible neo-Hookean tetrahedra from '
        'the skin inward to this body\'s measured soft-tissue depth, the deeper core carried '
        'rigidly by the segment. It deforms under contact against a frictionless half-space '
        'and returns the force and moment it transmits to its segment. It is NOT in the native '
        'plant: the engine still integrates the rigid scaffold and its own contact, and nothing '
        'here feeds back into that integration unless a caller applies the returned load. The '
        'core is depth-defined, not bone (skin and bone are not co-registered in this bundle); '
        'muscle is inside the rigid core; segments are independent, with a seam at every '
        'boundary; the surface is voxelised at the cell size. Per-solve cost is seconds, not '
        'milliseconds: see docs/SOFT_BODY.md before calling it real-time.',
}


def declared_ranges(root, model=MODEL):
    """Rotational coordinate ranges as the source model declares them.

    Read out of the model rather than typed. CLAUDE.md records a day lost to a knee
    whose range was mirror-imaged between two models, and a number nobody
    transcribed cannot be transcribed wrong.
    """
    root = Path(root)
    ranged, unranged = {}, []
    for coordinate in ET.parse(root / model).getroot().iter('Coordinate'):
        name = coordinate.get('name')
        element = coordinate.find('range')
        if not name or element is None or element.text is None:
            continue
        low, high = (float(v) for v in element.text.split())
        if name in _TRANSLATIONS:
            continue
        if high - low > _RANGE_SENTINEL_RAD:
            unranged.append(name)
            continue
        ranged[name] = (low, high)
    return ranged, sorted(unranged)


def joint_stops(root, profile=DEFAULT_JOINT_STOP_PROFILE, model=MODEL):
    """The plant's coordinate-limit rows, one per coordinate that declares a range."""
    if profile not in JOINT_STOP_PROFILES:
        raise ValueError('Unknown joint stop profile')
    constants = JOINT_STOP_PROFILES[profile]
    ranged, unranged = declared_ranges(root, model)
    if not ranged:
        raise ValueError('Model declares no usable rotational ranges')
    rows = [dict(coordinate=name, lower_rad=low, upper_rad=high, **constants)
            for name, (low, high) in sorted(ranged.items())]
    return rows, unranged


def _verified_bundle(root, relative, schema):
    path = (Path(root) / relative).resolve()
    if not path.is_relative_to(Path(root).resolve()):
        raise ValueError('Owned bundle required')
    manifest = json.loads((path / 'manifest.json').read_bytes())
    if manifest.get('schema') != schema:
        raise ValueError('Bundle schema is not the one this resolver understands: ' + str(relative))
    return manifest


def resolve_fidelity(root, value=None, *, environment='supine'):
    """Resolve a client's named mechanical fidelity into plant keywords.

    Returns `(kwargs, selection)`. `kwargs` go straight to NativeMechanicalStream;
    `selection` is the disclosure, and is carried on the session so no caller can
    read a result without having been handed what it rests on.

    `None` is the historical plant, exactly: no stops, COM-sphere contact, no
    tissue. That default is deliberate -- an existing measurement must not change
    because this resolver arrived.
    """
    if value is None:
        return {}, {'joint_stops': None, 'segment_contact': None, 'tissue_ligaments': None,
                    'soft_tissue': None,
                    'basis': 'Historical plant: no coordinate limits, inertia-ellipsoid COM '
                             'sphere contact, no tissue force elements.'}
    if not isinstance(value, dict) or set(value) - {'joint_stops', 'segment_contact', 'tissue_ligaments',
                                                    'soft_tissue'}:
        raise ValueError('Unknown mechanical fidelity configuration')

    kwargs, selection = {}, {}

    stops = value.get('joint_stops')
    if stops in (None, False):
        selection['joint_stops'] = None
    else:
        profile = DEFAULT_JOINT_STOP_PROFILE if stops is True else stops
        if not isinstance(profile, str):
            raise ValueError('Joint stop selection must be a profile name or a boolean')
        rows, unranged = joint_stops(root, profile)
        kwargs['coordinate_limits'] = rows
        selection['joint_stops'] = {
            'profile': profile, 'constants': JOINT_STOP_PROFILES[profile],
            'damping_note': JOINT_STOP_DAMPING_NOTE,
            'swept': profile == 'measured_soft' or 'stiffness swept in scripts/crawl.py at this '
                     'damping and transition; only measured_soft was adopted from that sweep',
            'coordinates': [r['coordinate'] for r in rows],
            'unranged_coordinates': unranged,
            'unranged_basis': 'Declared +-10 rad, which is the absence of a range; left free '
                              'rather than stopped at a limit the model never gave.',
            'disclosure': DISCLOSURE['joint_stops']}

    contact = value.get('segment_contact')
    if contact is None:
        selection['segment_contact'] = None
    else:
        if contact not in SEGMENT_CONTACT_BUNDLES:
            raise ValueError('Unknown segment contact bundle')
        if environment != 'upright':
            raise ValueError('Segment contact meshes require the upright environment')
        spec = SEGMENT_CONTACT_BUNDLES[contact]
        manifest = _verified_bundle(root, spec['path'], 'ihm.segment-contact-meshes.v1')
        if manifest.get('layer') != spec['layer']:
            raise ValueError('Segment contact bundle layer does not match its selection')
        kwargs['segment_contact_meshes'] = spec['path']
        kwargs['segment_contact_replaces_source_feet'] = spec['replaces_source_feet']
        # E, Poisson ratio and thickness come from the canonical skin-layer entities the
        # bundle itself recorded, so the arm measures what THIS body's declared skin does
        # rather than what a default number does. A bundle carrying a per-segment layer
        # map declares its own per record, and the plant refuses a caller override there.
        material = manifest.get('skin_material') or {}
        declared = {k: material[k] for k in ('youngs_modulus_pa', 'poissons_ratio', 'layer_thickness_m')
                    if k in material}
        if declared and not manifest.get('per_record_material'):
            kwargs['segment_contact_material'] = declared
        selection['segment_contact'] = {
            'id': contact, 'label': spec['label'], 'layer': spec['layer'],
            'replaces_source_feet': spec['replaces_source_feet'],
            'meshes': manifest.get('meshes'), 'total_faces': manifest.get('total_faces'),
            'watertight_meshes': manifest.get('watertight_meshes'),
            'concave_meshes': manifest.get('concave_meshes'),
            'refused': manifest.get('refused') or [],
            'bundle_basis': manifest.get('basis'),
            'material': declared or None,
            'material_basis': 'Declared by the bundle from this body\'s own canonical skin-layer '
                              'entities; not a tuned contact stiffness.' if declared else None,
            'caveat': spec.get('caveat'),
            'disclosure': DISCLOSURE['segment_contact_skin' if spec['layer'] == 'skin'
                                     else 'segment_contact_bone']}

    tissue = value.get('tissue_ligaments')
    if tissue is None:
        selection['tissue_ligaments'] = None
    else:
        if tissue not in TISSUE_BUNDLES:
            raise ValueError('Unknown tissue ligament bundle')
        spec = TISSUE_BUNDLES[tissue]
        manifest = _verified_bundle(root, spec['path'], 'ihm.tissue-force-elements.v1')
        kwargs['tissue_ligaments'] = spec['path']
        kwargs['tissue_ligament_classes'] = list(spec['classes'])
        kwargs['tissue_ligament_admissible_only'] = spec['admissible_only']
        selection['tissue_ligaments'] = {
            'id': tissue, 'label': spec['label'], 'classes': list(spec['classes']),
            'admissible_only': spec['admissible_only'],
            'bundle_scope': manifest.get('scope'),
            'measured_cost': None if spec['admissible_only'] else
                'docs/TISSUE_MECHANICS.md measured this unfiltered set as WORSE than no '
                'tissue on every drop tested. It is offered to reproduce that result.',
            'disclosure': DISCLOSURE['tissue_ligaments']}

    soft = value.get('soft_tissue')
    if soft is None:
        selection['soft_tissue'] = None
    else:
        if soft not in SOFT_TISSUE_LAYERS:
            raise ValueError('Unknown soft tissue layer')
        if environment not in SOFT_TISSUE_SUPPORT:
            raise ValueError('A soft tissue layer needs a support to press against; '
                             'the ' + str(environment) + ' environment has none')
        from .soft_tissue_layer import BUNDLE_SCHEMA, LAYER_MODULUS_MAPPINGS
        spec = SOFT_TISSUE_LAYERS[soft]
        manifest = _verified_bundle(root, spec['bundle'], BUNDLE_SCHEMA)
        if manifest.get('layer') != 'skin' or 'layer_map' not in manifest:
            raise ValueError('Soft tissue bundle carries no measured per-segment depth')
        bodies = [r['body'] for r in manifest['records']]
        # NOTHING is added to kwargs: NativeMechanicalStream has no such keyword, and the
        # plant it builds is exactly the plant it would have built without this key.
        selection['soft_tissue'] = {
            'id': soft, 'label': spec['label'], 'bundle': spec['bundle'],
            'mapping': spec['mapping'], 'mapping_basis': LAYER_MODULUS_MAPPINGS[spec['mapping']],
            'spacing_m': spec['spacing_m'],
            'segments': [b for b in bodies if b not in SOFT_TISSUE_UNANCHORED],
            'unanchored_segments': list(SOFT_TISSUE_UNANCHORED),
            'unanchored_basis': 'skin patch thinner than its own measured depth everywhere at '
                                'this spacing: no rigid core, nothing to carry the layer',
            'refused_by_bundle': manifest.get('refused') or [],
            'support': SOFT_TISSUE_SUPPORT[environment],
            'in_native_plant': False,
            'builder': 'ihm.assembly.soft_tissue_layer.build_selected_layers',
            'disclosure': DISCLOSURE['soft_tissue_layer']}

    selection['basis'] = ('Server-owned bundles resolved by identity; a client never supplies '
                          'a path, a mesh or a material.')
    return kwargs, selection
