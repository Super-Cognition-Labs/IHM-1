// What the body is actually made of, this session -- as a value, with no DOM in it.
//
// Each of these is a capability the native plant has always accepted and that
// nothing serving a live body could ask for: ArticulatedBodyPlant did not forward
// the keywords, so the only callers were offline scripts. The live workbench
// therefore ran a body with no joint limits, standing on spheres inscribed in its
// inertia ellipsoids, carrying none of the derived tissue force elements, with the
// anatomy drawn in a frame that does not move it to where the body is.
// docs/WORKBENCH_AUTHENTICITY.md Tier 1 is that gap.
//
// Off is the historical body, exactly, so an unchanged control reproduces every
// earlier measurement -- which is why the note below states what each one being
// OFF means rather than congratulating the reader for turning it on.

export const CONTACT_OPTIONS = [
  { value: '', label: 'COM spheres (inertia ellipsoid)' },
  { value: 'skin', label: 'Skin surfaces' },
  { value: 'bone_all', label: 'Bone surfaces' },
];

export const POSE_OPTIONS = [
  { value: '', label: 'Force frame (does not move the anatomy)' },
  { value: 'opensim', label: 'Verified binding · OpenSim pivot' },
  { value: 'anatomical', label: 'Verified binding · anatomical pivot' },
];

/** The session configuration these controls mean. Absent keys are the historical plant. */
export function fidelityConfiguration({ contact = '', jointStops = false, tissue = false, displayPose = '' } = {}) {
  const fidelity = {};
  if (contact) fidelity.segment_contact = contact;
  if (jointStops) fidelity.joint_stops = true;
  if (tissue) fidelity.tissue_ligaments = 'admissible';
  return {
    ...(Object.keys(fidelity).length ? { mechanical_fidelity: fidelity } : {}),
    ...(displayPose ? { display_pose: displayPose } : {}),
  };
}

/** What is NOT true of this body, in the words the measurement supports. */
export function fidelityNote(state = {}) {
  const { contact = '', jointStops = false, tissue = false, displayPose = '' } = state;
  const missing = [];
  if (!jointStops) missing.push('no joint limits are enforced (the model declares ranges and nothing holds the body inside them)');
  if (!contact) missing.push('the body stands on spheres inscribed in its inertia ellipsoids, not on its skin');
  if (!tissue) missing.push('645 tissue structures exist as geometry and carry no force');
  if (!displayPose) missing.push('the anatomy is drawn in the force frame, which does not move it to where the body is');
  return missing.length
    ? 'This body: ' + missing.join('; ') + '.'
    : 'Joint stops, real surface contact, tissue forces and the verified anatomy pose are all on.';
}

/** Real segment surfaces are an upright capability; supine already has its bed foundation. */
export function contactAvailable(environment) {
  return environment === 'floor';
}
