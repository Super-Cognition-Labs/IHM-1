import test from 'node:test';
import assert from 'node:assert/strict';
import { fidelityConfiguration, fidelityNote, contactAvailable } from '../src/mechanical-fidelity.js';

test('the default controls ask for the historical plant, and send no keys for it', () => {
  // Off must be off by ABSENCE, not by an explicit null: a session that sends
  // nothing reproduces every earlier measurement to the float.
  assert.deepEqual(fidelityConfiguration(), {});
  assert.deepEqual(fidelityConfiguration({ contact: '', jointStops: false, tissue: false, displayPose: '' }), {});
});

test('each control maps to the server-owned identity the resolver understands', () => {
  assert.deepEqual(fidelityConfiguration({ contact: 'skin' }), { mechanical_fidelity: { segment_contact: 'skin' } });
  assert.deepEqual(fidelityConfiguration({ jointStops: true }), { mechanical_fidelity: { joint_stops: true } });
  // only the admissible subset is offered: the unfiltered set is measured WORSE
  // than no tissue on every drop (docs/TISSUE_MECHANICS.md).
  assert.deepEqual(fidelityConfiguration({ tissue: true }), { mechanical_fidelity: { tissue_ligaments: 'admissible' } });
  assert.deepEqual(fidelityConfiguration({ displayPose: 'opensim' }), { display_pose: 'opensim' });
});

test('the display pose travels separately from the mechanics, because it is not mechanics', () => {
  const all = fidelityConfiguration({ contact: 'skin', jointStops: true, tissue: true, displayPose: 'anatomical' });
  assert.deepEqual(all, {
    mechanical_fidelity: { segment_contact: 'skin', joint_stops: true, tissue_ligaments: 'admissible' },
    display_pose: 'anatomical',
  });
});

test('the note says what is NOT true of this body, not what is', () => {
  const bare = fidelityNote({});
  assert.match(bare, /no joint limits are enforced/);
  assert.match(bare, /spheres inscribed in its inertia ellipsoids/);
  assert.match(bare, /carry no force/);
  assert.match(bare, /does not move it to where the body is/);
  assert.ok(bare.startsWith('This body:'));
});

test('a fully specified body says so without listing absences', () => {
  // "fully specified" now includes a COUPLED soft-tissue layer: a body whose segments
  // are all rigid is still missing something, and the note must keep saying so.
  const full = fidelityNote({ contact: 'skin', jointStops: true, tissue: true, displayPose: 'opensim',
                              softTissue: 'layer_fitted_local_confined_coupled' });
  assert.doesNotMatch(full, /no joint limits/);
  assert.doesNotMatch(full, /^This body:/);
  assert.match(full, /all on/);
});

test('one control off is named, and the others are not claimed', () => {
  const note = fidelityNote({ contact: 'skin', jointStops: true, tissue: true, displayPose: '',
                              softTissue: 'layer_fitted_local_confined_coupled' });
  assert.match(note, /does not move it to where the body is/);
  assert.doesNotMatch(note, /inertia ellipsoids/);
});

test('real segment surfaces are offered only where the plant accepts them', () => {
  // supine already carries a measured skin foundation against the bed; two contact
  // models on one body would be a different claim from either.
  assert.equal(contactAvailable('floor'), true);
  assert.equal(contactAvailable('bed'), false);
  assert.equal(contactAvailable('studio'), false);
});

test('the soft tissue layer is selectable, and an uncoupled layer says it is not in the loop', () => {
  assert.deepEqual(fidelityConfiguration({ softTissue: 'layer_fitted_local_confined' }),
    { mechanical_fidelity: { soft_tissue: 'layer_fitted_local_confined' } });
  assert.deepEqual(fidelityConfiguration({ softTissue: 'layer_fitted_local_confined_coupled' }),
    { mechanical_fidelity: { soft_tissue: 'layer_fitted_local_confined_coupled' } });
  // the distinction the register spent a day establishing: a layer that deforms and a
  // layer whose reaction reaches the body are not the same claim.
  assert.match(fidelityNote({ softTissue: 'layer_fitted_local_confined' }), /NOT in the loop/);
  assert.doesNotMatch(fidelityNote({ softTissue: 'layer_fitted_local_confined_coupled' }), /NOT in the loop/);
  assert.match(fidelityNote({}), /every segment is rigid/);
});

test('soft tissue off sends no key at all, like every other control', () => {
  assert.deepEqual(fidelityConfiguration({ softTissue: '' }), {});
});
