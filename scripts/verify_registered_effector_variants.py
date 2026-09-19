"""The fail-closed effector catalog must read EVERY registered body, and only those.

`whole_body_effector_catalog` accepted `ihm.upperbody-registration.v1` only, so it
could load the 92-muscle default plant and not the 98-muscle variant that every
stance controller forces (`data/models/engineering_stance_v1`). `peripheral_coverage`
inherited the hole: the trunk effectors of the body the live stance path actually
integrates had no port audit at all.

These tests pin both halves of the fix:

  READS   all four registered manifests on disk load, with the muscle count the
          manifest declares, and produce a peripheral coverage report.
  REFUSES the two fields the variant schema adds -- `insert_path`/`insert_sha256`
          and `base_model_path` -- are actually verified, not merely read. Each
          mutation below is a control that CAN fail: `base_model_path` is swapped
          for a real file that exists and parses and is simply not one of this
          manifest's verified sources, so a loader that only checked existence
          would pass it.

    OPENBLAS_NUM_THREADS=1 prlimit --as=2147483648 -- nice -n 10 \
      .venv/bin/python -m scripts.verify_registered_effector_variants
"""
import json
import unittest
from copy import deepcopy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

REGISTERED = {
    'data/derived/mechanics/whole_body_arm26_v2/registration.json': ('ihm.upperbody-registration.v1', 92),
    'data/derived/mechanics/whole_body_lumbar_current/registration.json': ('ihm.lumbar-muscle-variant.v1', 98),
    'data/models/engineering_stance_v1/registration.json': ('ihm.lumbar-muscle-variant.v1', 98),
    'data/models/engineering_supported_rest_v1/registration.json': ('ihm.lumbar-muscle-variant.v1', 98),
}
VARIANT = 'data/models/engineering_stance_v1/registration.json'
# A real, parseable model on disk that is NOT among this manifest's sources.
FOREIGN_MODEL = 'data/derived/mechanics/whole_body_arm26_v2/subject_with_arms.osim'


def manifest(relative):
    return json.loads((ROOT / relative).read_text())


class RegisteredEffectorVariants(unittest.TestCase):
    def test_every_registered_body_loads_with_its_declared_count(self):
        from ihm.assembly.sensorimotor_catalog import whole_body_effector_catalog, REGISTERED_EFFECTOR_SCHEMAS
        for relative, (schema, count) in REGISTERED.items():
            with self.subTest(relative):
                record = manifest(relative)
                self.assertEqual(record['schema'], schema)
                self.assertIn(schema, REGISTERED_EFFECTOR_SCHEMAS)
                rows = whole_body_effector_catalog(ROOT, record)
                self.assertEqual(len(rows), count)
                self.assertEqual(len(rows), record['muscle_count'])
                self.assertEqual(len({r['id'] for r in rows}), count)

    def test_trunk_effectors_appear_only_in_the_variant(self):
        from ihm.assembly.sensorimotor_catalog import whole_body_effector_catalog
        base = {r['id'] for r in whole_body_effector_catalog(
            ROOT, manifest('data/derived/mechanics/whole_body_arm26_v2/registration.json'))}
        variant = {r['id'] for r in whole_body_effector_catalog(ROOT, manifest(VARIANT))}
        self.assertEqual(sorted(variant - base),
                         ['gait2392_ercspn_l', 'gait2392_ercspn_r', 'gait2392_extobl_l',
                          'gait2392_extobl_r', 'gait2392_intobl_l', 'gait2392_intobl_r'])
        self.assertEqual(base - variant, set())
        self.assertEqual(len([r for r in base if r.startswith('arm26_')]), 12)

    def test_peripheral_coverage_reaches_the_variant_body(self):
        from ihm.assembly.peripheral_coverage import build_peripheral_coverage
        for relative, (_, count) in REGISTERED.items():
            with self.subTest(relative):
                report = build_peripheral_coverage(ROOT, registration=manifest(relative))
                self.assertEqual(len(report['effectors']), count)
                json.dumps(report, allow_nan=False)
        report = build_peripheral_coverage(ROOT, registration=manifest(VARIANT))
        trunk = [e for e in report['effectors'] if e['body_group'] == 'lumbar_trunk']
        self.assertEqual(len(trunk), 6)

    def test_variant_fields_are_verified_not_merely_read(self):
        from ihm.assembly.sensorimotor_catalog import whole_body_effector_catalog
        record = manifest(VARIANT)
        mutations = [
            ('schema', 'ihm.invented.v1'),
            ('insert_sha256', '0' * 64),
            ('insert_path', '../outside.xml'),
            ('insert_path', '/etc/hostname'),
            # Exists, parses, is a genuine registered model -- and is not one of
            # THIS manifest's verified sources. An existence check passes it.
            ('base_model_path', FOREIGN_MODEL),
            ('base_model_path', None),
            ('model_sha256', '0' * 64),
            ('catalog_sha256', '0' * 64),
            ('muscle_count', 97),
        ]
        for field, value in mutations:
            with self.subTest(field=field, value=value):
                bad = deepcopy(record)
                bad[field] = value
                with self.assertRaises(ValueError):
                    whole_body_effector_catalog(ROOT, bad)

    def test_a_dropped_source_is_still_refused(self):
        from ihm.assembly.sensorimotor_catalog import whole_body_effector_catalog
        record = manifest(VARIANT)
        row_sources = {r['source_path'] for r in json.loads((ROOT / record['catalog_path']).read_text())}
        self.assertTrue(row_sources)
        for source in sorted(row_sources):
            with self.subTest(source):
                bad = deepcopy(record)
                bad['sources'] = {k: v for k, v in bad['sources'].items() if k != source}
                with self.assertRaises(ValueError):
                    whole_body_effector_catalog(ROOT, bad)

    def test_catalog_is_a_function_of_its_manifest(self):
        """Call it twice at the same input (IHM-1 CLAUDE.md)."""
        from ihm.assembly.sensorimotor_catalog import whole_body_effector_catalog
        record = manifest(VARIANT)
        self.assertEqual(whole_body_effector_catalog(ROOT, record),
                         whole_body_effector_catalog(ROOT, deepcopy(record)))


if __name__ == '__main__':
    unittest.main()
