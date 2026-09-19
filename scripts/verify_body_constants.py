"""The seam between the two bodies is declared in one place. Prove it, twice.

    .venv/bin/python scripts/verify_body_constants.py
    .venv/bin/python scripts/verify_body_constants.py --self-test
    .venv/bin/python scripts/verify_body_constants.py --json

IHM-1 runs two bodies of different stature and different mass.  That is a
modelling fact and this script does not touch it.  What it enforces is that the
numbers describing them are *declared once*, in ``ihm/body_constants.py``, and
that each declaration still matches the artefact it came from.

Four checks, and the fourth is the one that makes the first three worth reading:

A. **No raw literal under ``ihm/``.**  Every ``.py`` file under ``ihm/`` is
   parsed with ``ast`` and every numeric literal compared against the declared
   constants.  A literal is caught however it is spelled -- ``77.6122029``,
   ``7.76122029e1`` and ``77.61220290`` are the same float and grep sees three
   different strings.
B. **No raw digit string under ``ihm/`` either.**  Prose drifts too: a docstring
   that still says 1.7973 m after the constant moved is a reader misled.  This
   is a text scan, and it catches the rounded forms an ``ast`` scan cannot see.
C. **Every constant still matches its source.**  The source-model mass is
   re-added from the ``.osim``, the mechanical stature is re-measured from its
   markers, and the anatomical pair is re-read from ``profile.json`` and from
   the composition ledger underneath it.  A constant that stops matching its
   evidence is a silent change of subject.
D. **The checks can fail.**  ``--self-test`` writes a literal, and a wrong
   mirror, into a scratch copy and asserts that A, B and the mirror check each
   report them.  A control whose pass looks identical to its absence is not
   evidence until it has been made to fail on purpose (IHM-1 ``CLAUDE.md``).

Plus an idempotence check: the scan is run twice over the same tree and the two
results must be equal.  It is one line and it is not a known-answer test -- a
known answer tests the value, idempotence tests whether the instrument is a
function at all.

**What this does NOT police, and why.**

* ``data/derived/**``, ``data/research/**``, ``out/``, ``logs/`` and dated
  entries in ``docs/``.  Those carry the numbers as a record of what a
  particular run was handed.  Rewriting them would falsify the record.
* ``scripts/**``.  22 live literals in 20 scripts still pass the mechanical mass
  by hand (counted 2026-09-18); each is a behaviour-preserving one-line
  replacement, and they are counted in ``docs/BODY_PARAMETERS.md`` rather than
  changed here.  Four of those files did not carry it that morning, so the
  number is still SPREADING through scripts/ -- extending this guard to cover
  them is the next thing to do.
* ``SUPERSEDED_BIOGEARS_STANDARDMALE_MASS_KG`` (77.1107029).  It appears in
  ``ihm/assembly/profile.py`` inside the argument for why it was rejected, and
  in the ``supersedes`` block profile.json records.  Policing a number that must
  appear in the case against it would only force a worse spelling of it.
* ``ANATOMY_REGISTRATION_SCALE`` (0.963) is not in the global scan: three
  significant figures is not distinctive enough to ban repository-wide without
  false alarms.  It is pinned instead by the mirror check on
  ``ihm/assembly/anatomy_pose.py``.
"""
from pathlib import Path
import argparse, ast, json, re, shutil, sys, tempfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from ihm import body_constants as bc                                # noqa: E402

DECLARED_IN = 'ihm/body_constants.py'

#: Values banned as a raw literal anywhere under ``ihm/``.  Each is distinctive
#: to at least seven significant figures, so an unrelated coefficient cannot
#: collide with one by accident.
POLICED = (
    ('MECHANICAL_TARGET_MASS_KG', bc.MECHANICAL_TARGET_MASS_KG),
    ('MECHANICAL_SOURCE_MASS_KG', bc.MECHANICAL_SOURCE_MASS_KG),
    ('MECHANICAL_STATURE_M', bc.MECHANICAL_STATURE_M),
    ('ANATOMICAL_MASS_KG', bc.ANATOMICAL_MASS_KG),
    ('ANATOMICAL_STATURE_M', bc.ANATOMICAL_STATURE_M),
)

#: Digit strings that identify a policed constant in prose, including the
#: rounded forms the docs use.  ``ast`` cannot see these; a reader can.
DIGIT_FORMS = (
    ('MECHANICAL_TARGET_MASS_KG', r'77\.6122'),
    ('MECHANICAL_SOURCE_MASS_KG', r'85\.2698'),
    ('MECHANICAL_STATURE_M', r'1\.79727|1\.7973(?![0-9])'),
    ('ANATOMICAL_MASS_KG', r'70\.7713'),
    ('ANATOMICAL_STATURE_M', r'1\.71947|1\.7195(?![0-9])'),
)

MIRROR_FILES = {path for path, _, _, _ in bc.MIRRORS}


def _py_files(root):
    return sorted(p for p in (root / 'ihm').rglob('*.py')
                  if '__pycache__' not in p.parts)


def non_python_under_ihm(root):
    """Both scans read ``.py`` only, so the claim 'nothing under ihm/ carries a
    literal' depends on ihm/ being entirely Python. It is today (0 non-``.py``
    files). This returns any file that would break that, rather than letting a
    JSON config quietly sit outside the scan.
    """
    return sorted(p.relative_to(root).as_posix()
                  for p in (root / 'ihm').rglob('*')
                  if p.is_file() and p.suffix != '.py'
                  and '__pycache__' not in p.parts)


def scan_literals(root):
    """A: numeric literals under ihm/ equal to a policed constant."""
    hits = []
    for path in _py_files(root):
        rel = path.relative_to(root).as_posix()
        if rel == DECLARED_IN:
            continue
        tree = ast.parse(path.read_text(), filename=rel)
        allowed = _mirror_allowance(rel, tree)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Constant) or isinstance(node.value, bool):
                continue
            if not isinstance(node.value, (int, float)):
                continue
            for name, value in POLICED:
                if node.value == value and allowed.get(node.lineno) != value:
                    hits.append({'file': rel, 'line': node.lineno,
                                 'constant': name, 'value': repr(node.value)})
    return hits


def _mirror_allowance(rel, tree):
    """``{line: value}`` this file is permitted to carry, from its MIRRORS rows.

    Keyed by the file being scanned, so a mirror declared for one file cannot
    excuse a literal in another, and by value, so the allowance covers exactly
    the constant the mirror declares and nothing else on that line.
    """
    allowed = {}
    for path, target, constant, _ in bc.MIRRORS:
        if path != rel:
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and any(
                    isinstance(t, ast.Name) and t.id == target for t in node.targets):
                allowed[node.value.lineno] = getattr(bc, constant)
    return allowed


def scan_text(root):
    """B: digit strings under ihm/, outside the declaration and its mirrors."""
    hits = []
    for path in _py_files(root):
        rel = path.relative_to(root).as_posix()
        if rel == DECLARED_IN or rel in MIRROR_FILES:
            continue
        for number, line in enumerate(path.read_text().splitlines(), 1):
            for name, pattern in DIGIT_FORMS:
                if re.search(pattern, line):
                    hits.append({'file': rel, 'line': number, 'constant': name,
                                 'text': line.strip()[:120]})
    return hits


def check_mirrors(root):
    """The declared copies must still equal their constant, exactly."""
    rows = []
    for path, target, constant, reason in bc.MIRRORS:
        source = root / path
        found = None
        if source.exists():
            for node in ast.walk(ast.parse(source.read_text(), filename=path)):
                if isinstance(node, ast.Assign) and any(
                        isinstance(t, ast.Name) and t.id == target
                        for t in node.targets):
                    found = ast.literal_eval(node.value)
        expected = getattr(bc, constant)
        rows.append({'file': path, 'name': target, 'constant': constant,
                     'declared': repr(expected), 'mirrored': repr(found),
                     # exact float equality on purpose: a mirror is a copy, not
                     # an approximation, and a tolerance here would hide a typo
                     'ok': found is not None and found == expected,
                     'why_not_an_import': reason})
    return rows


def check_sources(root):
    """C: every constant against the artefact it came from."""
    import xml.etree.ElementTree as ET
    from ihm.native.model_scaling import head_marker_height_m

    rows = []

    def row(name, declared, measured, tol, source):
        rows.append({'constant': name, 'declared': repr(declared),
                     'measured': repr(measured), 'tolerance': tol,
                     'source': source,
                     'ok': measured is not None and abs(declared - measured) <= tol})

    model = ET.parse(root / bc.MECHANICAL_MODEL).getroot()
    row('MECHANICAL_SOURCE_MASS_KG', bc.MECHANICAL_SOURCE_MASS_KG,
        sum(float(b.findtext('mass')) for b in model.iter('Body')), 1e-9,
        bc.MECHANICAL_MODEL + ' (sum of 22 <Body><mass>)')
    row('MECHANICAL_STATURE_M', bc.MECHANICAL_STATURE_M,
        head_marker_height_m(root / bc.MECHANICAL_MODEL), 1e-9,
        bc.MECHANICAL_MODEL + ' (head_marker_height_m)')

    profile_path = root / bc.ANATOMICAL_PROFILE
    if profile_path.exists():
        profile = json.loads(profile_path.read_text())
        row('ANATOMICAL_STATURE_M', bc.ANATOMICAL_STATURE_M,
            float(profile['height_m']), 1e-7, bc.ANATOMICAL_PROFILE)
        row('ANATOMICAL_MASS_KG', bc.ANATOMICAL_MASS_KG,
            float(profile['mass_kg']), 1e-7, bc.ANATOMICAL_PROFILE)
    else:
        row('ANATOMICAL_STATURE_M', bc.ANATOMICAL_STATURE_M, None, 1e-7,
            bc.ANATOMICAL_PROFILE + ' (ABSENT)')
        row('ANATOMICAL_MASS_KG', bc.ANATOMICAL_MASS_KG, None, 1e-7,
            bc.ANATOMICAL_PROFILE + ' (ABSENT)')

    ledger_path = root / 'data/derived/interstitial-composition-prior-v1/ledger.json'
    if ledger_path.exists():
        ledger = json.loads(ledger_path.read_text())
        row('ANATOMICAL_MASS_KG (ledger)', bc.ANATOMICAL_MASS_KG,
            float(ledger['verdict']['composed_total_body_mass_kg']), 5e-5,
            str(ledger_path.relative_to(root)))

    # MECHANICAL_TARGET_MASS_KG: re-derived two independent ways, because it was
    # recorded as having no derivation until 2026-09-18 and a prose claim about
    # provenance is worth nothing that a run cannot reproduce.
    #
    # (i) the arithmetic identity, in EXACT float -- 170 lb plus the stomach
    #     BioGears seeds at t=0. Tolerance 0, on purpose: an identity that needs
    #     a tolerance is a coincidence.
    identity = (bc.BIOGEARS_STANDARDMALE_WEIGHT_LB * bc.POUND_KG
                + bc.BIOGEARS_INITIAL_STOMACH_KG)
    rows.append({'constant': 'MECHANICAL_TARGET_MASS_KG (identity)',
                 'declared': repr(bc.MECHANICAL_TARGET_MASS_KG),
                 'measured': repr(identity), 'tolerance': 0,
                 'source': '%r lb x %r kg/lb + %r kg stomach, exact float'
                           % (bc.BIOGEARS_STANDARDMALE_WEIGHT_LB, bc.POUND_KG,
                              bc.BIOGEARS_INITIAL_STOMACH_KG),
                 'ok': identity == bc.MECHANICAL_TARGET_MASS_KG})
    rows.append({'constant': 'SUPERSEDED_BIOGEARS_STANDARDMALE_MASS_KG',
                 'declared': repr(bc.SUPERSEDED_BIOGEARS_STANDARDMALE_MASS_KG),
                 'measured': repr(bc.BIOGEARS_STANDARDMALE_WEIGHT_LB * bc.POUND_KG),
                 'tolerance': 0,
                 'source': '170 lb in kg, exact float',
                 'ok': (bc.BIOGEARS_STANDARDMALE_WEIGHT_LB * bc.POUND_KG
                        == bc.SUPERSEDED_BIOGEARS_STANDARDMALE_MASS_KG)})

    # (ii) the artefact it is written in, read as XML rather than grepped.
    state = root / bc.BIOGEARS_STANDARDMALE_STATE
    found = None
    if state.exists():
        node = ET.parse(state).getroot().find(
            './/{uri:/mil/tatrc/physiology/datamodel}Weight')
        if node is not None and node.get('unit') == 'kg':
            found = float(node.get('value'))
    row('MECHANICAL_TARGET_MASS_KG (state file)', bc.MECHANICAL_TARGET_MASS_KG,
        found, 0, bc.BIOGEARS_STANDARDMALE_STATE + ' <Weight unit="kg">')
    return rows


def check_reexports():
    """ihm.body_parameters must still hand out the same objects, not copies."""
    from ihm import body_parameters as bp
    rows = []
    for name, value in POLICED:
        rows.append({'name': name, 'ok': getattr(bp, name, None) == value,
                     'value': repr(getattr(bp, name, None))})
    for name in ('ANATOMY_REGISTRATION_SCALE', 'STATURE_DISAGREEMENT',
                 'MASS_DISAGREEMENT', 'MECHANICAL_MODEL', 'ANATOMICAL_PROFILE'):
        rows.append({'name': name,
                     'ok': getattr(bp, name, None) == getattr(bc, name),
                     'value': repr(getattr(bp, name, None))})
    return rows


def self_test():
    """D: make each check fail on purpose, in a scratch copy.

    Without this the report 'no raw literals found' is indistinguishable from
    'the scanner looked at nothing'.
    """
    outcomes = []
    with tempfile.TemporaryDirectory() as tmp:
        fake = Path(tmp) / 'tree'
        (fake / 'ihm' / 'assembly').mkdir(parents=True)
        shutil.copy(ROOT / DECLARED_IN, fake / DECLARED_IN)

        planted = fake / 'ihm' / 'planted.py'
        planted.write_text(
            '# a literal spelled so that grep would miss it\n'
            'TARGET = 7.76122029e1\n')
        literal_hits = scan_literals(fake)
        outcomes.append(('A catches an exponent-spelled literal',
                         any(h['file'] == 'ihm/planted.py' for h in literal_hits)))

        planted.write_text('"""the body is 1.7973 m tall."""\n')
        text_hits = scan_text(fake)
        outcomes.append(('B catches a rounded number in a docstring',
                         any(h['file'] == 'ihm/planted.py' for h in text_hits)))

        planted.unlink()
        outcomes.append(('A and B pass on a clean tree',
                         scan_literals(fake) == [] and scan_text(fake) == []))

        mirror_rel = bc.MIRRORS[0][0]
        target = bc.MIRRORS[0][1]
        (fake / mirror_rel).write_text('%s = 1.8\n' % target)
        rows = check_mirrors(fake)
        outcomes.append(('mirror check catches a drifted copy',
                         any(r['file'] == mirror_rel and not r['ok'] for r in rows)))

        (fake / mirror_rel).write_text('%s = %r\n' % (target, bc.ANATOMICAL_STATURE_M))
        rows = check_mirrors(fake)
        outcomes.append(('mirror check passes on an exact copy',
                         any(r['file'] == mirror_rel and r['ok'] for r in rows)))
    return outcomes


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--json', action='store_true')
    ap.add_argument('--self-test', action='store_true',
                    help='prove each check can fail, then run normally')
    args = ap.parse_args()

    failures = []

    if args.self_test:
        print('SELF-TEST -- each check made to fail on purpose')
        for label, passed in self_test():
            print('  %-46s %s' % (label, 'yes' if passed else 'NO'))
            if not passed:
                failures.append('self-test: ' + label)
        print()

    literals = scan_literals(ROOT)
    # Call it twice at the same input: a scan that is not a function of the tree
    # cannot be trusted about the tree.
    if scan_literals(ROOT) != literals or scan_text(ROOT) != scan_text(ROOT):
        failures.append('the scan is not idempotent; it is not a function of the tree')
    texts = scan_text(ROOT)
    unscanned = non_python_under_ihm(ROOT)
    mirrors = check_mirrors(ROOT)
    sources = check_sources(ROOT)
    reexports = check_reexports()

    report = {'schema': 'ihm.body-constants-guard.v1',
              'declared_in': DECLARED_IN,
              'files_scanned': len(_py_files(ROOT)),
              'raw_literals': literals, 'raw_digit_strings': texts,
              'non_python_under_ihm': unscanned,
              'mirrors': mirrors, 'sources': sources, 'reexports': reexports,
              'two_bodies': {
                  'stature_disagreement': bc.STATURE_DISAGREEMENT,
                  'mass_disagreement': bc.MASS_DISAGREEMENT,
                  'note': bc.KNOWN_SEAM}}

    if args.json:
        print(json.dumps(report, indent=2))
    else:
        print('A. raw literals under ihm/ outside %s: %d (of %d files scanned)'
              % (DECLARED_IN, len(literals), report['files_scanned']))
        for hit in literals:
            print('   FAIL %(file)s:%(line)s carries %(constant)s = %(value)s' % hit)
        print('B. raw digit strings under ihm/: %d' % len(texts))
        for hit in texts:
            print('   FAIL %(file)s:%(line)s %(text)s' % hit)
        print('B2. files under ihm/ the scans cannot read (non-.py): %d'
              % len(unscanned))
        for name in unscanned:
            print('   FAIL %s is outside both scans' % name)
        print('C. declared mirrors (copies that are checked, not removed):')
        for row in mirrors:
            print('   %-4s %s:%s == %s (%s)'
                  % ('ok' if row['ok'] else 'FAIL', row['file'], row['name'],
                     row['constant'], row['mirrored']))
        print('D. constants against their sources:')
        for row in sources:
            print('   %-4s %-28s declared %s  measured %s'
                  % ('ok' if row['ok'] else 'FAIL', row['constant'],
                     row['declared'], row['measured']))
        print('E. ihm.body_parameters re-exports: %d checked, %d wrong'
              % (len(reexports), sum(1 for r in reexports if not r['ok'])))
        for row in reexports:
            if not row['ok']:
                print('   FAIL %(name)s is %(value)s' % row)
        print()
        print('The two bodies remain two bodies: %.2f%% apart in stature, '
              '%.1f%% in mass. Unifying them is a modelling decision and this '
              'guard does not ask for it.'
              % (100 * bc.STATURE_DISAGREEMENT, 100 * bc.MASS_DISAGREEMENT))

    if literals:
        failures.append('%d raw literal(s) under ihm/' % len(literals))
    if texts:
        failures.append('%d raw digit string(s) under ihm/' % len(texts))
    if unscanned:
        failures.append('%d non-.py file(s) under ihm/, outside both scans'
                        % len(unscanned))
    failures += ['mirror %s:%s' % (r['file'], r['name'])
                 for r in mirrors if not r['ok']]
    failures += ['source %s' % r['constant'] for r in sources if not r['ok']]
    failures += ['re-export %s' % r['name'] for r in reexports if not r['ok']]

    if failures:
        print('\nverify_body_constants: FAIL -- ' + '; '.join(failures),
              file=sys.stderr)
        return 1
    # --json must stay machine-parseable: the verdict goes to stderr there.
    print('verify_body_constants: pass', file=sys.stderr if args.json else sys.stdout)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
