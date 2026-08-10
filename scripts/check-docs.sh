#!/usr/bin/env bash
# Specification checks for docs/ — see docs/20-testing.md §8.
# Run from the repository root. Exits non-zero if any check fails.

set -uo pipefail
cd "$(dirname "$0")/.."

CHECK_DOCS_TMP=$(mktemp -d "${TMPDIR:-/tmp}/pdbiox-docs.XXXXXX") || exit 1
cleanup() { rm -rf "$CHECK_DOCS_TMP"; }
trap cleanup EXIT
trap 'exit 130' HUP INT TERM

FAIL=0
pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAIL=1; }
note() { printf '        %s\n' "$1"; }

DOCS=(README.md CONTRIBUTING.md docs/*.md docs/adr/*.md)
DOC_MANIFEST=docs/MANIFEST.txt

echo "== 1. completeness =="
grep -vE '^[[:space:]]*(#|$)' "$DOC_MANIFEST" | sort -u > "$CHECK_DOCS_TMP/expected-docs"
printf '%s\n' "${DOCS[@]}" | sort -u > "$CHECK_DOCS_TMP/discovered-docs"
if diff -u "$CHECK_DOCS_TMP/expected-docs" "$CHECK_DOCS_TMP/discovered-docs" \
    > "$CHECK_DOCS_TMP/document-diff"; then
  pass "document set matches $DOC_MANIFEST"
else
  fail "document set differs from $DOC_MANIFEST"
  sed 's/^/        /' "$CHECK_DOCS_TMP/document-diff"
fi

# TODO markers, excluding the ADR template and lines that merely name the marker
if grep -rn 'TODO\|TBD\|FIXME' "${DOCS[@]}" 2>/dev/null \
   | grep -v '0000-template' | grep -v '`TODO`' | grep -q .; then
  fail "TODO/TBD/FIXME in a normative section"
  grep -rn 'TODO\|TBD\|FIXME' "${DOCS[@]}" | grep -v '0000-template' | grep -v '`TODO`' | sed 's/^/        /'
else
  pass "no TODO/TBD/FIXME in normative sections"
fi

echo "== 2. requirement traceability =="
grep -o '\*\*\(FR\|NFR\)-[0-9]\+\*\*' docs/01-requirements.md | tr -d '*' | sort -u > "$CHECK_DOCS_TMP/defined"
grep -rho '\(FR\|NFR\)-[0-9]\{3,4\}' "${DOCS[@]}" | sort -u > "$CHECK_DOCS_TMP/cited"
grep -rho '\(FR\|NFR\)-[0-9]\{3,4\}' README.md CONTRIBUTING.md docs/*.md docs/adr/*.md \
  --exclude=01-requirements.md | sort -u > "$CHECK_DOCS_TMP/elsewhere"

UNDEF=$(comm -13 "$CHECK_DOCS_TMP/defined" "$CHECK_DOCS_TMP/cited")
[ -z "$UNDEF" ] && pass "every cited requirement is defined" \
                || { fail "cited but undefined:"; note "$UNDEF"; }

UNCITED=$(comm -23 "$CHECK_DOCS_TMP/defined" "$CHECK_DOCS_TMP/elsewhere")
[ -z "$UNCITED" ] && pass "every requirement is cited by a subsystem document" \
                  || { fail "defined but never cited:"; note "$UNCITED"; }

echo "== 3. error codes =="
grep -rho 'PDBIOX-[EW][0-9]\{4\}' "${DOCS[@]}" | sort -u > "$CHECK_DOCS_TMP/error-cited"
grep -o '`[EW][0-9]\{4\}`' docs/18-errors-and-diagnostics.md | tr -d '`' \
  | sed 's/^/PDBIOX-/' | sort -u > "$CHECK_DOCS_TMP/error-registered"
UNREG=$(comm -23 "$CHECK_DOCS_TMP/error-cited" "$CHECK_DOCS_TMP/error-registered")
[ -z "$UNREG" ] && pass "every cited error code is registered" \
               || { fail "unregistered codes:"; note "$UNREG"; }

echo "== 4. ADR references =="
grep -rho 'ADR-[0-9]\{4\}' "${DOCS[@]}" | sort -u > "$CHECK_DOCS_TMP/adr-cited"
ls docs/adr/ | grep -o '^[0-9]\{4\}' | sed 's/^/ADR-/' | sort -u > "$CHECK_DOCS_TMP/adr-present"
BAD=$(comm -23 "$CHECK_DOCS_TMP/adr-cited" "$CHECK_DOCS_TMP/adr-present")
[ -z "$BAD" ] && pass "every ADR reference resolves" || { fail "dangling:"; note "$BAD"; }

echo "== 5. decision closure =="
# Hedging is permitted only in adr/ (Rejected alternatives), the GLOSSARY rule that
# quotes the banned forms, and the Loose severity definition.
HEDGE=$(grep -rniE '\b(we could|one option|might want|perhaps|maybe|not sure|to be decided)\b' \
        README.md CONTRIBUTING.md docs/*.md 2>/dev/null \
        | grep -v 'GLOSSARY.md.*does not belong' || true)
[ -z "$HEDGE" ] && pass "no hedging in normative prose" || { fail "hedging found:"; note "$HEDGE"; }

echo "== 6. naming =="
STRAY=$(grep -rnoE '\b(pdbio|pbiox)\b' README.md CONTRIBUTING.md docs/*.md docs/adr/*.md 2>/dev/null \
        | grep -vE '(00-vision-and-scope\.md|README\.md|20-testing\.md)' || true)
[ -z "$STRAY" ] && pass "retired names appear only where they are discussed" \
               || { fail "stray retired name:"; note "$STRAY"; }

echo "== 7. reference paths =="
grep -rhoE '(inspo|bio)/[A-Za-z0-9_./-]+' "${DOCS[@]}" | sed 's/[.,;:)`]*$//' | sort -u \
| while read -r p; do
    case "$p" in
      inspo/*)
        if [ ! -d inspo ]; then
          echo "inspo/" >> "$CHECK_DOCS_TMP/skipped-roots"
          continue
        fi
        full="$p"
        ;;
      bio/*)
        relative="${p#bio/}"
        project="${relative%%/*}"
        if [ ! -d "../$project" ]; then
          echo "bio/$project/" >> "$CHECK_DOCS_TMP/skipped-roots"
          continue
        fi
        full="../$relative"
        ;;
    esac
    [ -e "$full" ] || echo "$p"
  done > "$CHECK_DOCS_TMP/missing-paths"
if [ -s "$CHECK_DOCS_TMP/missing-paths" ]; then
  fail "cited reference paths do not exist:"; note "$(cat "$CHECK_DOCS_TMP/missing-paths")"
else
  pass "every cited path in an installed external reference root exists"
fi
if [ -s "$CHECK_DOCS_TMP/skipped-roots" ]; then
  note "external reference roots not installed: $(sort -u "$CHECK_DOCS_TMP/skipped-roots" | tr '\n' ' ')"
fi

echo "== 8. internal links =="
BROKEN=""
for l in $(grep -rhoE '\]\(([A-Za-z0-9_./-]+\.md)\)' "${DOCS[@]}" | sed -E 's/^\]\(//; s/\)$//' | sort -u); do
  found=0
  for base in . docs docs/adr; do [ -e "$base/$l" ] && found=1 && break; done
  [ $found -eq 0 ] && BROKEN="$BROKEN $l"
done
[ -z "$BROKEN" ] && pass "every internal link resolves" || { fail "broken links:"; note "$BROKEN"; }

echo "== 9. parity accounting =="
python3 - <<'PY'
import re, collections, sys
text = open('docs/PARITY.md').read()
body, _, summary = text.partition('## Summary by phase')

# Count capability rows only — the summary table is not part of the register.
rows = collections.Counter()
for line in body.splitlines():
    if not line.startswith('|'): continue
    cells = [c.strip() for c in line.strip().strip('|').split('|')]
    if len(cells) < 3 or set(cells[0]) <= set('-: '): continue
    for c in reversed(cells):
        if re.fullmatch(r'[1-8]', c): rows[int(c)] += 1; break
        if c == '–': rows['-'] += 1; break

claimed = {}
for line in summary.splitlines():
    m = re.match(r'\|\s*([1-8])\s*\|\s*(\d+)\s*\|', line)
    if m: claimed[int(m.group(1))] = int(m.group(2))

ok = True
for p in range(1, 9):
    if rows[p] != claimed.get(p):
        print(f'  \033[31mFAIL\033[0m  phase {p}: table has {rows[p]} rows, summary claims {claimed.get(p)}')
        ok = False
if ok:
    print(f'  \033[32mPASS\033[0m  summary matches {sum(rows[p] for p in range(1,9))} phased rows')
sys.exit(0 if ok else 1)
PY
[ $? -ne 0 ] && FAIL=1

echo
[ $FAIL -eq 0 ] && echo "All specification checks passed." || echo "Specification checks FAILED."
exit $FAIL
