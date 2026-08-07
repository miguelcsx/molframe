#!/usr/bin/env bash
# Specification checks for docs/ — see docs/20-testing.md §8.
# Run from the repository root. Exits non-zero if any check fails.

set -uo pipefail
cd "$(dirname "$0")/.."

FAIL=0
pass() { printf '  \033[32mPASS\033[0m  %s\n' "$1"; }
fail() { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAIL=1; }
note() { printf '        %s\n' "$1"; }

DOCS=(README.md CONTRIBUTING.md docs/*.md docs/adr/*.md)

echo "== 1. completeness =="
N=$(ls README.md CONTRIBUTING.md docs/*.md docs/adr/*.md | wc -l)
[ "$N" -eq 46 ] && pass "46 documents present" || fail "expected 46 documents, found $N"

# TODO markers, excluding the ADR template and lines that merely name the marker
if grep -rn 'TODO\|TBD\|FIXME' "${DOCS[@]}" 2>/dev/null \
   | grep -v '0000-template' | grep -v '`TODO`' | grep -q .; then
  fail "TODO/TBD/FIXME in a normative section"
  grep -rn 'TODO\|TBD\|FIXME' "${DOCS[@]}" | grep -v '0000-template' | grep -v '`TODO`' | sed 's/^/        /'
else
  pass "no TODO/TBD/FIXME in normative sections"
fi

echo "== 2. requirement traceability =="
grep -o '\*\*\(FR\|NFR\)-[0-9]\+\*\*' docs/01-requirements.md | tr -d '*' | sort -u > /tmp/pdbiox_defined
grep -rho '\(FR\|NFR\)-[0-9]\{3,4\}' "${DOCS[@]}" | sort -u > /tmp/pdbiox_cited
grep -rho '\(FR\|NFR\)-[0-9]\{3,4\}' README.md CONTRIBUTING.md docs/*.md docs/adr/*.md \
  --exclude=01-requirements.md | sort -u > /tmp/pdbiox_elsewhere

UNDEF=$(comm -13 /tmp/pdbiox_defined /tmp/pdbiox_cited)
[ -z "$UNDEF" ] && pass "every cited requirement is defined" \
                || { fail "cited but undefined:"; note "$UNDEF"; }

UNCITED=$(comm -23 /tmp/pdbiox_defined /tmp/pdbiox_elsewhere)
[ -z "$UNCITED" ] && pass "every requirement is cited by a subsystem document" \
                  || { fail "defined but never cited:"; note "$UNCITED"; }

echo "== 3. error codes =="
grep -rho 'PDBIOX-[EW][0-9]\{4\}' "${DOCS[@]}" | sort -u > /tmp/pdbiox_ecited
grep -o '`[EW][0-9]\{4\}`' docs/18-errors-and-diagnostics.md | tr -d '`' \
  | sed 's/^/PDBIOX-/' | sort -u > /tmp/pdbiox_ereg
UNREG=$(comm -23 /tmp/pdbiox_ecited /tmp/pdbiox_ereg)
[ -z "$UNREG" ] && pass "every cited error code is registered" \
               || { fail "unregistered codes:"; note "$UNREG"; }

echo "== 4. ADR references =="
grep -rho 'ADR-[0-9]\{4\}' "${DOCS[@]}" | sort -u > /tmp/pdbiox_adrcited
ls docs/adr/ | grep -o '^[0-9]\{4\}' | sed 's/^/ADR-/' | sort -u > /tmp/pdbiox_adrhave
BAD=$(comm -23 /tmp/pdbiox_adrcited /tmp/pdbiox_adrhave)
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
rm -f /tmp/pdbiox_missing_paths /tmp/pdbiox_skipped_roots
grep -rhoE '(inspo|bio)/[A-Za-z0-9_./-]+' "${DOCS[@]}" | sed 's/[.,;:)`]*$//' | sort -u \
| while read -r p; do
    case "$p" in
      inspo/*)
        if [ ! -d inspo ]; then
          echo "inspo/" >> /tmp/pdbiox_skipped_roots
          continue
        fi
        full="$p"
        ;;
      bio/*)
        relative="${p#bio/}"
        project="${relative%%/*}"
        if [ ! -d "../$project" ]; then
          echo "bio/$project/" >> /tmp/pdbiox_skipped_roots
          continue
        fi
        full="../$relative"
        ;;
    esac
    [ -e "$full" ] || echo "$p"
  done > /tmp/pdbiox_missing_paths
if [ -s /tmp/pdbiox_missing_paths ]; then
  fail "cited reference paths do not exist:"; note "$(cat /tmp/pdbiox_missing_paths)"
else
  pass "every cited path in an installed external reference root exists"
fi
if [ -s /tmp/pdbiox_skipped_roots ]; then
  note "external reference roots not installed: $(sort -u /tmp/pdbiox_skipped_roots | tr '\n' ' ')"
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
