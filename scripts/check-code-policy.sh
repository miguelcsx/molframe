#!/usr/bin/env bash
# Checks repository-wide source layout rules that do not require compilation.

set -uo pipefail

cd "$(dirname "$0")/.."

fail=0
pass() { printf '\033[32mPASS\033[0m  %s\n' "$1"; }
die()  { printf '\033[31mFAIL\033[0m  %s\n' "$1"; fail=$((fail + 1)); }

oversize="$({
    find crates python -type f \( \
        -name '*.rs' -o -name '*.py' -o -name '*.pyi' -o -name '*.toml' \
    \) -exec wc -l {} +
    find . -maxdepth 1 -type f -name '*.toml' -exec wc -l {} +
} | awk '$2 != "total" && $1 > 500 {print}')"
if [ -n "$oversize" ]; then
    die "source or configuration files over the 500-line ceiling:"
    printf '        %s\n' "$oversize"
else
    pass "every source and configuration file is under 500 lines"
fi

implementation_pattern='^\s*(pub(\([^)]*\))?\s+)?(async\s+)?(unsafe\s+)?(fn|struct|enum|union|trait|type|const|static|macro_rules!)\b|^\s*impl\b'
module_implementation="$(find crates -type f \( -name lib.rs -o -name mod.rs \) -print0 \
    | xargs -0 rg -n "$implementation_pattern" 2>/dev/null || true)"
if [ -n "$module_implementation" ]; then
    die "implementation found in lib.rs or mod.rs:"
    printf '%s\n' "$module_implementation" | sed 's/^/        /'
else
    pass "lib.rs and mod.rs contain declarations and re-exports only"
fi

allows="$(rg -n '#\s*\[\s*allow' crates --glob '*.rs' 2>/dev/null || true)"
if [ -n "$allows" ]; then
    die "source-level lint bypasses found:"
    printf '%s\n' "$allows" | sed 's/^/        /'
else
    pass "no source-level lint bypasses"
fi

absence_helpers="$(rg -n '\.(unwrap|expect)(_[A-Za-z0-9_]+)?\(' crates \
    --glob '*.rs' --glob '!**/*_tests.rs' --glob '!**/tests/**' 2>/dev/null || true)"
if [ -n "$absence_helpers" ]; then
    die "panic or hidden-absence helpers used outside tests:"
    printf '%s\n' "$absence_helpers" | sed 's/^/        /'
else
    pass "no panic or hidden-absence helpers outside tests"
fi

exit "$fail"
