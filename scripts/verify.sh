#!/usr/bin/env bash
# Verifies the pdbiox working tree.
#
# This is the single definition of "green" for the project. The harness
# dispatcher at bio/scripts/verify.sh runs this only when pdbiox is dirty; if
# you are working from pdbiox directly, run it yourself. Every check has a
# reason — see Rules.md §"Checking" for the rationale behind the ordering.
#
# Adjusting what this script runs means changing what "green" means. Don't.

set -uo pipefail

cd "$(dirname "$0")/.."

fail=0
pass() { printf '\033[32mPASS\033[0m  %s\n' "$1"; }
die()  { printf '\033[31mFAIL\033[0m  %s\n' "$1"; fail=$((fail + 1)); }

# 1. Format. Mechanical, fast, catches drift introduced by hand-editing.
if cargo fmt --all --check; then
    pass "cargo fmt --all --check"
else
    die "cargo fmt --all --check (run 'cargo fmt --all')"
fi

# 2. Lints. The workspace permits exactly the warns/allows declared in
#    Cargo.toml; everything else is a regression. -D warnings is the floor.
if cargo clippy --workspace --all-targets -- -D warnings 2>&1 | tee /tmp/pdbiox_clippy; then
    pass "cargo clippy --workspace --all-targets -- -D warnings"
else
    die "cargo clippy --workspace --all-targets -- -D warnings"
fi

# 3. Tests. The whole claim is that this compiles, runs, and is correct on
#    every fixture we care about. cargo test --workspace is the assertion.
if cargo test --workspace 2>&1 | tee /tmp/pdbiox_test; then
    pass "cargo test --workspace"
else
    die "cargo test --workspace"
fi

# 4. Spec checks. Only meaningful when docs/ is present. When this repo is
#    checked out without docs/ (the default for a standalone clone), skip
#    with a note — failing here would block legitimate work, not catch bugs.
if [ -d docs ] && [ -f docs/01-requirements.md ]; then
    if [ -x scripts/check-docs.sh ]; then
        if scripts/check-docs.sh; then
            pass "scripts/check-docs.sh"
        else
            die "scripts/check-docs.sh"
        fi
    fi
else
    printf '\033[33mnote\033[0m:  docs/ not present — skipping scripts/check-docs.sh\n'
fi

# 5. No unwrap outside tests. The whole point of the rule in Rules.md §3.
if grep -rn "unwrap" crates/ 2>/dev/null | grep -v "_tests.rs" > /tmp/pdbiox_unwrap; then
    if [ -s /tmp/pdbiox_unwrap ]; then
        die "unwrap used outside tests:"
        sed 's/^/        /' /tmp/pdbiox_unwrap
    fi
else
    pass "no unwrap outside tests"
fi

# 6. File size ceiling. Rules.md §2 — files cap at roughly 400–500 lines.
oversize="$(find crates -name "*.rs" -exec wc -l {} + | awk '$1>500 {print}')"
if [ -n "$oversize" ]; then
    die "files over the 500-line ceiling:"
    printf '        %s\n' $oversize
else
    pass "every file under the 500-line ceiling"
fi

printf '\n'
if [ "$fail" -eq 0 ]; then
    echo "pdbiox: green"
    exit 0
else
    echo "pdbiox: $fail check(s) failed"
    exit 1
fi
