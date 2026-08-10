#!/usr/bin/env bash
# Verifies the pdbiox working tree.
#
# This is the single definition of "green" for the project. The harness
# dispatcher at bio/scripts/verify.sh runs this only when pdbiox is dirty; if
# you are working from pdbiox directly, run it yourself. Every check has a
# reason — see RULES.md §"Checking" for the rationale behind the ordering.
#
# Adjusting what this script runs means changing what "green" means. Don't.

set -uo pipefail

cd "$(dirname "$0")/.."

VERIFY_TMP=$(mktemp -d "${TMPDIR:-/tmp}/pdbiox-verify.XXXXXX") || exit 1
cleanup() { rm -rf "$VERIFY_TMP"; }
trap cleanup EXIT
trap 'exit 130' HUP INT TERM

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
if cargo clippy --workspace --all-targets -- -D warnings 2>&1 | tee "$VERIFY_TMP/clippy"; then
    pass "cargo clippy --workspace --all-targets -- -D warnings"
else
    die "cargo clippy --workspace --all-targets -- -D warnings"
fi

# 3. Tests. The whole claim is that this compiles, runs, and is correct on
#    every fixture we care about. cargo test --workspace is the assertion.
if cargo test --workspace 2>&1 | tee "$VERIFY_TMP/test"; then
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

# 5. REUSE compliance. `reuse` is optional Python tooling — when it is not
#    installed the check is skipped with a note rather than blocking the Rust
#    loop. Compliance itself is defined by LICENSES/ + REUSE.toml.
if command -v reuse >/dev/null 2>&1; then
    if scripts/check-reuse.sh; then
        pass "scripts/check-reuse.sh"
    else
        die "scripts/check-reuse.sh"
    fi
else
    printf '\033[33mnote\033[0m:  `reuse` not installed — skipping scripts/check-reuse.sh\n'
fi

# 6. The Python boundary is mechanical and fully typed.
if scripts/check-python.py; then
    pass "scripts/check-python.py"
else
    die "scripts/check-python.py"
fi

# 7. Source policy. Layout, size and explicit-error rules apply to the whole
#    user-facing surface, not only to Rust implementation files.
if scripts/check-code-policy.sh; then
    pass "scripts/check-code-policy.sh"
else
    die "scripts/check-code-policy.sh"
fi

printf '\n'
if [ "$fail" -eq 0 ]; then
    echo "pdbiox: green"
    exit 0
else
    echo "pdbiox: $fail check(s) failed"
    exit 1
fi
