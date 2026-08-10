# Rules

House rules for code in this repository. They are not style preferences with
exceptions — they are constraints, and the build is configured to keep them.

`CONTRIBUTING.md` covers *what* to build and what a change must be traceable to.
This covers *how* the code is written.

---

## 1. Comments and docstrings never cite the specification

No requirement numbers, no document names, no section references, no internal
record numbers — nothing that points at a document living outside the code. Not
`REQ-306`, not `storage-layout.md §7`, not `DEC-0009`.

Explain the reasoning in its own terms instead. Docstrings should still be
thorough — the rule removes a pointer, not the explanation it pointed at.

**Why.** The specification is a separate artefact with its own lifecycle. It
will be renumbered, split and rewritten, and every citation in the code becomes
a lie the moment it is. A comment that says *why* stays true; a comment that
says *where it was decided* does not.

**In practice.**

```rust
// Wrong — cites the spec.
/// Both namespaces are stored (REQ-101, see identifiers.md §3).

// Right — says the thing.
/// Both identifier namespaces are stored. Choosing between them at read time is
/// what makes a result impossible to trace back to either the file or the paper,
/// so nothing here chooses.
```

Write the docstring as if the specification did not exist and the reader has only
the code in front of them.

---

## 2. Files cap at roughly 400–500 lines

When a module approaches the cap, split it into a directory module before adding
more. Do not let one file grow past it and plan to tidy later.

**Why.** A file that outgrows a screenful of structure has usually outgrown its
single responsibility too. The cap forces the split at the point where the seam
is still obvious.

**In practice.** `foo.rs` becomes `foo/` with `mod.rs` plus one file per
responsibility. See `crates/pdbiox-core/src/topology/` and
`crates/pdbiox-cif/src/lower/` for the shape.

---

## 3. No `unwrap`, `expect`, or any `unwrap_*` variant outside tests

This includes `unwrap_or`, `unwrap_or_default` and `unwrap_or_else`. Handle
absence with an explicit `match` or `let ... else`.

**Why.** Writing the absent case out forces you to decide what it means *here*,
and puts that decision where a reader will find it. `unwrap_or_default()` hides
the decision behind a type's default, which is rarely the answer the domain
wants — an absent occupancy is not zero occupancy.

**In practice.**

```rust
// Wrong.
let name = atom.name().unwrap_or_default();

// Right — and the empty string is now a visible choice.
let name = match atom.name() {
    Some(name) => name,
    None => "",
};
```

Clippy pushes the other way, so two of its lints are set to `allow` in the
workspace manifest specifically to keep this rule enforceable:

```toml
manual_unwrap_or_default = "allow"
manual_unwrap_or         = "allow"
```

**Do not "fix" those allows.** They exist because of this rule.

Tests are exempt: a test that cannot unwrap is a test whose failure message got
worse. Prefer `let ... else { panic!("...") }` there anyway, so the failure says
what was expected.

---

## 4. Avoid clones and allocations; performance is first-class

Reach for a fixed-size buffer, an index or a borrowed slice before a `String` or
a `Vec`. Nothing allocates per atom, per row or per iteration.

**Why.** This library's whole claim is that a representation choice makes breadth
affordable. An allocation in a per-atom path forfeits the claim regardless of how
correct the result is.

**In practice.**

- Interned identifiers, not strings, in anything a kernel touches.
- One arena for many small strings, not one allocation each.
- A reserved sentinel instead of `Option` where the column is per-atom or
  per-residue — see `crates/pdbiox-core/src/optional.rs`.
- Compare squared distances against a squared cutoff; take the square root only
  when a caller wants a length.
- Read the data-design reasoning in the specification before inventing a new
  layout. The chunked columnar store, the dictionary encoding and late
  materialisation are there for reasons that are written down.

---

## 5. Minimal asymptotic cost, and say what it is

Choose the algorithm whose growth matches the workload, and state the cost in the
module docstring where it is not obvious.

**Why.** A quadratic loop over atoms is invisible on a test fixture and fatal on
a ribosome. Stating the cost makes a regression something a reviewer can catch by
reading rather than by profiling.

**In practice.** Prefer a binary search over a scan where the data is ordered;
prefer one pass with an accumulator over two passes; bound any loop that could
otherwise depend on malformed input. Where a structure carries a summary that
lets work be skipped entirely — per-chunk statistics, bounding boxes, element
sets — use it before reading rows.

---

## 6. `lib.rs` and `mod.rs` contain no implementation

Only `mod` declarations, `pub use` re-exports, and the module docstring. A
`#[cfg(test)] mod fixture;` declaration is fine.

**Why.** It makes the shape of a crate readable in one file, and it means moving
an item between modules never touches the file that publishes it.

---

## 7. Tests live in a sibling file

Declared from the source file:

```rust
#[cfg(test)]
#[path = "thing_tests.rs"]
mod tests;
```

so `thing.rs` is paired with `thing_tests.rs`.

**Why.** It keeps the implementation file at its real length, and it keeps the
tests' access to private items — which a `tests/` directory would not.

Test names are sentences that state the expected behaviour, so a failure is
readable without opening the file:

```rust
#[test]
fn a_leading_space_distinguishes_c_alpha_from_calcium() { ... }
```

---

## 8. Deterministic, complete, and DRY

**Deterministic.** The same input produces the same output, including the order
of any findings, with no dependence on hash iteration order, thread count or
allocator. Anything that fingerprints or sorts must do so on a stable key.

**Complete.** Prefer capabilities that can actually be finished and verified over
ones that need data or an environment that is not here. A crate that compiles,
passes its tests and does what its docstrings claim is worth more than three that
are half-written.

**DRY.** One implementation per idea. If two places need the same thing, extract
it — but extract the *idea*, not a superficial similarity. The symmetric
eigensolver in `pdbiox-geom` serves both the inertia tensor and the superposition
matrix because they genuinely want the same computation; two format writers that
happen to both emit text do not.

---

## Checking

Everything above is checked by:

```bash
cargo clippy --workspace --all-targets    # must be zero warnings
cargo fmt --all --check
cargo test --workspace
./scripts/check-docs.sh

grep -rn "unwrap" crates/ | grep -v "_tests.rs"          # must be empty
find crates -name "*.rs" | xargs wc -l | awk '$1>500'    # must be empty
```
