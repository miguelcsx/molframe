# Contributing to pdbiox

pdbiox is **specification-first**: the specification is written before the code, and it is normative. If you are about to write code that the specification does not describe, the specification change comes first.

Before your first contribution, read the vision-and-scope (what pdbiox is and, importantly, is not), the architecture (the crate graph and the inward-dependency rule) and the requirements (the contract every change is traceable to).

---

## The eight gates

Every pull request passes all eight. They are listed in the order they are usually failed.

### 1. Traceable to a requirement

Your change implements an `FR-nnn` or `NFR-nnn` from the requirements, or it adds one. Cite it in the PR description.

A change that is not traceable is either scope creep or a missing requirement, and both need discussion before code.

### 2. `PARITY.md` updated in the same change

**NFR-608.** Adding or removing a capability updates its row in `PARITY.md` — status, owning crate, phase. Not in a follow-up.

A `✓` requires a golden workflow (gate 3). A `≈` requires a phase for going native. A `–` requires a written argument.

### 3. A golden workflow

New capabilities get an entry in the benchmarks catalogue — executable, measured, doubling as a correctness fixture. This is what makes a `✓` in `PARITY.md` mean something.

### 4. A differential test

If a reference library does this, compare against it under matched policies (testing conventions §5). Classify every divergence:

- **pdbiox is wrong** — fix it, add a regression test
- **the reference is wrong** — document it, report upstream, assert pdbiox's behaviour
- **a legitimate policy difference** — document in `PARITY.md`, verify pdbiox reproduces the reference under the matching policy
- **numerical tolerance** — document the tolerance and why it is acceptable

**An unexplained divergence blocks the merge.** This is the gate that turns "we're different" into either a fix or a finding.

### 5. Documentation with a runnable example

**NFR-503, NFR-504.** Every public item is documented and carries an example that compiles and runs. `missing_docs` is a CI error.

### 6. No logic in the binding layer

**ADR-0013.** `pdbiox-py` bodies are argument conversion, one call into a `pdbiox-*` crate, result conversion. Nothing else. `python/pdbiox/` contains only stubs and re-exports.

If a Python convenience is worth having, put it in Rust, where the CLI and Rust users get it too.

CI enforces this with a dependency allowlist and a Python-source check.

### 7. Licence check on adapted code

**NFR-607.** Adapting third-party code? Name the source, its licence and the attribution in the PR.

**MDAnalysis is GPL-2.0-or-later. Its code must not enter pdbiox.** Read it for design, cite it, write your own. The same applies to MDTraj (LGPL) and OpenStructure (LGPL); gemmi is MPL-2.0 and file-level copyleft, so reimplement rather than adapt.

pdbtbx (MIT), Biopython and Biotite (BSD) are adaptable with attribution.

Full table in the governance document §1.1.

### 8. An ADR for architectural decisions

**NFR-606.** Changing how something fundamental works? Write an ADR first, using the ADR template.

ADRs are immutable once accepted — a change of mind is a *new* ADR that supersedes the old one. And they must include **rejected alternatives with reasons**; an ADR listing only benefits is incomplete and will be sent back.

---

## Development

```bash
git clone <repo> && cd pdbiox

cargo build --workspace
cargo test  --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check

# Python
maturin develop
pytest python/tests
mypy --strict python/

# benchmarks (dedicated machine; shared runners give noise, not data)
cargo bench

# specification checks (requires the spec tree)
./scripts/check-docs.sh
```

## Code conventions

Taken from `bio/atpts/`, which is the house style.

**Modules:** one responsibility each, with a `//!` header stating the scientific intent. Not "this module handles PDB files" — *what it computes and why*.

**Errors:** the crate error type with a stable code, a location, a cause and a remedy. Never a bare string. Never a panic on data.

**Tests:** `#[cfg(test)] mod tests` at the bottom of every module, with sentence-style names that state the expected behaviour:

```rust
#[test]
fn residue_boundaries_fall_back_to_file_order_when_label_seq_id_repeats() { ... }
```

A test name should make a CI failure readable without opening the file.

**Comments:** explain *why*, not *what*. The most valuable comments in this codebase will be the ones explaining why a reference implementation behaves the way it does, and why pdbiox agrees or diverges.

**Allocation:** none in hot loops. `SmallVec` for small collections. Reuse buffers across iterations.

**Kernels:** no string comparison (dictionary ids), no `sqrt` (compare squared distances), no unconditional parallelism (threshold-gate it).

Full anti-pattern list: see the performance notes before writing kernels.

## Where things go

| Change | Crate |
|---|---|
| Types, storage, policy, provenance | `pdbiox-core` |
| Selection language | `pdbiox-query` |
| A format | `pdbiox-cif` / `-bcif` / `-pdb` / `-modelcif` / `-traj` |
| Chemistry, CCD, bonds | `pdbiox-chem` |
| A geometric kernel | `pdbiox-geom` |
| Neighbour search | `pdbiox-spatial` |
| SASA, SES, buried surface | `pdbiox-surface` |
| An analysis | `pdbiox-analysis` |
| A comparison metric | `pdbiox-compare` |
| A validation check | `pdbiox-validate` |
| Assemblies, symmetry, maps | `pdbiox-xtal` |
| Alignment, phylogeny | `pdbiox-seq` |
| Tensor or graph export | `pdbiox-ml` |
| A library bridge | `pdbiox-adapters` |

Unsure? The layering rules usually answer it: put it in the lowest layer that can hold it without adding a dependency.

## Reporting bugs

Include: pdbiox version, platform, a minimal input file (or its PDB ID), the exact code, what you expected, what happened, and the full diagnostic including the `PDBIOX-` code.

**A wrong number is more serious than a crash.** If pdbiox produced a plausible but incorrect result, say so prominently — that is the failure mode this project exists to prevent, and it goes to the front of the queue.

## Proposing capabilities

1. Check `PARITY.md` — it may already have a row and a phase.
2. Check the scope — it may be deliberately another domain.
3. Open a discussion with: the requirement it satisfies, the crate it belongs in, the prior art, and the golden workflow that would demonstrate it.

Under ADR-0014 the default answer to "should pdbiox do this?" is **yes, if it is structural bioinformatics**. The scope boundary is domain, not effort.

## Reviewing

Reviewers check the eight gates, and then:

- Does it hold the complexity budget?
- Is it deterministic — independent of thread count, hash order and allocator?
- Does it compute `Coverage` honestly?
- Does it record what it decided on the user's behalf?
- Does it fail loudly rather than silently losing data?

That last one is the project's character. A library that truncates `AA` to `A` and carries on is worse than one that refuses.
