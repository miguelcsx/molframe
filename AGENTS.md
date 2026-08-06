# AGENTS.md — pdbiox

A batteries-included structural bioinformatics engine in Rust. This file extends
`bio/AGENTS.md` (the harness contract) — it does not restate it. Anything not
overridden here falls back to the outer file.

---

## 1. The crate map

| Crate | Responsibility |
|---|---|
| `pdbiox-core` | Chunked columnar storage, topology tables, `Structure`/`StructureView`, hierarchy handles, `AtomSelection`, diagnostics, `AnalysisPolicy`/`Analysis<T>`/`Provenance`, format dispatch. `#![forbid(unsafe_code)]`. |
| `pdbiox-cif` | mmCIF lexer, `Document` (lossless, order-preserving), lowering with residue-boundary rules, canonical and document-preserving writers. |
| `pdbiox-pdb` | PDB fixed-column reader (hybrid-36, insertion codes, multi-model) and a writer that refuses rather than truncates. |
| `pdbiox-geom` | Distances, angles, torsions, centroid/centre of mass, radius of gyration, inertia tensor, principal axes, asphericity, rigid transforms, quaternion superposition. The cyclic-Jacobi eigensolver serves both the inertia tensor (3×3) and the superposition matrix (4×4). |
| `pdbiox` | Facade crate re-exporting the public surface. |
| `pdbiox-cli` | `info`, `convert`, `validate`, `measure`, `rmsd`, `policy`. |

`Rules.md` is the binding style guide for code in these crates. Read it before
editing.

---

## 2. What lives here, what does not

| Path | In git? | Why |
|---|---|---|
| `crates/`, `Cargo.toml`, `Cargo.lock`, `README.md`, `CONTRIBUTING.md`, `Rules.md`, `rust-toolchain.toml`, `scripts/` | yes | The published surface. |
| `docs/`, `idea/`, `inspo/` | no — `info/exclude` | Spec, design notes, and reference libraries. Re-publish deliberately, not by accident. |
| `target/` | no — `info/exclude` | Build output. |

A standalone clone of this repo has its own `AGENTS.md`, `.githooks/`, and
`scripts/verify.sh` — vendored from `bio/` by `scripts/sync-hooks.sh`. No
outside dependency is needed for the loop.

---

## 3. Scopes

The conventional-commit scope is a crate name or one of the cross-cutting
sections below. A scope present in a commit subject must be in
`.githooks/scopes`.

```
core cif pdb geom cli py query spatial chem xtal surface traj compare validate
spec repo ci
```

- `py` — Python bindings (PyO3 + NumPy zero-copy views + `python/` stubs).
- `query` — selection DSL.
- `spatial` — cell list, k-d tree.
- `chem` — chemistry data (atoms, bonds, residues, polymer types).
- `xtal` — crystallography (cell, symmetry, assembly).
- `surface` — solvent-accessible / molecular surfaces.
- `traj` — trajectory / multi-frame handling.
- `compare` — structural alignment, RMSD, similarity.
- `validate` — diagnostic registry, schema validation.
- `spec` — specification documents at `docs/`.
- `repo` — workspace manifests, CI, developer-facing setup.
- `ci` — CI configuration only.

---

## 4. Verification

A working tree is green when:

```bash
cd pdbiox
./scripts/verify.sh
```

returns zero. That runs (in order):

1. `cargo fmt --all --check`
2. `cargo clippy --workspace --all-targets -- -D warnings`
3. `cargo test --workspace`
4. `scripts/check-docs.sh` — the spec checks. **Requires `docs/`.** When this
   repo is checked out standalone with `docs/` excluded, this step is skipped
   with a note rather than failing.
5. `grep -rn "unwrap" crates/ | grep -v "_tests.rs"` — must be empty.
6. `find crates -name "*.rs" | xargs wc -l | awk '$1>500'` — must be empty.

---

## 5. State on 2026-08-06

**Phase 1 of 2 done** — the 33 of 44 phase-1 capability rows are implemented
across the six crates above. Green: 361 tests, zero clippy-pedantic warnings,
`cargo fmt --check` clean. Outstanding, in unblocking order: `pdbiox-py`,
edit overlay and altloc resolution, mmap/streaming input, golden workflows and
property tests, the differential harness, `pdbiox-query`, `pdbiox-spatial`,
`pdbiox-bcif`.

`docs/ROADMAP.md` §Phase 1 carries the dated met/not-met table. Keep it
truthful — the table is the single source of "are we done" for the phase.
