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
| `pdbiox-geom` | Distances, angles, torsions, centroid/centre of mass, radius of gyration, inertia tensor, principal axes, asphericity, per-atom fluctuation (RMSF), rigid transforms, quaternion superposition. The cyclic-Jacobi eigensolver serves both the inertia tensor (3×3) and the superposition matrix (4×4). |
| `pdbiox-bcif` | BinaryCIF codecs, lazy document, reader and writer with deterministic encoding selection. |
| `pdbiox-chem` | Chemistry data: 118-element properties, named vdW/covalent/ionic radius sets, CCD providers, bond construction with provenance, atom equivalence, residue classification. |
| `pdbiox-xtal` | Assemblies, operator expressions, lazy instances, all 230 space groups over 530 Hall settings, symmetry, fractional↔Cartesian, crystal neighbour search, NCS. |
| `pdbiox-ic` | Internal coordinates (hedra/dihedra) and deterministic Cartesian rebuilding. |
| `pdbiox-query` | Selection DSL: lexer, AST, typed builder, logical/physical plans, glob, macros, connectivity, geometric predicates. |
| `pdbiox-spatial` | Cell list, k-d tree, brute force, neighbour list, planner, PBC (orthorhombic/triclinic), minimum image. |
| `pdbiox-surface` | Solvent-accessible surface (Shrake–Rupley), buried surface, deterministic sphere sampling. |
| `pdbiox-analysis` | Contacts, contact maps, native contacts (Q), hydrogen bonds, salt bridges, chain interfaces. |
| `pdbiox-validate` | Steric clashes, cis-peptide detection, occupancy/B-factor sanity checks. |
| `pdbiox-compare` | Superposition-free lDDT; TM-score, GDT-TS and GDT-HA over the shared superposition. |
| `pdbiox-seq` | Pairwise alignment (global, local, semi-global) with affine gap costs. |
| `pdbiox-ml` | Apache Arrow interop (C-stream export, extension types). |
| `pdbiox-mmap` | The single audited `unsafe` boundary for OS memory mapping. |
| `pdbiox-py` | PyO3 bindings: zero-copy NumPy views, Arrow C-stream, scoped coordinate mutation. |
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

## 5. State on 2026-08-07

Phases 1–3 are implemented at `β` across the crates above; phase-4/6 native
analysis has begun. Every implemented capability sits at `β` (native + tested)
rather than `✓`, because a `✓` needs a row-specific golden workflow and
differential evidence against an external corpus that is not in this checkout.

Recently landed native, deterministic capabilities (all `β`, all with unit and
property tests): RMSF (`pdbiox-geom`); solvent-accessible surface and buried
surface (`pdbiox-surface`); contacts, contact maps, native contacts/Q, hydrogen
bonds, salt bridges, chain interfaces (`pdbiox-analysis`); clashes, cis-peptide
detection, occupancy/B-factor checks (`pdbiox-validate`); pairwise alignment
(`pdbiox-seq`); lDDT, TM-score, GDT-TS/HA (`pdbiox-compare`).

Outstanding deterministic work, roughly in dependency order: the rest of each of
those crates (Lee–Richards/SES/cavities; secondary structure, π-stacking, water
bridges; bond/angle deviation, planarity, chirality, rotamer/valence/stereo,
completeness; MSA, phylogenetics, k-mers, sequence formats, substitution
matrices; chain/atom mapping, DockQ/CAD/QS, CE alignment); `pdbiox-query`
altloc/entity/assembly and chirality selectors; PEOE charges in `pdbiox-chem`;
and the format breadth (trajectories, density maps, remaining structural
formats). Track B (`pdbiox-audit`, `pdbiox-fx`) and adapters are separate.

`docs/PARITY.md` is the row-by-row source of truth; update a capability's row in
the same change that adds it. `docs/ROADMAP.md` §Phase 1 carries the dated
met/not-met table — keep it truthful.
