# AGENTS.md — molframe

A batteries-included structural bioinformatics engine in Rust. This file is
self-contained: it is the whole contract for this repository.

---

## 1. The crate map

| Crate | Responsibility |
|---|---|
| `molframe-core` | Chunked columnar storage, topology tables, `Structure`/`StructureView`, hierarchy handles, `AtomSelection`, diagnostics, `AnalysisPolicy`/`Analysis<T>`/`Provenance`, format dispatch. `#![forbid(unsafe_code)]`. |
| `molframe-cif` | mmCIF lexer, `Document` (lossless, order-preserving), lowering with residue-boundary rules, canonical and document-preserving writers. |
| `molframe-pdb` | PDB fixed-column reader (hybrid-36, insertion codes, multi-model) and a writer that refuses rather than truncates. |
| `molframe-geom` | Distances, angles, torsions, centroid/centre of mass, radius of gyration, inertia tensor, principal axes, asphericity, per-atom fluctuation (RMSF), rigid transforms, quaternion superposition. The cyclic-Jacobi eigensolver serves both the inertia tensor (3×3) and the superposition matrix (4×4). |
| `molframe-bcif` | BinaryCIF codecs, lazy document, reader and writer with deterministic encoding selection. |
| `molframe-chem` | Chemistry data: 118-element properties, named vdW/covalent/ionic radius sets, CCD providers, bond construction with provenance, atom equivalence, residue classification. |
| `molframe-xtal` | Assemblies, operator expressions, lazy instances, all 230 space groups over 530 Hall settings, symmetry, fractional↔Cartesian, crystal neighbour search, NCS. |
| `molframe-ic` | Internal coordinates (hedra/dihedra) and deterministic Cartesian rebuilding. |
| `molframe-query` | Selection DSL: lexer, AST, typed builder, logical/physical plans, glob, macros, connectivity, geometric predicates. |
| `molframe-spatial` | Cell list, k-d tree, brute force, neighbour list, planner, PBC (orthorhombic/triclinic), minimum image. |
| `molframe-surface` | Solvent-accessible surface (Shrake–Rupley), buried surface, deterministic sphere sampling. |
| `molframe-analysis` | Contacts, contact maps, native contacts (Q), hydrogen bonds, salt bridges, chain interfaces. |
| `molframe-validate` | Steric clashes, cis-peptide detection, occupancy/B-factor sanity checks. |
| `molframe-compare` | Superposition-free lDDT; TM-score, GDT-TS and GDT-HA over the shared superposition. |
| `molframe-seq` | Pairwise alignment (global, local, semi-global) with affine gap costs. |
| `molframe-interop` | Apache Arrow interop (C-stream export, extension types). |
| `molframe-modelcif` | ModelCIF metadata and confidence metrics. |
| `molframe-audit` | Bounded policy-sensitivity audits for analyses. |
| `molframe-fx` | Declarative functional-geometry evaluation. |
| `molframe-adapters` | Neutral topology transfer and verified byte acquisition. |
| `molframe-bench` | Shared real-structure fixtures for benchmarks. Not published. |
| `molframe-resource-bench` | Allocation and resident-memory measurement. Not published. |
| `molframe-mmap` | The single audited `unsafe` boundary for OS memory mapping. |
| `molframe-py` | PyO3 bindings: zero-copy NumPy views, Arrow C-stream, scoped coordinate mutation. |
| `molframe` | Facade crate re-exporting the public surface. |
| `molframe-cli` | `info`, `convert`, `validate`, `measure`, `rmsd`, `policy`. |

`RULES.md` is the binding style guide for code in these crates. Read it before
editing.

**A known, accepted upstream duplication.** `molframe-traj` depends on
`hoomd-gsd` (from `hoomd-rs`), which pulls in `hoomd-utility`; that crate hard-pins
`parquet = "58.0.0"` for an internal `ParquetLogger` convenience type, while our
own `molframe-interop` wants `parquet = "59"`. These are two incompatible majors, so a
build enabling both `traj` and `interop` compiles `parquet` (and its `parquet_derive`
proc-macro, and a duplicate `syn 2.0.119`) twice. `hoomd-utility`'s use of
`parquet` (`RecordWriter`, `SerializedFileWriter`, `WriterProperties`) is
unchanged between 58.4.0 and 59.2.0, so the fix is a version-bump request filed
upstream against `glotzerlab/hoomd-rs`, not a vendored/patched copy of
`hoomd-utility` in this repo — this repo carries no patched third-party crates,
full stop. In practice this only bites builds that deliberately opt into both
`traj` and `interop`; neither is in `molframe`'s `default` feature.

---

## 2. What lives here, what does not

| Path | In git? | Why |
|---|---|---|
| `crates/`, `Cargo.toml`, `Cargo.lock`, `README.md`, `CONTRIBUTING.md`, `RULES.md`, `rust-toolchain.toml` | yes | The published surface. |
| `docs/`, `idea/`, `inspo/` | no — `info/exclude` | Spec, design notes, and reference libraries. Re-publish deliberately, not by accident. |
| `HARNESS.md` | no — `info/exclude` | What the machine-local checks used to be, and how to re-adopt them. |
| `target/` | no — `info/exclude` | Build output. |

**No generated files and no local scripts.** Every capability is written by
hand; there is no codegen step to re-run. Nothing that belongs to one machine —
check scripts, git hooks, agent config — is tracked. Anything of that kind lives
outside the repo, listed in `.git/info/exclude`. No outside dependency is needed
for the loop.

Two artifacts are committed data rather than code, and stay:
`crates/molframe-chem/src/element_data.rs` (the 118-element table, consumed by
`use crate::element_data::ELEMENTS`) and `crates/molframe-xtal/data/space-groups.bin`
(the 530-Hall-setting catalogue, loaded by `include_bytes!`). Their provenance is
recorded in the `data/README.md` beside each; there is no script to regenerate
them from.

---

## 3. Scopes

The conventional-commit scope is a crate name or one of the cross-cutting
sections below. The list here is the definition: a scope in a commit subject
that is not a crate name must appear below.

```
core cif pdb geom cli py query spatial chem xtal surface interop traj compare
validate spec repo ci
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

Green means all of these pass:

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo test -p molframe --doc --features full

# The facade compiles with any single feature and with none. This is what keeps
# a gate on a module, a re-export or an enum variant matching the features that
# actually make the item exist (RULES §9).
#
# Every single feature, not a sample: the three pass-through ones at the end
# (`gzip`, `zstd`, `mmap`) gate nothing in `molframe` itself, but they change how
# `molframe-core` is built — `mmap` is what turns on its audited `unsafe` boundary —
# and nothing else in this list reaches that configuration. `full` is still
# covered, the same way it always was: `molframe-py` depends on it explicitly
# (`features = ["full"]`, `default-features = false`), so `cargo test --workspace`
# unifies the whole build back up to `full` regardless of what `molframe`'s own
# default is. `default` is no longer `full` — it's its own distinct combination
# (`mmcif` + `pdb`, not equal to any single entry in the loop below), so it gets
# its own explicit, bare check.
cargo check -p molframe
for f in "" pdb mmcif bcif modelcif geometry ic query spatial chemistry interop crystal \
         surface analysis validation sequence compare trajectory audit motif adapters \
         gzip zstd mmap; do
  if [ -z "$f" ]; then cargo check -p molframe --no-default-features || exit 1; \
  else cargo check -p molframe --no-default-features --features "$f" || exit 1; fi
done

# No unwrap or expect outside test code. A test that cannot unwrap is a test
# whose failure message got worse. Both exclusions matter: unit tests are
# `*_tests.rs`, integration and golden tests live under `tests/`.
grep -rnE '\.(unwrap|expect)(_[A-Za-z0-9_]+)?\(' crates/ --include="*.rs" \
  | grep -v '_tests.rs' | grep -v '/tests/'   # must be empty

# The file cap, .rs only. Nothing else in the tree is source.
find crates -name "*.rs" -print0 | xargs -0 wc -l | awk '$1>500 && $2 != "total" {print}'   # must be empty
```

`-D warnings` is the floor: the workspace permits exactly the lints declared in
`Cargo.toml`, and everything else is a regression.

No check reads `docs/` or a benchmark corpus, so a standalone clone passes the
whole list without either.

---

## 5. Current state

Every implemented capability sits at `β` (native + tested) rather than `✓`,
because a `✓` needs a row-specific golden workflow and differential evidence
against an external corpus that is not in this checkout.

Recently landed native, deterministic capabilities (all `β`, all with unit and
property tests): RMSF (`molframe-geom`); solvent-accessible surface and buried
surface (`molframe-surface`); contacts, contact maps, native contacts/Q, hydrogen
bonds, salt bridges, chain interfaces (`molframe-analysis`); clashes, cis-peptide
detection, occupancy/B-factor checks (`molframe-validate`); pairwise alignment
(`molframe-seq`); lDDT, TM-score, GDT-TS/HA (`molframe-compare`).

Outstanding deterministic work, roughly in dependency order: the rest of each of
those crates (Lee–Richards/SES/cavities; secondary structure, π-stacking, water
bridges; bond/angle deviation, planarity, chirality, rotamer/valence/stereo,
completeness; MSA, phylogenetics, k-mers, sequence formats, substitution
matrices; chain/atom mapping, DockQ/CAD/QS, CE alignment); `molframe-query`
altloc/entity/assembly and chirality selectors; PEOE charges in `molframe-chem`;
and the format breadth (trajectories, density maps, remaining structural
formats).

Track B (`molframe-audit`, `molframe-fx`) and `molframe-adapters` have landed
natively alongside the rest and are held to the same bar; what they still lack is
the row-specific golden workflow a `✓` would need.

**Where the golden and benchmark registers live.** They are the crates, not a
directory of manifests: the workflows are `crates/molframe/tests/golden_workflows.rs`
over `crates/molframe/tests/golden/`, the facade benchmarks are
`crates/molframe/benches/{facade,golden}.rs` over `crates/molframe/benches/golden/`,
and the fixtures and the measurement harnesses are `molframe-bench` and
`molframe-resource-bench`. Nothing under §2's excluded paths carries a register,
so a standalone clone runs the whole green list without one.

**An open question, recorded rather than built.** `PyGraph` holds its node and
edge tables as `Py<PyArray2<..>>` — zero-copy NumPy views, which is the format a
caller already consumes — and it exports no `__arrow_c_stream__`. A stream carries
one schema, and the node table (N×F) and the edge table (E×F′) are different
shapes; concatenating them is the O(n) copy that `ArrowStream` exists to avoid. If
a consumer appears that needs the stream, the shape of it should be decided then.
