# pdbiox

**A batteries-included structural bioinformatics engine: semantics-preserving molecular data, high-performance analysis, and reproducible scientific evaluation.**

pdbiox is one library for structural bioinformatics. It reads what the field actually deposits — PDBx/mmCIF, BinaryCIF, ModelCIF, PDB-IHM, legacy PDB, and the trajectory and topology formats simulation produces — into a single columnar data model, and runs the whole workflow on it: selections, geometry, surfaces, contacts, assemblies, trajectories, structural comparison, validation, and machine-learning export.

It exists because the alternative is what everyone does today:

```python
structure  = Bio.PDB.MMCIFParser().get_structure(...)   # hierarchy, but slow
array      = biotite.structure.io.pdbx.get_structure(...)  # arrays, but no crystallography
universe   = MDAnalysis.Universe(...)                   # selections and trajectories
st         = gemmi.read_structure(...)                  # assemblies and symmetry
traj       = mdtraj.load(...)                           # trajectory kernels
ent        = ost.io.LoadPDB(...)                        # complex comparison
```

Six data models, six sets of defaults, and a conversion between each pair. Every conversion is a place where altlocs get collapsed, `auth_` and `label_` identifiers get confused, insertion codes get dropped and assemblies get silently replaced by the asymmetric unit.

## The three claims

**1. Coverage.** pdbiox aims to be the single library that does everything the six reference libraries do. Not a subset of "high-value workflows" — everything, tracked row by row in [`docs/PARITY.md`](docs/PARITY.md). Adapters are a staging state with a scheduled replacement, not a destination.

**2. A representation that makes coverage fast.** Chunked columnar storage with per-chunk statistics, dictionary-encoded identifiers, hierarchy as offsets rather than pointers, adaptive selection sets, and late materialization. Coordinates live as `N × 3 f32` and hand straight to NumPy without a copy; annotations export as Arrow without a copy. See [`docs/04-storage-layout.md`](docs/04-storage-layout.md).

**3. Analysis contracts.** This is the part no other library has. Every analysis takes an explicit `AnalysisPolicy` and returns a value *plus* its coverage, its ambiguity status, its warnings and its full provenance. Decisions that are currently silent defaults — which assembly, which model, which conformer, which identifier namespace, what to do about missing atoms — become data you can inspect, hash, vary and publish. See [`docs/10-analysis-contracts.md`](docs/10-analysis-contracts.md).

```python
result = structure.contacts("chain A", "chain B", cutoff=4.5, policy=policy)

result.value        # the contacts
result.status       # Complete | Partial | Ambiguous | Indeterminate
result.coverage     # what fraction of the intended atoms were actually available
result.warnings     # what pdbiox had to decide for you
result.provenance   # input hash, policy hash, CCD version, dictionary version, ...
```

`Indeterminate` is a real answer. A library that cannot decline to produce a number will always produce a defensible-looking wrong one.

## Status

**Phase 1, in progress.** The specification is complete and `docs/` remains normative — it is what gets built. Implementation has started:

```
crates/pdbiox-core   storage, topology, structure, views, selections,
                     diagnostics, AnalysisPolicy, Analysis<T>, provenance
crates/pdbiox-cif    PDBx/mmCIF: lexer, Document, lowering, both writers
crates/pdbiox-geom   distances, angles, torsions, moments, superposition
crates/pdbiox-pdb    legacy PDB reader and refusing writer
crates/pdbiox        facade
crates/pdbiox-cli    info · convert · validate · measure · rmsd · policy
```

```bash
cargo run -p pdbiox-cli -- info    inspo/pdbtbx/example-pdbs/1ubq.cif
cargo run -p pdbiox-cli -- convert inspo/pdbtbx/example-pdbs/1ubq.pdb out.cif
cargo run -p pdbiox-cli -- measure inspo/pdbtbx/example-pdbs/1ubq.pdb
```

`pdbiox-py` and the rest of phase 1 are outstanding. What is and is not met is tabulated in [`docs/ROADMAP.md`](docs/ROADMAP.md) §Phase 1; per-capability status is the `pdbiox` column of [`docs/PARITY.md`](docs/PARITY.md), where `β` means implemented and tested but not yet backed by a golden workflow. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for how to work on it.

## Documentation map

| | |
|---|---|
| **Start here** | [`00-vision-and-scope.md`](docs/00-vision-and-scope.md) · [`GLOSSARY.md`](docs/GLOSSARY.md) · [`01-requirements.md`](docs/01-requirements.md) |
| **Architecture** | [`02-architecture.md`](docs/02-architecture.md) · [`03-data-model.md`](docs/03-data-model.md) · [`04-storage-layout.md`](docs/04-storage-layout.md) · [`05-identifiers.md`](docs/05-identifiers.md) |
| **Subsystems** | [`06-formats`](docs/06-formats.md) · [`07-chemistry`](docs/07-chemistry.md) · [`08-selection-language`](docs/08-selection-language.md) · [`09-spatial-and-geometry`](docs/09-spatial-and-geometry.md) · [`10-analysis-contracts`](docs/10-analysis-contracts.md) · [`11-assemblies-symmetry`](docs/11-assemblies-symmetry.md) · [`12-trajectories`](docs/12-trajectories.md) · [`13-comparison-validation`](docs/13-comparison-validation.md) · [`14-interop`](docs/14-interop.md) · [`24-sequence`](docs/24-sequence.md) · [`25-analysis-catalogue`](docs/25-analysis-catalogue.md) |
| **Surfaces** | [`15-api-rust`](docs/15-api-rust.md) · [`16-api-python`](docs/16-api-python.md) · [`17-cli`](docs/17-cli.md) · [`18-errors-and-diagnostics`](docs/18-errors-and-diagnostics.md) · [`19-performance`](docs/19-performance.md) |
| **Process** | [`20-testing`](docs/20-testing.md) · [`21-benchmarks`](docs/21-benchmarks.md) · [`22-research-programme`](docs/22-research-programme.md) · [`23-governance`](docs/23-governance.md) · [`PARITY.md`](docs/PARITY.md) · [`ROADMAP.md`](docs/ROADMAP.md) |
| **Decisions** | [`docs/adr/`](docs/adr/) — fourteen accepted records, each with its rejected alternatives |

## Repository layout

```
pdbiox/
├── docs/          normative specification — this is what gets built
├── idea/          historical design deliberation (Spanish, NOT normative)
└── inspo/         reference implementations, for design study and test oracles
```

`idea/` is the original brainstorming that produced this project. It is kept for provenance and is explicitly **not** authoritative: where `docs/` and `idea/` disagree, `docs/` wins.

`inspo/` holds clones of Biopython, Biotite, MDAnalysis and pdbtbx. They are read for *design* and used as differential-test oracles. **Check the licence before adapting any code** — MDAnalysis is GPL-2.0-or-later and its code must not enter pdbiox. See [`23-governance.md`](docs/23-governance.md).

## Naming

`pdbiox` throughout: the crate, the PyPI distribution, the Python import, the CLI and the GitHub repository. Earlier drafts used `pdbio` and `pbiox`; both are retired. `pdbio` in particular is unusable — it is taken on PyPI by an unrelated package and is also the name of a Biopython class.

```
crates:  pdbiox, pdbiox-core, pdbiox-cif, ...
PyPI:    pdbiox
python:  import pdbiox
native:  pdbiox._native
CLI:     pdbiox
```

## Licence

Dual MIT / Apache-2.0.
