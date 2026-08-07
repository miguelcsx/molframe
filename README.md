# pdbiox

**A batteries-included structural bioinformatics engine: semantics-preserving molecular data, high-performance analysis, and reproducible scientific evaluation.**

pdbiox reads what the field actually deposits — PDBx/mmCIF, BinaryCIF, ModelCIF, PDB-IHM, legacy PDB, and the trajectory and topology formats simulation produces — into a single columnar data model, and runs the whole workflow on it: selections, geometry, surfaces, contacts, assemblies, trajectories, structural comparison, validation, and machine-learning export.

It exists because the alternative is what everyone does today:

```python
structure  = Bio.PDB.MMCIFParser().get_structure(...)     # hierarchy, but slow
array      = biotite.structure.io.pdbx.get_structure(...) # arrays, but no crystallography
universe   = MDAnalysis.Universe(...)                     # selections and trajectories
st         = gemmi.read_structure(...)                    # assemblies and symmetry
traj       = mdtraj.load(...)                             # trajectory kernels
ent        = ost.io.LoadPDB(...)                          # complex comparison
```

Six data models, six sets of defaults, and a conversion between each pair. Every conversion is a place where altlocs get collapsed, `auth_` and `label_` identifiers get confused, insertion codes get dropped and assemblies get silently replaced by the asymmetric unit.

## Features

- **Formats that match the field.** PDBx/mmCIF, BinaryCIF, ModelCIF, PDB-IHM and legacy PDB, read losslessly and written without ever truncating data to fit a fixed column.
- **One columnar data model.** Chunked columnar storage with per-chunk statistics, dictionary-encoded identifiers, and hierarchy stored as offsets rather than pointers.
- **Zero-copy interop.** Coordinates live as `N × 3 f32` and hand straight to NumPy; annotations export as Apache Arrow without a copy.
- **Honest analysis contracts.** Every analysis takes an explicit `AnalysisPolicy` and returns its value *plus* its coverage, ambiguity status, warnings and full provenance — so a number is never silently produced from partial data.
- **Rust, Python and CLI** on the same engine. The Python package is a thin, typed binding over the Rust core; the CLI covers `info`, `convert`, `validate`, `measure`, `rmsd` and `policy`.

## Getting started

### Rust

Add `pdbiox` to your `Cargo.toml`, then:

```rust
use pdbiox::prelude::*;

let (structure, findings) = pdbiox::read_with_diagnostics("1abc.pdb")?;

println!("{} atoms in {} chains", structure.atom_count(), structure.chain_count());
for finding in &findings {
    eprintln!("{}", Rendered::new(finding));
}
```

`read_with_diagnostics` is the honest form: the findings say what was wrong with the file, and a file that parses is not the same as a file that is right. `pdbiox::read` discards them for convenience.

### Python

```python
import pdbiox

structure = pdbiox.read("1abc.pdb")
xyz = structure.xyz()                      # N x 3 float32, zero-copy NumPy view

for chain in structure.chains():
    print(chain.label(), len(chain.residues()))
```

### Command line

```bash
pdbiox info     1ubq.pdb
pdbiox convert  1ubq.pdb  out.cif
pdbiox validate 1ubq.pdb
pdbiox measure  1ubq.pdb    # extent, radius of gyration, centre
pdbiox rmsd     model1.pdb model2.pdb
```

## The three claims

**1. Coverage.** pdbiox aims to be the single library that does everything the reference libraries do — everything, tracked row by row, with adapters as a staging state and a scheduled replacement, not a destination.

**2. A representation that makes coverage fast.** Chunked columnar storage with per-chunk statistics, dictionary-encoded identifiers, hierarchy as offsets rather than pointers, adaptive selection sets, and late materialization. Coordinates live as `N × 3 f32` and hand straight to NumPy without a copy; annotations export as Arrow without a copy.

**3. Analysis contracts.** This is the part no other library has. Every analysis takes an explicit `AnalysisPolicy` and returns a value *plus* its coverage, its ambiguity status, its warnings and its full provenance. Decisions that are currently silent defaults — which assembly, which model, which conformer, which identifier namespace, what to do about missing atoms — become data you can inspect, hash, vary and publish.

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

Phase 1 of the specification is implemented across the core crates: storage, topology, structure views and selections, diagnostics and analysis contracts, PDBx/mmCIF and legacy-PDB readers and writers, geometry, internal coordinates, the CLI and the Python bindings. The remaining Phase 1 work is evidence: full-corpus differential classification and row-specific golden workflows. See `CONTRIBUTING.md` and `Rules.md` for how to work on it.

## Repository layout

```
crates/     the workspace: pdbiox-core, pdbiox-cif, pdbiox-pdb, pdbiox-geom,
            pdbiox-ic, pdbiox-cli and the pdbiox facade
python/     Python bindings (PyO3 + NumPy zero-copy views) and typed stubs
scripts/    verification, data generation and Python boundary checks
LICENSES/   licence texts (MIT, Apache-2.0, BSD-3-Clause)
REUSE.toml  machine-readable licence annotations (REUSE-compliant)
CITATION.cff  how to cite pdbiox
```

## Naming

`pdbiox` throughout: the crate, the PyPI distribution, the Python import, the CLI and the GitHub repository. Earlier drafts used `pdbio` and `pbiox`; both are retired. `pdbio` in particular is unusable — it is taken on PyPI by an unrelated package and is also the name of a Biopython class.

```
crates:  pdbiox, pdbiox-core, pdbiox-cif, ...
PyPI:    pdbiox
python:  import pdbiox
native:  pdbiox._native
CLI:     pdbiox
```

## Development

```bash
cargo build --workspace
cargo test  --workspace
./scripts/verify.sh   # the single definition of "green"
```

## License and citation

pdbiox is dual-licensed under **MIT OR Apache-2.0** — pick either. Both are permissive: you may use pdbiox in private, commercial or proprietary work, and you are never required to open-source the code that uses it. What the licences require is that you retain the copyright and license notices on any copies or substantial portions of the Software.

If you use pdbiox in research, please cite it. Full machine-readable metadata lives in [`CITATION.cff`](CITATION.cff):

```yaml
title: "pdbiox"
authors:
  - family-names: "Cárdenas"
    given-names: "Miguel"
license: "MIT OR Apache-2.0"
```

The repository is [REUSE](https://reuse.software)-compliant: every file carries machine-readable licence information (see `LICENSES/` and `REUSE.toml`). The bundled reference data keeps the licences of its upstream sources — ionic radii from [mendeleev](https://github.com/lmmentel/mendeleev) (MIT) and the space-group catalogue from [spglib](https://spglib.github.io/spglib/) (BSD-3-Clause).
