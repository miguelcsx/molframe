# pdbiox

**A batteries-included structural bioinformatics engine: semantics-preserving molecular data, high-performance analysis, and reproducible scientific evaluation.**

[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](#license-and-citation)
[![Built with Rust](https://img.shields.io/badge/built%20with-Rust-orange.svg)](https://www.rust-lang.org)
[![REUSE compliant](https://img.shields.io/badge/REUSE-compliant-green.svg)](https://reuse.software)

pdbiox reads structural and simulation data into a single columnar data model
and runs selections, geometry, surfaces, contacts, assemblies, trajectories,
structural comparison, validation and machine-learning export over that model.
The live format inventory is tracked row by row in `docs/PARITY.md`; formats not
yet marked `β` there are commitments, not shipped claims.

## Why pdbiox

Today the same analysis is stitched together from libraries that each own one slice of the problem:

```python
structure = Bio.PDB.MMCIFParser().get_structure(...)     # hierarchy, but slow
array     = biotite.structure.io.pdbx.get_structure(...) # arrays, but no crystallography
universe  = MDAnalysis.Universe(...)                     # selections and trajectories
st        = gemmi.read_structure(...)                    # assemblies and symmetry
traj      = mdtraj.load(...)                             # trajectory kernels
ent       = ost.io.LoadPDB(...)                          # complex comparison
```

Six data models, six sets of defaults, and a conversion between each pair. Every conversion is a place where altlocs get collapsed, `auth_` and `label_` identifiers get confused, insertion codes get dropped, and assemblies get silently replaced by the asymmetric unit. pdbiox is one engine that keeps those semantics intact from parse to result.

## Features

- **Formats that match the field.** PDBx/mmCIF, BinaryCIF, legacy PDB, PQR and
  PDBQT are wired through the Rust facade; document-preserving CIF keeps unknown
  extension categories losslessly, and fixed-column writers refuse truncation.
- **One columnar data model.** Chunked columnar storage with per-chunk statistics, dictionary-encoded identifiers, and hierarchy stored as offsets rather than pointers.
- **Zero-copy interop.** Coordinates live as `N × 3 f32` and hand straight to NumPy; annotations export as Apache Arrow without a copy.
- **Honest analysis contracts.** `AnalysisPolicy`, `Analysis<T>`, coverage,
  ambiguity status and provenance are native core types. Uniformly wrapping
  every public kernel in that contract is still a measured cross-cutting gate;
  direct kernels do not silently acquire an undocumented policy.
- **Rust, Python, and CLI on one engine.** The Python package is a thin, typed binding over the Rust core. The CLI exposes inspection, conversion, selection, geometry, surfaces, validation, comparison, trajectories, system containers, audit, functional evaluation and deterministic batch execution without duplicating kernels.

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

`read_with_diagnostics` is the honest form: the findings say what was wrong with the file, because a file that parses is not the same as a file that is right. `pdbiox::read` discards them for convenience.

### Python

```python
import pdbiox

structure = pdbiox.read("1abc.pdb")
xyz = structure.xyz                        # N x 3 float32, zero-copy NumPy view

for index in range(len(structure.chains)):
    chain = structure.chains[index]
    print(chain.label, len(chain.residues))
```

### Command line

```bash
pdbiox info     1ubq.pdb
pdbiox convert  1ubq.pdb  out.cif
pdbiox validate 1ubq.pdb
pdbiox measure  1ubq.pdb    # extent, radius of gyration, centre
pdbiox rmsd     model1.pdb model2.pdb
```

## Analysis contracts

This is the part no other library has. The core can carry an explicit
`AnalysisPolicy` and return a value together with its coverage, ambiguity
status, warnings, and full provenance. Decisions that are usually silent
defaults (which assembly, which model, which conformer, which identifier
namespace, what to do about missing atoms) become data you can inspect, hash,
vary, and publish. The repository does not claim this criterion complete until
every public analysis entry point is wired to the enriched result.

`Indeterminate` is a real answer. A library that cannot decline to produce a number will always produce a defensible-looking wrong one.

## Status

All 302 deterministic product rows in `docs/PARITY.md` have native
implementations and tests at `β`. The only phased `·` is the external PDB-scale
audit study, which is evidence work rather than a library capability. Promotion
to `✓` still requires the row-specific golden and differential evidence defined
by the specification; `β` is not a claim that this evidence already exists.

## Repository layout

```
crates/       the Rust workspace (pdbiox facade plus the format, geometry,
              chemistry, analysis, fx, audit, ModelCIF, binding and CLI crates)
python/       Python bindings (PyO3 + NumPy zero-copy views) and typed stubs
scripts/      verification, data generation, and Python boundary checks
LICENSES/     license texts (MIT, Apache-2.0, BSD-3-Clause)
REUSE.toml    machine-readable license annotations (REUSE-compliant)
CITATION.cff  how to cite pdbiox
```

## Contributing

```bash
cargo build --workspace
cargo test  --workspace
./scripts/verify.sh   # the single definition of "green"
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the workflow and [RULES.md](RULES.md) for the binding style guide.

## License and citation

pdbiox is dual-licensed under **MIT OR Apache-2.0**; pick either. Both are permissive: you may use pdbiox in private, commercial, or proprietary work, and you are never required to open-source the code that uses it. What the licenses require is that you retain the copyright and license notices on any copies or substantial portions of the Software.

If you use pdbiox in research, please cite it. Full machine-readable metadata lives in [CITATION.cff](CITATION.cff):

```yaml
title: "pdbiox"
authors:
  - family-names: "Cárdenas"
    given-names: "Miguel"
license: "MIT OR Apache-2.0"
```

The repository is [REUSE](https://reuse.software)-compliant: every file carries machine-readable license information (see `LICENSES/` and `REUSE.toml`). Bundled reference data keeps the licenses of its upstream sources: ionic radii from [mendeleev](https://github.com/lmmentel/mendeleev) (MIT) and the space-group catalogue from [spglib](https://spglib.github.io/spglib/) (BSD-3-Clause).
