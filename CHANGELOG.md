# Changelog

All notable changes to MolFrame are recorded here. The project follows
[Semantic Versioning](https://semver.org/); before 1.0 a minor release may
change the public API.

## Unreleased

### Fixed

- **Deposited secondary structure from mmCIF and BinaryCIF.** The direct
  readers dropped `struct_conf` and `struct_sheet_range` (and BinaryCIF also
  `atom_site_anisotrop`) before lowering, so every residue read from those
  formats had no helix or strand and the fallback assignment drew the whole
  cartoon. The HELIX and SHEET ranges are now read as the PDB reader already did.
- **Automatic secondary structure follows Kabsch–Sander.** The amide hydrogen
  was placed from the donor's own carbonyl with the energy's sign reversed, so
  almost no hydrogen bond formed and helices were inflated by a distance rule.
  Helices are now two overlapping n-turns, strands are ladders of consecutive
  bridges, proline and chain breaks have no amide hydrogen, and the first
  alternate location is read. Agreement with the deposited annotation is 80%
  on crambin, 84% on ubiquitin and haemoglobin and 99% on GroEL.

## 0.4.0

### Added

- **Python surface.** The namespaces that were placeholders now carry the Rust
  facade: `surface`, `compare`, `sequence`, `chemistry` (component-dictionary
  annotation, radii), `validation`, `spatial`, `trajectory` (reading, RMSD),
  `formats` (writers), `motif` (declarative functional-geometry evaluation), and
  the governed `analysis` kernels with `AnalysisPolicy`, `Analysis` and `Coverage`
  results, plus `Table` for tabular results.
- **Crystallography.** X-ray form factors and structure factors over the full
  space group with occupancy and isotropic or anisotropic displacement;
  `crystal.assembly`/`assemblies` as placements; resolution shells
  (`ResolutionBins`) and normalized amplitudes (`normalizers`); Niggli cell
  reduction (`reduce_cell`). Each is checked against Gemmi 0.7.5; see the
  evidence ledger for the agreement reached.
- A PDB `CRYST1` space group is kept and resolved to symmetry operations.
- `CHANGELOG.md`.

### Changed

- Selection by chain, residue name, residue number and similar predicates decides
  once per residue instead of once per atom (`chain A` on 58,870 atoms: 1.05 ms to
  0.014 ms), and thresholded `bfactor`/`occupancy` read each chunk directly.
- Default bond perception and secondary-structure assignment share one dense cell
  grid; reading 1AON through the facade drops from 114 ms to 40 ms.
- BinaryCIF atom-row access indexes columns by field position.
- Module roots in `molframe-py` and `molframe-traj` declare and re-export only.

### Fixed

- A conjunction of an atom-local and a non-local predicate could be reordered by
  cost and change its result.
- A build with only the `pdb` feature failed to compile.
- The `unwrap_or` family was banned by `RULES.md` but not enforced; six uses are
  now explicit matches, and `molframe-py` inherits the workspace lints.

### Known limits

- BinaryCIF read to a `Structure` reaches about 0.3 GiB/s of encoded input, below
  the 1 GB/s target; container parsing alone exceeds it.
- Radius-query and cell-list-build targets are not measured separately; a release
  Python import time is not measured.
- Anomalous dispersion, neutron and electron scattering tables, density
  calculation from a model and FFT are not implemented.
