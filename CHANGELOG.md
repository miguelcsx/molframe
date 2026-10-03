# Changelog

All notable changes to MolFrame are recorded here. The project follows
[Semantic Versioning](https://semver.org/); before 1.0 a minor release may
change the public API.

## Unreleased

### Added

- Secondary query macros, typed Rust/Python selectors, the public Python
  SecondaryStructure vocabulary and native DSSP table expose all modern states.
  Canonical CIF/PDB secondary records round-trip PPII class 10, named helices,
  other helices and strands without collapsing unknown classes into alpha.
  PDB ranges resolve insertion-aware negative/hybrid-36 endpoint identities in
  chain order and isolate model-local annotations. Canonical secondary output
  refuses ambiguous or unrepresentable endpoint identities before streaming;
  CIF label ranges never match author aliases.
- Structure-aligned partial charges prefer complete finite file columns or use
  versioned CCD PEOE chemistry with explicit linkage and hydrogen projection.
  Source provenance is retained; unresolved chemistry fails explicitly.
- Affine screened contact-potential grids in kT/e at 298 K, with dielectric 4r,
  sub-ångström softening, finite 12 Å cutoff and shared spatial/worker execution.
- MRC descriptor/density voxel-to-world affine matrices share the triclinic
  sampling and origin conventions of the existing Cartesian sampler.
- Assembly covalent links use the canonical bond-perception predicate and
  stable chain-instance IDs; metal coordination and same-instance pairs are
  excluded. Rust and Python expose the resulting links.

- Carbohydrate chemistry: curated MIT-attributed Mol* SNFG metadata, topology-only
  five/six-member saccharide rings, finite ring geometry, and provenance-preserving
  glycosidic/protein attachments. Optional vacant-site inference is bounded to
  2 Å and refuses competing sites; Rust and Python use the chemistry namespace.
  Bundled 1HZH/CCD fixtures check every sugar and all 16 deposited branch links;
  the two root NAG–ASN contacts exceed 2 Å and are correctly left unlinked.

- Hydrogen query predicates `polar_hydrogen` and `nonpolar_hydrogen`, with
  Rust typed builders, Python `sel` constructors and completion. Polarity
  follows full bond topology even in hydrogen-only views: any N/O/S neighbour
  is polar; all other hydrogen, including unbonded H, is nonpolar. Unavailable
  connectivity reports `E4003` instead of classifying by guesswork.

- **Comparison across renamed chains.** `mapped_dockq` and `mapped_qs_score`
  (and `CompareExt` methods) match chains by sequence, residues by alignment and
  atoms by name, with CCD-equivalent atoms permuted to the nearest fit over the
  atoms both residues contain. Chains of equal sequence identity, such as the
  halves of a homodimer, are resolved by scoring every assignment of the
  requested pair. The CLI runs it as `compare --map-chains` with every control
  explicit. On the haemoglobin alpha/beta chains of 4HHB, renaming the chains
  leaves DockQ bit-identical.
- **Reads.** Category filters (`ReadOptions.categories`, CLI `--skip-category`);
  SDF, MOL2 and core CIF as one-residue structures, with SDF output; typed PDB
  headers and `REMARK 350` assemblies; the ten-state secondary-structure
  vocabulary with its source; an `assembly N` selector.
- **Superposition.** Shape parameters and a quaternion characteristic-polynomial
  (QCP) superposition with a Jacobi fallback. Its benchmark is not yet recorded.
- **Reference inputs.** `read_reference_library`, `read_rotamer_profile`,
  `read_plane_restraints` and `read_tls_groups` load the data the validation
  kernels refuse to carry, from JSON or TOML with strict schemas and typed
  errors; a selection that matches nothing is refused.
- **Command line.** `interactions`, `contact-map`, `native-contacts`, `nma`,
  `chem peoe|smarts`, `crystal mtz-info|map-stats|map-correlation|mates` and
  `trajectory rmsf` are new; `validate` gains valence, ligand, chirality,
  cis-peptide, nucleic, Ramachandran, rotamer, plane-restraint and TLS checks,
  and `compare` gains QS, CAD, CE, contact similarity and interface, pocket and
  ligand RMSD. Every scientific parameter is a required flag.
- `named_ligand_rmsd` compares ligands by atom name over the observed fragment
  with chemical symmetry; `GaussianNetworkModel::fluctuations` and
  `CrystalNeighbor::is_symmetry_mate` are new.
- `sse` requires a caller-authored polymer role profile; the CLI never guesses
  backbone atom names.
- The Python surface is gated by ruff (all rules), pyright (strict) and
  `pyright --verifytypes`.

### Changed

- Automatic enrichment and explicit secondary-structure analysis share one native
  DSSP 4 classifier, including cross-chain sheets, beta-bulge spans, peptide-break
  guards, donor top-two energies and stretch-level helix precedence. Genuine
  `PolyProline` appends stable code 10 without changing `OtherHelix` or
  `BetaBridge`. PPII uses three consecutive phi/psi windows and preserves turns
  and bends. Exact-state comparison with mkdssp 4.5.0 matched all 8015 1AON and
  1544 7QPD reference residues; native operation requires no external executable.

- **Read and analysis throughput.** BinaryCIF to `Structure` on 1AON went from
  13.96 ms to 7.53 ms (369 MiB/s) and mmCIF to `Structure` on 4HHB from 3.55 ms
  to 3.20 ms (228 MiB/s): rows are fed to one monomorphised lowerer, dictionary
  strings are interned once per column, `struct_conn` endpoints are resolved
  only in the residues they name, and the lexer tracks offsets instead of
  positions. Default perception on 1AON fell from 31.3 ms to 11.0 ms (parallel
  bond search); `element ZN` over 100,000 atoms from 0.97 ms to 105 ns (chunk
  pruning); Shrake–Rupley over 100,000 atoms at 960 points from 233 ms to
  122 ms, bit-identical, through one streaming engine with SIMD occlusion.
- **Known limit.** BinaryCIF reading remains below the 1 GB/s specification
  target (369 MiB/s on 1AON); the remaining cost is per-row lowering. The
  evidence ledger lists the top frames and the attempts that were reverted.
- `molframe::perceive_in` and `molframe_chem::perceive_bonds_in` take an
  `ExecutionContext`; `FrameTransform::apply` takes `&mut self`, and
  `ChainedReader`/`PipelineReader` support bounded reads.
- A distance matrix under a memory budget now needs one matrix of headroom
  rather than two.

### Fixed

- Atom-selection materialization compacts secondary states and their source
  provenance with the residue hierarchy, so deleting earlier residues cannot
  change a retained helix class. Unavailable columns remain unavailable, empty
  selections discard all residue rows, and dense models retain one shared
  compacted hierarchy rather than duplicating its residues per frame.

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
