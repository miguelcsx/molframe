# Changelog

All notable changes to MolFrame are recorded here. The project follows
[Semantic Versioning](https://semver.org/); before 1.0 a minor release may
change the public API.

## 0.6.1 — 2026-10-15

- **PyPI now carries the 0.6.0 content.** The 0.6.0 wheels uploaded on
  2026-10-10 predate the changes below and the crate published on
  2026-10-15, so a registry version cannot be reused; this release ships
  them together to both registries.
- **`compare.lddt` is now a real, residue-aware lDDT; the coordinate kernel
  is `compare.distance_difference`.** The old `lddt(model, reference)` took
  bare coordinates, so it counted distances inside a residue, could not
  recognise a model that named the two carboxylate oxygens the other way
  round, and could not count a missing atom against the model. It is renamed
  `distance_difference` (Rust `distance_difference`, `DistanceDifferenceOptions`;
  its governed provenance name is `distance-difference`) and documented as a
  coordinate kernel. The new `lddt(model, reference, *, components,
  components_version, ...)` scores structures: reference distances within
  15 A between atoms of different residues, tolerances 0.5/1/2/4 A, a missing
  model atom or residue counts as not preserved, each pair is tallied from
  both of its atoms, and per-residue scores are reported. Interchangeable
  atoms (ASP/GLU oxygens, PHE/TYR ring carbons, ARG nitrogens, phosphate
  oxygens) are renamed to the better naming using the bond graph with bond
  orders ignored and a planar-centre guard, so the two methyls of valine or
  leucine are never swapped. Residues are paired by chain and sequence
  position; `symmetry=False` scores names as given. The result carries the
  parameters and conventions it was computed under, and Rust has
  `governed_residue_lddt`.
- **A wrong array dtype or rank is a clear `MolframeTypeError`.** Passing a
  `float64` array where `float32` is expected used to fail with
  `'ndarray' object is not an instance of 'ndarray'`, because PyO3's own
  check can only name the type, not the dtype. Every array parameter now goes
  through one extraction that states the dtype and rank expected and received
  (`expected a float32 array, received dtype float64; convert it once with
  array.astype(numpy.float32)`) and names the parameter in the exception
  note. Nothing is converted silently, since the copy is the size of the
  input and would narrow the values.
- **A failed native BCIF encoding reports its real cause.** The structure
  capsule (ABI version 3) appends a failure callback returning the
  diagnostic number and message, so a consumer such as molgfx raises the
  matching `molframe.errors` exception with its `code` instead of a generic
  `ValueError`. Consumers still accept version 2 producers, which have no
  failure channel.
- **The k-d tree is parallel and exact, and `spatial.knn` is first class.**
  The k-d tree is now an implicit, bucketed tree built and queried in
  parallel; on degenerate layouts (needle, collinear, offset, 1e6 A outliers)
  it no longer degrades, and periodic search visits images nearest first
  instead of sweeping 125 shifts. `molframe.spatial.knn(points, k, query=None,
  max_distance=None, cell=None, exclude_self=None)` returns `(n_query, k)`
  `uint32` indices and `float32` distances ordered by distance then index,
  padding with `4294967295` and `inf`. `interop.graph(edges="k_nearest")` uses
  it and honours `backend`. `neighbor_pairs` and `cross_pairs` hand NumPy
  columns over without extra copies, `neighbor_pairs_chunks` yields bounded
  chunks, and `auto` moves to the tree on crowded clouds. Pair sets are
  unchanged; distances are float32, so a pair within float32 rounding of the
  cutoff may differ from a float64 library.
- **Selection no longer holds the GIL, and `name CA` is no longer slow.**
  `Structure.select` and `Query.select` run through the execution context: the
  GIL is released for the native evaluation, Ctrl-C raises `KeyboardInterrupt`,
  and both accept `context=` so `ExecutionContext.cancel()` stops them.
  `resolve_altlocs` and `Selection.to_coordinates` also release the GIL.
  Single-predicate selections on a per-atom symbol column (`name`, `auth_name`,
  `altloc`) compare symbol ids chunk by chunk instead of resolving the chunk,
  residue and chain of every atom, which took `name CA` on a 2.4M-atom
  structure from about 0.9 s to a few tens of milliseconds. Results are
  unchanged.
- **Input limits are bounded by default.** `Limits::default()` now caps
  decompressed size at 16 GiB, declared rows per category at 2^28 and
  dictionary entries at 2^24 (previously unbounded). `Limits::unbounded()`
  restores the old behaviour for trusted input; in Python pass a larger
  `max_*` to `ReadOptions`. The compression-ratio and nesting ceilings are
  unchanged.

- **Pinned differential evidence against three reference libraries.** GW-041 to
  GW-045 promote the crystallography, cell-transform, TM-score and element-data
  rows using three new files under `crates/molframe/tests/reference/`: all 530
  Hall settings against gemmi 0.7.5, the atomic weights against the CIAAW 2021
  conventional table (gemmi still carries pre-2021 Ge and Se), and the TM-score
  against TM-align on the 1UBQ C-alpha trace. The TM reference records that the
  fixed-correspondence fit reads about 0.04 lower than TM-align on a hinge,
  because TM-align's search can also re-pair residues.

### Known limits found while validating

- `molframe.spatial.neighbor_pairs` and `cross_pairs` take no unit cell, so the
  Python surface cannot ask for a periodic neighbour search even though
  `molframe-spatial` implements the triclinic minimum image and the Rust
  `SpatialSearchOptions` carry it. The periodic path is exercised only from Rust.
  Closing this means adding a `cell=` argument to the two bindings, which changes
  a published signature and is left for a deliberate API decision.

### Added

- `molframe_geom::superpose_weighted`: the weighted quaternion rigid fit, the
  minimiser of a per-pair weighted sum of squared deviations. An empty weight
  slice reproduces `superpose`; a weight of zero removes a pair from the fit.

### Fixed

- **A component coordinate that cannot be represented as `f32` saturates rather
  than becoming an infinity.** `small_cif_to_structure` narrowed each site with
  an unchecked cast, so a magnitude beyond the `f32` range became an infinity
  and every distance, angle and fit downstream inherited it. The new
  `molframe_chem::numeric::f64_to_f32` saturates into the finite range, carries a
  recorded NaN or infinity through as itself, and reserves saturation for a
  finite value that genuinely exceeded the type.
- **QS-score now uses the published normalisation.** The numerator is
  `w(min(d_ref, d_mod)) * (1 - |d_ref - d_mod| / 12)` summed over shared residue
  pairs, and the denominator is that shared weight plus the weight of every
  contact present on one side only. The previous squared-weight-mass ratio
  scored a contact moved 6 Å as 0.130 where the published score is 0.5, and
  penalised an extra model contact by 1/3 where the published score is 1/2.
  The published contact test is inclusive at the cutoff, so a pair exactly at
  12 Å now counts.
- **Python QS defaults moved from 5 Å to the published 12 Å** contact cutoff in
  `qs_score` and `mapped_qs_score`. The `DockQ` 5 Å default is a different
  metric and is unchanged. A non-default `contact_distance` selects a
  nonpublished variant: it changes which pairs are contacts, while the 12 Å
  distance penalty stays tied to the definition.
- **Deposited connection classifications survive storage and export.**
  `BondRecord` carries the source's `_struct_conn.conn_type_id` as an interned
  identifier, and the canonical mmCIF and BinaryCIF writers emit each bond's
  own classification instead of one caller-supplied type. 4HHB's four `metalc`
  Fe–histidine connections are written as `metalc`, not as `covale`. A bond
  whose source recorded no classification is written `?` rather than asserted
  to be covalent. The global `CifWriteOptions::with_connection_type_id`, its
  error and the CLI `--cif-connection-type` flag are gone; the kind is remapped
  across dictionaries by `Structure::merge` and the adapter and Arrow bond
  exports carry it as a nullable UTF-8 column.
- **The connection class survives the provider bond projection.** The chunked
  bond payload carried only the order and the provenance, so a consumer reading a
  bond across a chunk boundary could not tell a metal coordination from a
  covalent link — the distinction the class column exists for. `BondChunkRecord`
  now carries the interned class as well.
- **A semantic structure difference compares bond classes by their meaning.**
  `changed_bonds` compared the raw `BondRecord`, whose connection class is an
  interned identifier numbered per dictionary. Two independently read copies of
  one file numbered the same class differently and were reported as different,
  and two files whose classes happened to land in the same slot were reported as
  identical. The comparison now resolves the class against each structure's own
  dictionary, as it already did for entity and chain fields.
- **A ragged PDB ensemble keeps its atoms reachable from Python.** Such an
  ensemble holds its models as independent structures, so the outer object's
  counts are zero — and the binding exposed no way to reach the models, meaning a
  Python caller saw a file of several hundred atoms as an empty molecule with no
  diagnostic. `Structure.is_ragged_ensemble` and `Structure.ragged_models` are now
  bound, and the facade gained the matching accessors.
- **Half-sphere exposure reports coverage over its own domain.** The governed
  kernel counted every residue as intended, so a structure with bound water
  reported an analysis that measured all 70 of its protein residues as 70 of 134,
  or `partial`, and `missing_atoms` would reject a complete answer. The domain is
  now the residues the role profile identifies as protein, read from the roles
  rather than from the chain kind, because a PDB read has no `_entity_poly.type`
  to classify a chain from.
- **A coordinate row that will not parse is not reported as an empty file.** A
  serial that outgrew its five columns pushes a digit into the record name, so
  `ATOM 100000` reads `ATOM 1` and matched no record at all; the file was then
  refused with "the file contains no coordinate records", which describes a
  different problem and leaves the caller nothing to act on. Such rows are now
  counted and reported as `E1201`, a field-layout fault, while a file that truly
  holds no coordinate rows still reports `E1001`.
- **A connection to another crystal image is refused instead of flattened.**
  `_struct_conn` lowering preflights `ptnr1_symmetry` and `ptnr2_symmetry`; a
  partner named through a nonzero cell translation, or through an operator the
  block's own operator table does not define as the identity, is reported as
  `E3016` (breaking) rather than lowered to a bond between the atoms that share
  those identifiers. Lossless `Document` parsing is unaffected. The
  `space_group_symop` and `symmetry_equiv` categories are retained whenever
  `struct_conn` is, so a category filter cannot turn an explicit nonidentity
  operator into the implicit local-identity assumption.
- **TM-score now searches superpositions instead of taking one least-squares
  fit.** The published score is a maximum over superpositions, and the TM sum at
  a single least-squares fit is a smaller number whenever the movement is not
  rigid: on the 1UBQ C-alpha trace, a least-squares fit scored 0.2169 where the
  reweighted search scores 0.5237, which also beats every single-window fit of
  the same trace (0.5236). The fit is now the deterministic reweighted search
  whose weights are the squared TM kernel, keeping the best score seen.
  `molframe-geom` gains `superpose_weighted` for the weighted quaternion fit it
  is built from. The search stays a local one and fits a fixed correspondence,
  so it remains below TM-align on motion that swings one part of a chain past
  another (0.524 against 0.568 on a half-and-half hinge, 0.665 against 0.680 on a
  displaced tail, exact for a rigid displacement); that limit is recorded in
  `tests/reference/tm_score_tmalign_0.3.0.json` rather than papered over.
- `molframe_bench::Sample::bcif` returns `Option<&'static [u8]>`: the 2M7C
  ensemble fixture ships PDB, and the accessor previously returned those bytes
  as BinaryCIF, which made a format sweep that trusted it measure PDB parsing
  under a BinaryCIF label.
- **Half-sphere exposure no longer returns a silent empty table.** A residue needs
  both an alpha and a beta carbon role — the beta carbon fixes the outward
  direction — so a residue without one takes no part. The governed kernel reported
  an empty result as `Complete`, which reads as "no residue is crowded" when the
  truth is "nothing was measured"; a caller whose profile supplied only backbone
  roles saw zero rows and no explanation. The skipped residues are now reported as
  missing coverage, so `missing_atoms` decides what a partial answer means, and
  the descriptor requires a side-chain direction, so an input with none is
  indeterminate rather than empty. Measured on a three-residue fixture with one
  glycine: `used=2, missing=1, intended=3`, quality `Partial`.

## 0.5.1 — 2026-10-05

- Synchronize MOL/SDF formal charges with chemistry annotations; omitted MDL charges are known neutral.
- Annotate explicit aromatic bonds and perceive conjugated Hückel circuits conservatively, preserving original bond orders.
- Match aromatic SMARTS consistently across aromatic and Kekulé SDF encodings, including implicit neutral nitrogen hydrogens.
- Read unversioned fixed-column SDF records produced by RCSB ModelServer.
- Preserve heavy atoms with hydrogen-like CCD names during annotation.
- Resolve alternate conformers and merge annotated snapshots through Python.
- Preserve molecular identifiers through native BinaryCIF transport.
- Expose budgeted cross-set neighbour queries with independent local indices.
- Describe native NumPy outputs with concrete scalar types for strict consumers.

## 0.5.0 — 2026-10-04

### Research validation and Python parity

- Python geometry now exposes Rust's internal-coordinate forest, Cartesian
  rebuilding and reusable bond-angle-torsion frames, retaining native atom
  indices, missing coordinates and radian units. In-process differential tests
  compare Python and Rust on ubiquitin, including TM-score, GDT and distances.
- Release wheels run the Python contract suite after installation on all five
  supported platform builds. Python formatting, linting and strict typing also
  gate publication.
- Audit runs refuse inconsistent estimands and check policy reads in every
  universe, with matching Python error behavior.
- Fixed-column V2000 parsing handles three-digit atom indices correctly;
  acquisition and study runners verify frozen input SHA-256 before analysis.
- Added reproducible, external-data studies against ProLIF, FreeSASA, Boba,
  Gemmi and EvoEF and a grouped SKEMPI2 alanine-scanning panel. Measured
  sensitivity and descriptive experimental associations do not establish
  predictive superiority or publication readiness. External data and tools
  remain separate from the library distribution.


### Changed

- **Breaking.** An indeterminate analysis carries no value. `Analysis<T>` holds an
  `Outcome<T>`, `Determinate(value)` or `Indeterminate(reason)`; the `Deref` and the
  public `value` and `status` fields are gone, replaced by `value()`, `outcome()`,
  `into_result()` and `status()`. `Quality` qualifies an answer that exists. In
  Python, `Analysis.value` raises `IndeterminateError` when there is none, with
  `is_determinate` and `indeterminacy` to ask first.
- **Breaking.** `ImpactEstimate` (none, low, moderate, high) is replaced by `Impact`:
  unmeasured, inert, or measured by an audit. The assigned labels claimed evidence
  nothing supplied.
- **Breaking.** Re-execution verifies input bytes with SHA-256 (`ContentDigest`,
  `Provenance::input_digest`, serialised as `input_sha256`) instead of the 64-bit
  `Fingerprint`, which stays as the fast identity of policies and cache entries.

### Fixed

- Sampling a density map at a Cartesian position that is exactly a grid point could report no
  value, because the cell transform put the coordinate a few ulps below the first grid point.

- `group_coordinate_variance` centred a group's samples on the group's grand mean, so a
  group of atoms at different places reported its extent in space as fluctuation; it is
  now the mean of each atom's variance about its own mean position.

### Added (policy engine and audit)

- The governed executor applies `assembly` (a biological assembly, or crystal
  contacts within a radius as placed whole chains) and `symmetry`, resolves
  alternate conformations and hydrogens over that system, drives copies from
  per-frame coordinates, and applies `hydrogens = exclude`. Contradictory or
  impossible systems are classified errors; `include_inferred` is refused.
  Descriptors declare the policy fields they read and provenance records them.
- `contacts_by_definition` (Rust `definition_contacts_kernel`): contacts under the
  policy's `contact_def` and `vdw_radii`.
- `molframe-audit` measures any outcome (`OutcomeMetric`: set overlap, scalar error,
  vector difference, Kendall tau, graph difference, conclusion flip), splits the
  variation into effects and pairwise interactions, attributes it by Shapley shares
  and by class of uncertainty, supports forbidden pairs and constrained plans,
  refuses decisions the analysis never read, and reports universes with no answer
  as a fraction. Python: `molframe.audit`. `AnalysisPolicy` names all fourteen
  decisions.
- Python `molframe.crystal`: space groups by Hermann-Mauguin name, Hall symbol or type with
  their symmetry operations, the space group a structure carries, density maps (MRC read
  and write, statistics, histograms, interpolation) and reflection tables (MTZ read and
  write). Symmetry operations print as CIF triplets that parse back to the same operation,
  and map, MRC and reflection errors gain diagnostic codes.
- Python `molframe.compare`: combinatorial-extension alignment (`ce_align`), contact areas
  and the CAD score, contact-map overlap, and sequence-based chain assignment and residue
  mapping.
- Python `molframe.trajectory`: pairwise fitted RMSD, Cartesian PCA, the Procrustes mean
  structure, diffusion maps, k-means, agglomerative and DBSCAN clustering, medoids, mean
  squared displacement, group variance, block convergence, harmonic and population
  ensemble similarity and path similarity, each with provenance. Their errors gain
  diagnostic codes, and the method words are read once in Rust.
- Python `molframe.chemistry`: `Molecule` (an MDL record: graph, header and data fields),
  `read_sdf`, `read_mol`, `write_sdf`, `molecule(structure)` and `smarts(structure, pattern)`.
  MOL, MOL2, SMARTS and automorphism errors gain diagnostic codes (`E1301` is new).
- Python `molframe.surface`: surface points with normals, the solvent-excluded surface as
  a `Mesh` (components, curvatures, edge distances, patches, OBJ), buried surface between
  two molecules, and atom depths. Rust: `IndexedSurfaceMesh::area`, `MoleculeRole::from_str`,
  and diagnostic codes for the mesh, burial and depth errors.
- `analysis.sasa` (Rust `sasa_kernel`): per-atom solvent-accessible surface area under
  the policy's `vdw_radii`, `hydrogens` and `assembly`, so each is a decision an audit
  can vary.
- `ReadOptions::digest_input` (Python `ReadOptions(digest_input=True)`) records the
  SHA-256 of the bytes parsed in the structure's metadata (`input_sha256`,
  `input_name`), and governed analyses carry it in their provenance, so an audit
  certificate can tie its results to the bytes. Off by default: hashing is work on
  every read.
- A descriptor can forbid a resolution, require information and state its estimand.
  `hydrogen_bonds`, `water_bridges` and base pairing refuse `hydrogens = exclude`
  (`E6103`) and are indeterminate, with the reason, on an input that models no
  hydrogens, where they used to return an empty table. `Analysis.estimand` (Rust
  `Provenance::estimand`) says what each analysis estimates.
- Python physical analyses: `leaflets`, `radial_distribution`,
  `coordination_numbers`, `linear_density`, `density_map`, `pore_profile`,
  `surface_contacts`.

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
  (QCP) superposition with a Jacobi fallback. It is performance-neutral
  against the Jacobi fit (the O(n) covariance pass dominates; 38.0 µs against
  37.9 µs for 10,000 atoms), kept for its rotation-free deviation and fallback.
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
- Accessible-area, cell-list pair, bond-perception and contact-potential
  blocks stream into their results, so the peak is the result plus a bounded
  window rather than the result twice; each is tested at 1, 2, 4 and 8 workers
  and under a budget that refuses what cannot fit.
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
