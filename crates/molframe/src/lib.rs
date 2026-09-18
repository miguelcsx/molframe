//! A batteries-included structural bioinformatics engine.
//!
//! This crate is the composition layer. It owns no kernel of its own: what lives
//! here is the facade that dispatches reading and writing to whichever format
//! crates are linked, the typed operation vocabulary over the analysis kernels
//! ([`Plan`] and its requests), the policy configuration, and the prelude that
//! puts all of it in scope. Everything else is re-exported from the crate that
//! has the logic. The default feature set is the complete user surface; a caller
//! that needs one format can disable default features and link just that format.
//!
//! # Reading, and saying what was wrong with the file
//!
//! ```
//! use molframe::prelude::*;
//!
//! # const PDB: &str = "\
//! # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
//! # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
//! # ATOM      3  C   ALA A   1      13.100   7.520  -6.504  1.00  0.00           C
//! # END
//! # ";
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let (structure, findings) =
//!         read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
//!
//!     println!("{} atoms in {} chains", structure.atom_count(), structure.chain_count());
//!     for finding in &findings {
//!         eprintln!("{}", Rendered::new(finding));
//!     }
//!     Ok(())
//! }
//! ```
//!
//! [`read_with_diagnostics`] is the honest form and is what a pipeline should
//! use: the findings say what was wrong with the file, and a file that parses is
//! not the same as a file that is right. [`read`] discards them for convenience.
//! Both return [`Findings`], which is an [`Error`](std::error::Error) — hence the
//! `?` above, and hence `anyhow::Result` or `Box<dyn Error>` at a `main` without
//! a match anywhere.
//!
//! # Composing an operation
//!
//! The same structure drives the typed operation vocabulary, so a pipeline is
//! data rather than a chain of calls:
//!
//! ```
//! use molframe::prelude::*;
//!
//! # const PDB: &str = "\
//! # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
//! # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
//! # END
//! # ";
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let (structure, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
//!
//!     let mut plan = Plan::new();
//!     plan.add_selection(
//!         "alpha_carbons",
//!         SelectionRequest::new("name CA", AnalysisPolicy::default())?,
//!     )?;
//!
//!     let result = plan.execute(
//!         PlanInput { structure: Some(&structure), ..PlanInput::default() },
//!         &ExecutionContext::default(),
//!     )?;
//!     assert_eq!(result.entries.len(), 1);
//!     Ok(())
//! }
//! ```
//!
//! # Examples are owed, not owed everywhere
//!
//! The in-memory paths above run as doctests. The remaining public items follow
//! the reference documentation of the crate that owns them rather than repeating
//! it here; `CONTRIBUTING.md`'s rule that every public item carry an example is
//! not yet met for the facade's format writers, which is a known debt.

#![forbid(unsafe_code)]

pub use molframe_core as core;

pub use molframe_core::annotation::{
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn,
    AtomAnnotation, AtomAnnotations, COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION,
    HBOND_ACCEPTOR_ANNOTATION, HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION,
    PLDDT_ANNOTATION, POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION,
    STEREO_CONFIGURATION_ANNOTATION,
};
pub use molframe_core::chunk::{
    AtomChunk, AtomChunkStats, AtomRecord, ChunkBuilder, ElementMask, Extremes, ParentMapping,
    TARGET_CHUNK_ATOMS,
};
pub use molframe_core::column::{BitVec, EncodedColumn, Presence, ValidityMask};
pub use molframe_core::contract::{
    AlgorithmId, AltlocPolicy, Analysis, AnalysisParameters, AnalysisPolicy, AssemblyChoice,
    Assumption, AssumptionSource, Coverage, DictionaryVersion, ImpactEstimate, MissingPolicy,
    ModelChoice, Namespace, ParameterValue, PolicyField, ProfileId, Provenance, SourceRef, Status,
    Tolerance,
};
pub use molframe_core::coords::{Aabb, CoordinateBlock, CoordinateGeneration};
pub use molframe_core::diagnostic::{
    Class, Code, ContextItem, Diagnostic, Diagnostics, Findings, Kind, Rendered, Severity,
    Strictness,
};
pub use molframe_core::element::Element;
pub use molframe_core::execution::ExecutionContext;
pub use molframe_core::index::{
    AtomIndex, BondIndex, ChainIndex, EntityIndex, InstanceId, ModelIndex, ResidueIndex,
};
pub use molframe_core::io::{
    AmbiguousResidueBoundaryPolicy, BatchContinuity, Compression, ContinuityLevel, Format,
    InputBuffer, InputKind, Limits, MissingElementPolicy, OutputOptions, ParseMode, ReadOptions,
    ReadResult, Reader, Select, SelectAll, StructureAtomRecord, StructureBatch,
    StructureBatchBuilder, StructureBatchError, collect_structure,
};
pub use molframe_core::provider::{
    AtomEndpoint, BondChunk, BondChunkProvider, BondChunkRecord, ChunkDescriptor, ChunkId,
    ChunkLayout, DatasetCatalog, DatasetDescriptor, DatasetId, FrameChunk, FrameChunkProvider,
    LocalRow, LogicalRow, PayloadKind, PropertyChunk, PropertyChunkProvider, PropertyKind,
    PropertyValue, ProviderError, StructureChunk, StructureChunkProvider, TARGET_CHUNK_BONDS,
};
pub use molframe_core::selection::AtomSelection;
pub use molframe_core::span::{ByteSpan, Position};
pub use molframe_core::structure::{
    AtomRef, ChainRef, ChainSequenceExt, CoordinateEditor, CoordinateStore, CountDifference,
    DifferenceError, EntryMetadata, ExtensionStore, MetadataDifference, MissingResidue, ModelRef,
    ReferenceAlignment, ReferenceSequence, ResidueRef, SEQUENCE_REFERENCES_EXTENSION,
    SequenceMapping, SequenceReferences, Structure, StructureData, StructureDifference,
    StructureDifferenceOptions, StructureEditor, StructureView, UnitCell, ValueDifference,
    structure_difference, validate,
};
pub use molframe_core::symbol::{AltId, Interner, SymbolId};
pub use molframe_core::topology::{EntityKind, PolymerKind, Topology};
pub use molframe_core::{
    BondAdjacency, BondOrder, BondProvenance, BondRecord, BondTable, BondTableBuilder,
};

#[cfg(feature = "adapters")]
pub use molframe_adapters as adapters;

#[cfg(feature = "audit")]
pub use molframe_audit as audit;
#[cfg(feature = "audit")]
pub use molframe_audit::{
    AuditPlan, AuditReport, AuditRun, DimensionSensitivity, PlanError, PolicyDimension,
    PolicySpace, PolicyValue, SensitiveItem, audit,
};

#[cfg(feature = "fx")]
pub use molframe_fx as fx;

#[cfg(feature = "ml")]
pub use molframe_ml as ml;
#[cfg(feature = "ml")]
pub use molframe_ml::{
    ArrowStream, AtomTable as AtomArrowTable, BondTable as BondArrowTable,
    ChainTable as ChainArrowTable, DLDataType, DLDevice, DLManagedTensor, DLTensor, Dataset,
    DatasetError, DatasetFilter, DatasetSplit, DatasetWarning, DlpackError, DlpackTensor,
    EdgeDirection, EdgeFeature, EdgeKind, ExportCost, Graph, GraphError, GraphOptions, LoadError,
    ManifestEntry, MissingFeaturePolicy, MolframeExtension, NodeFeature, NodeLevel,
    ResidueTable as ResidueArrowTable, SplitOptions, SplitRatios, SplitStrategy, TableFileError,
    extension_name, graph, write_atom_ipc, write_atom_ipc_with_metadata, write_atom_parquet,
    write_atom_parquet_with_metadata,
};

#[cfg(feature = "chem")]
pub use molframe_chem as chem;
#[cfg(feature = "chem")]
pub use molframe_chem::{
    AutomorphismLimit, ChemistryReport, CifProvider, Component, ComponentAtom, ComponentBond,
    ComponentCoverage, ComponentKind, ComponentProvider, ElementProperties, EquivalenceCache,
    EquivalenceClasses, IonicRadius, IonicSpin, MemoryProvider, PeoeAtom, PeoeAtomType, PeoeBond,
    PeoeError, PeoeOptions, PeoeParameterProfile, PolymerAtomRole, PolymerLinkPolicy,
    PolymerLinkRule, PolymerRoleProfile, PolymerRoleReport, PolymerRoleRule, RadiusSet,
    RadiusTable, SideChainDefinition, SideChainRoles, SmartsDataError, SmartsError, SmartsMatch,
    SmartsPattern, StereoConfiguration, apply_component_chemistry, apply_polymer_role_profile,
    automorphisms, component_coverage, component_peoe_charges, element_properties,
    equivalence_classes, ionic_radii, peoe_charges, side_chain_definition, vdw_radius,
};

#[cfg(feature = "mmcif")]
pub use molframe_cif as cif;
#[cfg(feature = "mmcif")]
pub use molframe_cif::{
    Category, CifValue, CifWriteError, CifWriteOptions, CifWriteToError, Column, DataBlock,
    Document, write_preserving, write_preserving_to,
};

#[cfg(feature = "modelcif")]
pub use molframe_modelcif as modelcif;
#[cfg(feature = "modelcif")]
pub use molframe_modelcif::{
    GlobalMetric, LocalMetric, MODEL_CIF_EXTENSION, MetricDefinition, ModelCategory, ModelCif,
    ModelCifExt, ModelDescription, ModelRow, PairwiseMetric, ProtocolStep, QualityMetrics,
    SoftwareGroup, Target, Template,
};

#[cfg(feature = "bcif")]
pub use molframe_bcif as bcif;
#[cfg(feature = "bcif")]
pub use molframe_bcif::{BcifReader, BinaryDocument};

#[cfg(feature = "geom")]
pub use molframe_geom as geom;
#[cfg(feature = "geom")]
pub use molframe_geom::{
    BackboneFrame, BackboneResidue, BackboneTorsions, BatchGeometryError, CircularSummary,
    Decomposition, DistanceMatrix, EigenError, EigenOptions, FluctuationError, HelixGeometry,
    MatrixError, PeriodicAngle, PeriodicError, Plane, Rigid, Rotation3, RotationError,
    RotationMeanOptions, RotationOptions, SuperposeError, SuperposeOptions, Superposition,
    TorusMetric, angle, angles_into, asphericity, asphericity_with_options, backbone_frames,
    backbone_torsions, best_fit_plane, best_fit_plane_with_options, centre_of_mass, centroid,
    circular_summary, cross, degrees, dihedral, displacement, distance, distance_matrix,
    distance_matrix_between, distance_squared, distances_into, dot, gyration_axes,
    gyration_axes_with_options, helix_geometry, helix_geometry_with_options, inertia_tensor, norm,
    normalise, path_torsions, plane_deviation, plane_deviation_with_options, principal_axes,
    principal_axes_with_options, radius_of_gyration, rmsd, rmsd_flat, rmsf, rotation_mean,
    rotation_mean_with_options, superpose, superpose_with_options, torsions_into, torus_summary,
};

#[cfg(feature = "ic")]
pub use molframe_ic as ic;
#[cfg(feature = "ic")]
pub use molframe_ic::{
    BatFrame, Dihedron, Hedron, InternalAtom, InternalCoordinates, internal_coordinates, place_atom,
};

#[cfg(feature = "pdb")]
pub use molframe_pdb as pdb;
#[cfg(feature = "pdb")]
pub use molframe_pdb::{
    PDB_HEADERS_EXTENSION, PdbHeaderRecord, PdbHeaders, PdbHeadersExt, PdbIdentifierNamespace,
    PdbOptions, write_mmtf, write_mmtf_to, write_pdbqt, write_pdbqt_to, write_pqr, write_pqr_to,
};

#[cfg(feature = "spatial")]
pub use molframe_spatial as spatial;
#[cfg(feature = "spatial")]
pub use molframe_spatial::{
    AutoBackendProfile, CellGridOptions, CellList, KdPeriodicOptions, KdTree, NeighborList,
    NeighborListOptions, NeighborPair, NeighborSkinProfile, PeriodicBox, PeriodicImage,
    SpatialBackend, SpatialError, SpatialOption, SpatialPlan, SpatialSearchOptions, pairs_within,
    pairs_within_with_options, within, within_with_options,
};

#[cfg(feature = "query")]
pub use molframe_query as query;
#[cfg(feature = "query")]
pub use molframe_query::{Builder as QueryBuilder, Evaluation, Groups, Query, col};

#[cfg(feature = "xtal")]
pub use molframe_xtal as xtal;
#[cfg(feature = "xtal")]
pub use molframe_xtal::{
    ASSEMBLIES_EXTENSION, AffineTransform, AssemblyDef, AssemblyExt, AssemblyNeighbor, AssemblySet,
    AssemblyView, AtomInstance, CellTransform, ChainInstance, CrystalImage, CrystalImageBatch,
    CrystalImageBatchOptions, CrystalNeighbor, CrystalNeighborBatch, CrystalNeighborOptions,
    DEFAULT_CRYSTAL_IMAGE_LIMIT, DEFAULT_INSTANCE_LIMIT, Generator, INSTANCE_ID_ANNOTATION,
    NCS_EXTENSION, NcsAtomInstance, NcsCode, NcsExt, NcsOperator, NcsSet, NcsView, OperExpression,
    Operator, Rational, SYMMETRY_EXTENSION, SpaceGroupSetting, SymmetryExt, SymmetryOperation,
    SymmetrySet, collect_crystal_neighbors, crystal_image_batches, crystal_neighbor_batches,
    lower_assemblies, lower_ncs, lower_symmetry, space_group_by_hall, space_group_setting,
    space_group_settings, visit_crystal_images, visit_crystal_neighbors,
};

// The analysis crates carry many small, related items, so they are
// re-exported under their own namespace rather than flattened into the root.
#[cfg(feature = "analysis")]
pub use molframe_analysis as analysis;
#[cfg(feature = "analysis")]
pub use molframe_analysis::{
    StreamlineDirection, StreamlineOptions, VectorFieldError, VectorFieldGrid,
    integrate_streamlines,
};
#[cfg(feature = "compare")]
pub use molframe_compare as compare;
#[cfg(feature = "seq")]
pub use molframe_seq as seq;
#[cfg(feature = "surface")]
pub use molframe_surface as surface;
#[cfg(feature = "surface")]
pub use molframe_surface::{
    SurfaceComponent, SurfaceComponentError, SurfaceComponentFilter, filter_surface_components,
    surface_components,
};
#[cfg(feature = "traj")]
pub use molframe_traj as traj;
#[cfg(feature = "traj")]
pub use molframe_traj::{
    TrajectoryInterpolation, TrajectoryInterpolationError, interpolate_trajectory_frames,
};
#[cfg(feature = "validate")]
pub use molframe_validate as validate;

mod facade;
#[cfg(all(feature = "analysis", feature = "geom"))]
mod operations;
mod policy_config;
pub mod prelude;
mod structure;

pub use facade::{
    read, read_bytes, read_with_diagnostics, read_with_options, write, write_with_options,
};
// The bounded batch reader needs a format crate to read with, so these two
// names exist under exactly the features that give `StructureBatchReader` a
// variant.
#[cfg(any(
    feature = "mmcif",
    feature = "pdb",
    feature = "bcif",
    feature = "modelcif"
))]
pub use facade::{StructureBatchReader, open_structure_batches};
pub use policy_config::{
    ApplicationConfiguration, ChemistryConfiguration, OutputConfiguration, PolicyConfigError,
    PolicyOverrides, read_configuration, read_policy, read_policy_overrides,
};

#[cfg(all(feature = "analysis", feature = "geom"))]
pub use operations::FloatInput;
// The comparison requests live in `operations`, which needs the two kernels
// its executor is built from, so selecting `compare` alone is not enough for
// the module to exist.
#[cfg(all(feature = "compare", feature = "analysis", feature = "geom"))]
pub use operations::{ComparisonMetric, ComparisonRequest, ComparisonResult};
#[cfg(all(feature = "analysis", feature = "geom"))]
pub use operations::{
    ContactsRequest, CoordinateInput, ExecutionPlanError, FrameInput, GeometryRequest,
    GeometryValue, IndexInput, PhysicalRequest, PhysicalValue, Plan, PlanInput, PlanOperation,
    PlanResult, PlanResultEntry, PlanValue, RmsdRequest, ScalarInput, SelectionRequest,
    SpatialRequest, SpatialValue, StructureRequest, StructureValue,
};
#[cfg(all(feature = "analysis", feature = "geom", feature = "surface"))]
pub use operations::{MaskInput, SurfaceRequest, SurfaceValue};
#[cfg(all(feature = "analysis", feature = "geom", feature = "traj"))]
pub use operations::{TrajectoryRequest, TrajectoryValue};

#[cfg(feature = "geom")]
pub use facade::transform;

#[cfg(all(feature = "geom", feature = "chem"))]
pub use structure::{
    BackboneTorsionRecord, ProteinAlphaTrace, SideChainTorsionRecord, SideChainTorsionReport,
    structure_backbone_torsions, structure_backbone_torsions_model, structure_protein_alpha_traces,
    structure_side_chain_torsions,
};

#[cfg(feature = "mmcif")]
pub use facade::{
    read_document, write_mmcif, write_mmcif_to, write_mmcif_to_with_options,
    write_mmcif_with_options,
};

#[cfg(feature = "bcif")]
pub use facade::{write_bcif, write_bcif_with_options};

#[cfg(feature = "chem")]
pub use facade::read_component_dictionary;

#[cfg(feature = "pdb")]
pub use facade::write_pdb;

#[cfg(feature = "query")]
pub use structure::QueryStructure;

#[cfg(all(feature = "chem", feature = "spatial"))]
pub use structure::{
    BondInference, BondInferenceReport, DEFAULT_BOND_RADIUS_SCALE, DEFAULT_MINIMUM_BOND_DISTANCE,
    infer_bonds,
};
