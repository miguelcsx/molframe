//! A batteries-included structural bioinformatics engine.
//!
//! This crate is a facade. It contains no logic of its own — only re-exports,
//! the prelude, and the feature flags that decide which formats get linked.
//! The default is the complete user surface. A caller that needs only one
//! format can disable default features and link just that format.
//!
//! # Reading a structure
//!
//! ```no_run
//! use pdbiox::prelude::*;
//!
//! # fn main() -> Result<(), Vec<pdbiox::Diagnostic>> {
//! let (structure, findings) = pdbiox::read_with_diagnostics("1abc.pdb")?;
//!
//! println!("{} atoms in {} chains", structure.atom_count(), structure.chain_count());
//! for finding in &findings {
//!     eprintln!("{}", Rendered::new(finding));
//! }
//! # Ok(())
//! # }
//! ```
//!
//! `read_with_diagnostics` is the honest form and is what a pipeline should use:
//! the findings say what was wrong with the file, and a file that parses is not
//! the same as a file that is right. [`read`] discards them for convenience.

#![forbid(unsafe_code)]

pub use pdbiox_core as core;

pub use pdbiox_core::annotation::{
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn,
    AtomAnnotation, AtomAnnotations, COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION,
    HBOND_ACCEPTOR_ANNOTATION, HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION,
    PLDDT_ANNOTATION, POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION,
    STEREO_CONFIGURATION_ANNOTATION,
};
pub use pdbiox_core::chunk::{
    AtomChunk, AtomChunkStats, AtomRecord, ChunkBuilder, ElementMask, Extremes, ParentMapping,
    TARGET_CHUNK_ATOMS,
};
pub use pdbiox_core::column::{BitVec, EncodedColumn, Presence, ValidityMask};
pub use pdbiox_core::contract::{
    AlgorithmId, AltlocPolicy, Analysis, AnalysisParameters, AnalysisPolicy, AssemblyChoice,
    Assumption, AssumptionSource, Coverage, DictionaryVersion, ImpactEstimate, MissingPolicy,
    ModelChoice, Namespace, ParameterValue, PolicyField, ProfileId, Provenance, SourceRef, Status,
    Tolerance,
};
pub use pdbiox_core::coords::{Aabb, CoordinateBlock, CoordinateGeneration};
pub use pdbiox_core::diagnostic::{
    Class, Code, ContextItem, Diagnostic, Diagnostics, Kind, Rendered, Severity, Strictness,
};
pub use pdbiox_core::element::Element;
pub use pdbiox_core::index::{
    AtomIndex, BondIndex, ChainIndex, EntityIndex, InstanceId, ModelIndex, ResidueIndex,
};
pub use pdbiox_core::io::{
    AmbiguousResidueBoundaryPolicy, Compression, Format, InputBuffer, InputKind, Limits,
    MissingElementPolicy, ParseMode, ReadOptions, ReadResult, Reader, Select, SelectAll,
};
pub use pdbiox_core::selection::AtomSelection;
pub use pdbiox_core::span::{ByteSpan, Position};
pub use pdbiox_core::structure::{
    AtomRef, ChainRef, ChainSequenceExt, CoordinateEditor, CoordinateStore, CountDifference,
    DifferenceError, EntryMetadata, ExtensionStore, MetadataDifference, MissingResidue, ModelRef,
    ReferenceAlignment, ReferenceSequence, ResidueRef, SEQUENCE_REFERENCES_EXTENSION,
    SequenceMapping, SequenceReferences, Structure, StructureData, StructureDifference,
    StructureDifferenceOptions, StructureEditor, StructureView, UnitCell, ValueDifference,
    structure_difference, validate,
};
pub use pdbiox_core::symbol::{AltId, Interner, SymbolId};
pub use pdbiox_core::topology::{EntityKind, PolymerKind, Topology};
pub use pdbiox_core::{
    BondAdjacency, BondOrder, BondProvenance, BondRecord, BondTable, BondTableBuilder,
};

#[cfg(feature = "adapters")]
pub use pdbiox_adapters as adapters;

#[cfg(feature = "audit")]
pub use pdbiox_audit as audit;
#[cfg(feature = "audit")]
pub use pdbiox_audit::{
    AuditPlan, AuditReport, AuditRun, DimensionSensitivity, PlanError, PolicyDimension,
    PolicySpace, PolicyValue, SensitiveItem, audit,
};

#[cfg(feature = "fx")]
pub use pdbiox_fx as fx;

#[cfg(feature = "ml")]
pub use pdbiox_ml as ml;
#[cfg(feature = "ml")]
pub use pdbiox_ml::{
    ArrowStream, AtomTable as AtomArrowTable, BondTable as BondArrowTable,
    ChainTable as ChainArrowTable, DLDataType, DLDevice, DLManagedTensor, DLTensor, Dataset,
    DatasetError, DatasetFilter, DatasetSplit, DatasetWarning, DlpackError, DlpackTensor,
    EdgeDirection, EdgeFeature, EdgeKind, ExportCost, Graph, GraphError, GraphOptions, LoadError,
    ManifestEntry, MissingFeaturePolicy, NodeFeature, NodeLevel, PdbioxExtension,
    ResidueTable as ResidueArrowTable, SplitOptions, SplitRatios, SplitStrategy, TableFileError,
    extension_name, graph, write_atom_ipc, write_atom_ipc_with_metadata, write_atom_parquet,
    write_atom_parquet_with_metadata,
};

#[cfg(feature = "chem")]
pub use pdbiox_chem as chem;
#[cfg(feature = "chem")]
pub use pdbiox_chem::{
    AutomorphismLimit, ChemistryReport, CifProvider, Component, ComponentAtom, ComponentBond,
    ComponentCoverage, ComponentKind, ComponentProvider, ElementProperties, EquivalenceCache,
    EquivalenceClasses, IonicRadius, IonicSpin, MemoryProvider, PeoeAtom, PeoeAtomType, PeoeBond,
    PeoeError, PeoeOptions, PeoeParameterProfile, PolymerAtomRole, PolymerLinkPolicy,
    PolymerLinkRule, PolymerRoleProfile, PolymerRoleReport, PolymerRoleRule, RadiusSet,
    RadiusTable, SideChainDefinition, SideChainRoles, SmartsDataError, SmartsError, SmartsMatch,
    SmartsPattern, StereoConfiguration, apply_component_chemistry, apply_polymer_role_profile,
    automorphisms, component_coverage, component_peoe_charges, element_properties,
    equivalence_classes, ionic_radii, peoe_charges, read_ccd, side_chain_definition, vdw_radius,
};

#[cfg(feature = "mmcif")]
pub use pdbiox_cif as cif;
#[cfg(feature = "mmcif")]
pub use pdbiox_cif::{
    Category, CifValue, CifWriteError, CifWriteOptions, Column, DataBlock, Document,
    write_preserving,
};

#[cfg(feature = "modelcif")]
pub use pdbiox_modelcif as modelcif;
#[cfg(feature = "modelcif")]
pub use pdbiox_modelcif::{
    GlobalMetric, LocalMetric, MODEL_CIF_EXTENSION, MetricDefinition, ModelCategory, ModelCif,
    ModelCifExt, ModelDescription, ModelRow, PairwiseMetric, ProtocolStep, QualityMetrics,
    SoftwareGroup, Target, Template,
};

#[cfg(feature = "bcif")]
pub use pdbiox_bcif as bcif;
#[cfg(feature = "bcif")]
pub use pdbiox_bcif::{BcifReader, BinaryDocument};

#[cfg(feature = "geom")]
pub use pdbiox_geom as geom;
#[cfg(feature = "geom")]
pub use pdbiox_geom::{
    BackboneFrame, BackboneResidue, BackboneTorsions, CircularSummary, Decomposition,
    DistanceMatrix, EigenError, EigenOptions, FluctuationError, HelixGeometry, MatrixError,
    PeriodicAngle, PeriodicError, Plane, Rigid, Rotation3, RotationError, RotationMeanOptions,
    RotationOptions, SuperposeError, SuperposeOptions, Superposition, TorusMetric, angle,
    asphericity, asphericity_with_options, backbone_frames, backbone_torsions, best_fit_plane,
    best_fit_plane_with_options, centre_of_mass, centroid, circular_summary, cross, degrees,
    dihedral, displacement, distance, distance_matrix, distance_matrix_between, distance_squared,
    dot, gyration_axes, gyration_axes_with_options, helix_geometry, helix_geometry_with_options,
    inertia_tensor, norm, normalise, path_torsions, plane_deviation, plane_deviation_with_options,
    principal_axes, principal_axes_with_options, radius_of_gyration, rmsd, rmsd_flat, rmsf,
    rotation_mean, rotation_mean_with_options, superpose, superpose_with_options, torus_summary,
};

#[cfg(feature = "ic")]
pub use pdbiox_ic as ic;
#[cfg(feature = "ic")]
pub use pdbiox_ic::{
    BatFrame, Dihedron, Hedron, InternalAtom, InternalCoordinates, internal_coordinates, place_atom,
};

#[cfg(feature = "pdb")]
pub use pdbiox_pdb as pdb;
#[cfg(feature = "pdb")]
pub use pdbiox_pdb::{
    PDB_HEADERS_EXTENSION, PdbHeaderRecord, PdbHeaders, PdbHeadersExt, PdbIdentifierNamespace,
    PdbOptions, read_mmtf, write_mmtf, write_pdbqt, write_pqr,
};

#[cfg(feature = "spatial")]
pub use pdbiox_spatial as spatial;
#[cfg(feature = "spatial")]
pub use pdbiox_spatial::{
    AutoBackendProfile, CellGridOptions, CellList, KdPeriodicOptions, KdTree, NeighborList,
    NeighborListOptions, NeighborPair, NeighborSkinProfile, PeriodicBox, PeriodicImage,
    SpatialBackend, SpatialError, SpatialOption, SpatialPlan, SpatialSearchOptions, pairs_within,
    pairs_within_with_options, within, within_with_options,
};

#[cfg(feature = "query")]
pub use pdbiox_query as query;
#[cfg(feature = "query")]
pub use pdbiox_query::{Builder as QueryBuilder, Evaluation, Groups, Query, col};

#[cfg(feature = "xtal")]
pub use pdbiox_xtal as xtal;
#[cfg(feature = "xtal")]
pub use pdbiox_xtal::{
    ASSEMBLIES_EXTENSION, AffineTransform, AssemblyDef, AssemblyExt, AssemblyNeighbor, AssemblySet,
    AssemblyView, AtomInstance, CellTransform, ChainInstance, CrystalNeighbor,
    DEFAULT_CRYSTAL_IMAGE_LIMIT, DEFAULT_INSTANCE_LIMIT, Generator, INSTANCE_ID_ANNOTATION,
    NCS_EXTENSION, NcsAtomInstance, NcsCode, NcsExt, NcsOperator, NcsSet, NcsView, OperExpression,
    Operator, Rational, SYMMETRY_EXTENSION, SpaceGroupSetting, SymmetryExt, SymmetryOperation,
    SymmetrySet, crystal_neighbors, crystal_neighbors_with_backend, crystal_neighbors_with_limit,
    lower_assemblies, lower_ncs, lower_symmetry, space_group_by_hall, space_group_setting,
    space_group_settings,
};

// The analysis-science crates carry many small, related items, so they are
// re-exported under their own namespace rather than flattened into the root.
#[cfg(feature = "analysis")]
pub use pdbiox_analysis as analysis;
#[cfg(feature = "analysis")]
pub use pdbiox_analysis::{
    StreamlineDirection, StreamlineOptions, VectorFieldError, VectorFieldGrid,
    integrate_streamlines,
};
#[cfg(feature = "compare")]
pub use pdbiox_compare as compare;
#[cfg(feature = "seq")]
pub use pdbiox_seq as seq;
#[cfg(feature = "surface")]
pub use pdbiox_surface as surface;
#[cfg(feature = "surface")]
pub use pdbiox_surface::{
    SurfaceComponent, SurfaceComponentError, SurfaceComponentFilter, filter_surface_components,
    surface_components,
};
#[cfg(feature = "traj")]
pub use pdbiox_traj as traj;
#[cfg(feature = "traj")]
pub use pdbiox_traj::{
    TrajectoryInterpolation, TrajectoryInterpolationError, interpolate_trajectory_frames,
};
#[cfg(feature = "validate")]
pub use pdbiox_validate as validate;

#[cfg(all(feature = "chem", feature = "spatial"))]
mod chemistry_api;
mod facade;
#[cfg(all(feature = "geom", feature = "chem"))]
mod geometry_api;
#[cfg(all(feature = "analysis", feature = "geom"))]
mod operations;
mod policy_config;
#[cfg(all(feature = "geom", feature = "chem"))]
mod polymer_roles;
pub mod prelude;
#[cfg(feature = "query")]
mod query_api;
#[cfg(all(feature = "geom", feature = "chem"))]
mod side_chain_api;

pub use facade::{
    default_limits, read, read_bytes, read_with_diagnostics, read_with_options, write,
};
pub use policy_config::{
    ApplicationConfiguration, ChemistryConfiguration, OutputConfiguration, PolicyConfigError,
    PolicyOverrides, read_configuration, read_policy, read_policy_overrides,
};

#[cfg(all(feature = "analysis", feature = "geom"))]
pub use operations::FloatInput;
#[cfg(feature = "compare")]
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
pub use geometry_api::{
    BackboneTorsionRecord, ProteinAlphaTrace, structure_backbone_torsions,
    structure_backbone_torsions_model, structure_protein_alpha_traces,
};

#[cfg(all(feature = "geom", feature = "chem"))]
pub use side_chain_api::{
    SideChainTorsionRecord, SideChainTorsionReport, structure_side_chain_torsions,
};

#[cfg(feature = "mmcif")]
pub use facade::{read_document, write_mmcif, write_mmcif_with_options};

#[cfg(feature = "bcif")]
pub use facade::{write_bcif, write_bcif_with_options};

#[cfg(feature = "chem")]
pub use facade::read_component_dictionary;

#[cfg(feature = "pdb")]
pub use facade::write_pdb;

#[cfg(feature = "query")]
pub use query_api::QueryStructure;

#[cfg(all(feature = "chem", feature = "spatial"))]
pub use chemistry_api::{
    BondInference, BondInferenceReport, DEFAULT_BOND_RADIUS_SCALE, DEFAULT_MINIMUM_BOND_DISTANCE,
    infer_bonds,
};
