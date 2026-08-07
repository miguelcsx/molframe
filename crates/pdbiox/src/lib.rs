//! A batteries-included structural bioinformatics engine.
//!
//! This crate is a facade. It contains no logic of its own — only re-exports,
//! the prelude, and the feature flags that decide which formats get linked. A
//! caller that reads one format does not pay for the others.
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

pub use pdbiox_core::annotation::{
    ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn, AtomAnnotation,
    AtomAnnotations, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION, PLDDT_ANNOTATION,
    SEGMENT_ID_ANNOTATION,
};
pub use pdbiox_core::chunk::{AtomChunk, AtomRecord, ChunkBuilder};
pub use pdbiox_core::column::{BitVec, EncodedColumn, Presence, ValidityMask};
pub use pdbiox_core::contract::DictionaryVersion;
pub use pdbiox_core::contract::{
    AltlocPolicy, Analysis, AnalysisPolicy, AssemblyChoice, Assumption, Coverage, MissingPolicy,
    ModelChoice, Namespace, PolicyField, ProfileId, Provenance, Status,
};
pub use pdbiox_core::coords::{Aabb, CoordinateBlock, CoordinateGeneration};
pub use pdbiox_core::diagnostic::{
    Class, Code, ContextItem, Diagnostic, Diagnostics, Rendered, Severity, Strictness,
};
pub use pdbiox_core::element::Element;
pub use pdbiox_core::index::{
    AtomIndex, BondIndex, ChainIndex, EntityIndex, InstanceId, ModelIndex, ResidueIndex,
};
pub use pdbiox_core::io::{
    Format, InputBuffer, InputKind, Limits, ParseMode, ReadOptions, ReadResult, Reader, Select,
    SelectAll,
};
pub use pdbiox_core::selection::AtomSelection;
pub use pdbiox_core::span::{ByteSpan, Position};
pub use pdbiox_core::structure::{
    AtomRef, ChainRef, ChainSequenceExt, CoordinateEditor, CoordinateStore, EntryMetadata,
    ExtensionStore, MissingResidue, ModelRef, ReferenceAlignment, ReferenceSequence, ResidueRef,
    SEQUENCE_REFERENCES_EXTENSION, SequenceMapping, SequenceReferences, Structure, StructureData,
    StructureEditor, StructureView, UnitCell, validate,
};
pub use pdbiox_core::symbol::{AltId, Interner, SymbolId};
pub use pdbiox_core::topology::{EntityKind, PolymerKind, Topology};
pub use pdbiox_core::{BondAdjacency, BondOrder, BondProvenance, BondRecord, BondTable};

#[cfg(feature = "arrow")]
pub use pdbiox_ml::{
    ArrowStream, AtomTable as AtomArrowTable, BondTable as BondArrowTable,
    ChainTable as ChainArrowTable, ExportCost, PdbioxExtension, ResidueTable as ResidueArrowTable,
    extension_name,
};

#[cfg(feature = "chem")]
pub use pdbiox_chem::{
    ChemistryReport, CifProvider, Component, ComponentAtom, ComponentBond, ComponentCoverage,
    ComponentKind, ComponentProvider, ElementProperties, EquivalenceCache, EquivalenceClasses,
    IonicRadius, IonicSpin, MemoryProvider, RadiusSet, RadiusTable, StereoConfiguration,
    apply_component_chemistry, component_coverage, element_properties, equivalence_classes,
    ionic_radii, read_ccd, vdw_radius,
};

#[cfg(feature = "mmcif")]
pub use pdbiox_cif::{CifValue, Document, write_preserving};

#[cfg(feature = "bcif")]
pub use pdbiox_bcif::{BcifReader, BinaryDocument};

#[cfg(feature = "geom")]
pub use pdbiox_geom::{
    BackboneResidue, BackboneTorsions, DistanceMatrix, Rigid, SuperposeError, Superposition, angle,
    backbone_torsions, centroid, dihedral, distance, distance_matrix, distance_matrix_between,
    radius_of_gyration, rmsd, superpose,
};

#[cfg(feature = "ic")]
pub use pdbiox_ic::{
    Dihedron, Hedron, InternalAtom, InternalCoordinates, internal_coordinates, place_atom,
};

#[cfg(feature = "pdb")]
pub use pdbiox_pdb::{PdbOptions, write_pdbqt, write_pqr};

#[cfg(feature = "spatial")]
pub use pdbiox_spatial::{
    CellList, KdTree, NeighborList, NeighborPair, PeriodicBox, SpatialBackend, SpatialError,
    pairs_within, within,
};

#[cfg(feature = "query")]
pub use pdbiox_query::{Builder as QueryBuilder, Evaluation, Groups, Query, col};

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

#[cfg(all(feature = "chem", feature = "spatial"))]
mod chemistry_api;
mod facade;
#[cfg(feature = "geom")]
mod geometry_api;
pub mod prelude;
#[cfg(feature = "query")]
mod query_api;
#[cfg(all(feature = "geom", feature = "chem"))]
mod side_chain_api;

pub use facade::{
    default_limits, read, read_bytes, read_with_diagnostics, read_with_options, write,
};

#[cfg(feature = "geom")]
pub use facade::transform;

#[cfg(feature = "geom")]
pub use geometry_api::{BackboneTorsionRecord, structure_backbone_torsions};

#[cfg(all(feature = "geom", feature = "chem"))]
pub use side_chain_api::{
    SideChainTorsionRecord, SideChainTorsionReport, structure_side_chain_torsions,
};

#[cfg(feature = "mmcif")]
pub use facade::{read_document, write_mmcif};

#[cfg(feature = "bcif")]
pub use facade::write_bcif;

#[cfg(feature = "chem")]
pub use facade::read_component_dictionary;

#[cfg(feature = "pdb")]
pub use facade::write_pdb;

#[cfg(feature = "query")]
pub use query_api::QueryStructure;

#[cfg(all(feature = "chem", feature = "spatial"))]
pub use chemistry_api::{BondInference, BondInferenceReport, infer_bonds};
