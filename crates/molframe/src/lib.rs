//! A batteries-included structural bioinformatics engine.
//!
//! This crate is the composition layer, and it is the whole 1.0 contract: what
//! is reachable from the root or from [`prelude`] is curated and stable, and
//! everything else — the raw storage engine — lives under [`engine`],
//! documented there as the uncurated, lower-stability door. The default feature
//! set includes mmCIF and PDB; callers opt into analysis domains independently.
//!
//! # What the root carries
//!
//! The root namespace is curated, not a mirror of the workspace. It holds the
//! curated structure handles ([`Structure`], [`Selection`], [`StructureEditor`],
//! the [`Chains`]/[`Residues`]/[`Atoms`]/[`Models`] collections), the reading and
//! writing verbs, and the policy and diagnostic vocabulary.
//!
//! Everything else is module-only. `molframe::analysis::hydrogen_bonds`, never a
//! flat `molframe::hydrogen_bonds`; `molframe::chemistry::vdw_radius`, never
//! `molframe::vdw_radius`; `molframe::interop::AtomTable`, never
//! `molframe::AtomTable`. `geometry`, `analysis`, `compare`, `sequence`,
//! `surface`, `validation`, `trajectory`, `crystal`, `motif`, `chemistry`,
//! `spatial`, `query`, `ic`, `interop`, `audit` and `adapters` are all reached
//! that way, so the root stays a page a reader can hold. The four structure
//! formats live under [`formats`].
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
//! # Navigating and editing
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
//!     for chain in structure.chains() {
//!         println!("chain {:?}: {} residues", chain.label(), chain.residues().count());
//!     }
//!
//!     let mut editor = structure.edit();
//!     editor
//!         .rename_chain(molframe::ChainIndex::new(0), "B")
//!         .map_err(Findings::from)?;
//!     let renamed = editor.finish().map_err(Findings::from)?;
//!     assert_eq!(renamed.chain_at(0).and_then(|chain| chain.label()), Some("B"));
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

pub use molframe_core::annotation::{
    AROMATIC_ATOM_ANNOTATION, ATOM_RADIUS_ANNOTATION, AUTODOCK_TYPE_ANNOTATION, AnnotationColumn,
    AtomAnnotation, AtomAnnotations, COMPONENT_KIND_ANNOTATION, FORMAL_CHARGE_ANNOTATION,
    HBOND_ACCEPTOR_ANNOTATION, HBOND_DONOR_ANNOTATION, PAE_ANNOTATION, PARTIAL_CHARGE_ANNOTATION,
    PLDDT_ANNOTATION, POLYMER_ATOM_ROLE_ANNOTATION, SEGMENT_ID_ANNOTATION,
    STEREO_CONFIGURATION_ANNOTATION,
};
pub use molframe_core::contract::{
    AlgorithmId, AltlocPolicy, Analysis, AnalysisParameters, AnalysisPolicy, AssemblyChoice,
    Assumption, AssumptionSource, Coverage, DictionaryVersion, ImpactEstimate, MissingPolicy,
    ModelChoice, Namespace, ParameterValue, PolicyField, ProfileId, Provenance, SourceRef, Status,
    Tolerance,
};
pub use molframe_core::coords::Aabb;
pub use molframe_core::diagnostic::{
    Class, Code, ContextItem, Diagnostic, Diagnostics, Findings, Kind, Rendered, Severity,
    Strictness,
};
pub use molframe_core::element::Element;
pub use molframe_core::execution::{ExecutionContext, MemoryBudgetError};
pub use molframe_core::index::{
    AtomIndex, BondIndex, ChainIndex, EntityIndex, InstanceId, ModelIndex, ResidueIndex,
};
pub use molframe_core::io::{
    AmbiguousResidueBoundaryPolicy, BatchContinuity, Compression, ContinuityLevel, Format,
    InputBuffer, InputKind, Limits, MissingElementPolicy, OutputOptions, ParseMode, ReadOptions,
    ReadResult, Reader, Select, SelectAll, StructureAtomRecord, StructureBatch,
    StructureBatchBuilder, StructureBatchError, collect_structure,
};
pub use molframe_core::span::{ByteSpan, Position};
pub use molframe_core::structure::{
    AtomRef, ChainRef, ChainSequenceExt, CoordinateEditor, CountDifference, DifferenceError,
    EntryMetadata, ExtensionStore, MetadataDifference, MissingResidue, ModelRef,
    ReferenceAlignment, ReferenceSequence, ResidueRef, SEQUENCE_REFERENCES_EXTENSION,
    SequenceMapping, SequenceReferences, StructureDifference, StructureDifferenceOptions, UnitCell,
    ValueDifference,
};
pub use molframe_core::symbol::{AltId, SymbolId};
pub use molframe_core::topology::{EntityKind, PolymerKind, Topology};
pub use molframe_core::{
    BondAdjacency, BondOrder, BondProvenance, BondRecord, BondTable, BondTableBuilder,
};
pub use molframe_engine::{
    CompiledWorkflow, Cost, Explanation, Input, Node, OperationMetadata, Output, PhysicalNode,
    Workflow, WorkflowBuildError, WorkflowError, WorkflowInputs, WorkflowResults,
};

#[cfg(feature = "adapters")]
pub use molframe_adapters as adapters;

#[cfg(feature = "audit")]
pub use molframe_audit as audit;

#[cfg(feature = "motif")]
pub use molframe_fx as motif;

#[cfg(feature = "interop")]
pub use molframe_interop as interop;

#[cfg(feature = "chemistry")]
pub use molframe_chem as chemistry;

pub mod formats;

#[cfg(feature = "geometry")]
pub use molframe_geom as geometry;

#[cfg(feature = "ic")]
pub use molframe_ic as ic;

#[cfg(feature = "spatial")]
pub use molframe_spatial as spatial;

#[cfg(feature = "query")]
pub use molframe_query as query;
#[cfg(feature = "query")]
pub use molframe_query::{Query, QueryFingerprint};

#[cfg(feature = "crystal")]
pub use molframe_xtal as crystal;

// The analysis crates carry many small, related items, so they are
// re-exported under their own namespace rather than flattened into the root.
#[cfg(feature = "analysis")]
pub use molframe_analysis as analysis;
#[cfg(feature = "compare")]
pub use molframe_compare as compare;
#[cfg(feature = "sequence")]
pub use molframe_seq as sequence;
#[cfg(feature = "surface")]
pub use molframe_surface as surface;
#[cfg(feature = "trajectory")]
pub use molframe_traj as trajectory;
#[cfg(feature = "trajectory")]
pub use molframe_traj::Trajectory;
#[cfg(feature = "validation")]
pub use molframe_validate as validation;

/// The uncurated, low-level door onto the storage engine.
///
/// Everything under [`engine::core`] is a direct re-export of `molframe-core`
/// internals: chunked columnar storage, provider machinery, and the snapshot
/// representation the curated [`Structure`] wraps. Declarative execution is
/// the stable typed [`Workflow`] API at the crate root; there is no second
/// low-level batch executor.
pub mod engine {
    pub mod core;
}

mod extensions;
mod facade;
mod policy_config;
pub mod prelude;
mod structure;

pub use facade::{
    WriteOptions, read, read_buffer, read_bytes, read_with_diagnostics, read_with_options, write,
    write_with_options,
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
    PolicyOverrides, read_configuration, read_policy,
};
pub use structure::{
    Atoms, AtomsIter, Chains, ChainsIter, Models, ModelsIter, Residues, ResiduesIter, Selection,
    Structure, StructureEditor, structure_difference,
};

// The method syntax for the domain kernels whose first argument is a
// structure. See `extensions` for why geometry, surface and the spatial
// planners have no trait here.
#[cfg(feature = "analysis")]
pub use extensions::AnalysisExt;
#[cfg(feature = "compare")]
pub use extensions::CompareExt;
#[cfg(feature = "validation")]
pub use extensions::ValidationExt;

#[cfg(feature = "geometry")]
pub use facade::transform;

#[cfg(feature = "mmcif")]
pub use facade::{read_document, write_mmcif, write_mmcif_with_options};

#[cfg(feature = "bcif")]
pub use facade::{write_bcif, write_bcif_with_options};

#[cfg(feature = "pdb")]
pub use facade::write_pdb;

#[cfg(all(feature = "geometry", feature = "chemistry"))]
pub use structure::{
    BackboneTorsionRecord, ProteinAlphaTrace, SideChainTorsionRecord, SideChainTorsionReport,
    structure_backbone_torsions, structure_backbone_torsions_model, structure_protein_alpha_traces,
    structure_side_chain_torsions,
};

#[cfg(feature = "chemistry")]
pub use facade::read_component_dictionary;

#[cfg(feature = "query")]
pub use structure::QueryStructure;

#[cfg(all(feature = "chemistry", feature = "spatial"))]
pub use structure::{
    BondInference, BondInferenceReport, DEFAULT_BOND_RADIUS_SCALE, DEFAULT_MINIMUM_BOND_DISTANCE,
    infer_bonds,
};
