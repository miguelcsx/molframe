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

pub use pdbiox_core::chunk::{AtomChunk, AtomRecord, ChunkBuilder};
pub use pdbiox_core::column::{BitVec, EncodedColumn, Presence, ValidityMask};
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
    Format, InputBuffer, Limits, ParseMode, ReadOptions, ReadResult, Reader, Select, SelectAll,
};
pub use pdbiox_core::selection::AtomSelection;
pub use pdbiox_core::span::{ByteSpan, Position};
pub use pdbiox_core::structure::{
    AtomRef, ChainRef, CoordinateStore, EntryMetadata, ModelRef, ResidueRef, Structure,
    StructureData, StructureView, UnitCell, validate,
};
pub use pdbiox_core::symbol::{AltId, Interner, SymbolId};
pub use pdbiox_core::topology::{EntityKind, PolymerKind, Topology};

#[cfg(feature = "mmcif")]
pub use pdbiox_cif::{CifValue, Document, write_preserving};

#[cfg(feature = "geom")]
pub use pdbiox_geom::{
    Rigid, SuperposeError, Superposition, angle, centroid, dihedral, distance, radius_of_gyration,
    rmsd, superpose,
};

#[cfg(feature = "pdb")]
pub use pdbiox_pdb::PdbOptions;

mod facade;

pub use facade::{default_limits, read, read_bytes, read_with_diagnostics, read_with_options};

#[cfg(feature = "mmcif")]
pub use facade::{read_document, write_mmcif};

#[cfg(feature = "pdb")]
pub use facade::write_pdb;

/// The names most programs want in scope.
pub mod prelude {
    pub use crate::{
        Analysis, AnalysisPolicy, AtomRef, AtomSelection, ChainRef, Code, Coverage, Diagnostic,
        Element, ModelIndex, ModelRef, ParseMode, ReadOptions, Rendered, ResidueRef, Status,
        Structure, StructureView,
    };
}

/// Capabilities this build does not include yet, and the crate each will come
/// from. Named here so the gap is visible from inside the library rather than
/// only from its plan.
pub mod roadmap {
    /// The selection language, from `pdbiox-query`.
    pub const SELECTION: &str = "pdbiox-query";
    /// Neighbour search, from `pdbiox-spatial`.
    pub const SPATIAL: &str = "pdbiox-spatial";
}
