//! The curated structure handle, plus structure-scoped operations layered on
//! the domain crates.
//!
//! [`Structure`], [`Selection`], [`StructureEditor`] and the collection types
//! are the facade's owned encapsulation of the engine's hierarchy handles —
//! see `handle.rs`, `selection.rs`, `editor.rs` and `collections.rs`. The
//! remaining submodules take a [`Structure`] and project it into chemistry,
//! geometry or query kernels. Each is a direct computation, gated on the
//! feature combination it needs and nothing more.

mod collections;
mod difference;
mod editor;
mod extensions;
mod handle;
mod selection;

#[cfg(all(feature = "geometry", feature = "chemistry"))]
mod backbone;
#[cfg(all(feature = "chemistry", feature = "spatial"))]
mod bonds;
#[cfg(feature = "query")]
mod query;
#[cfg(all(feature = "geometry", feature = "chemistry"))]
mod roles;
#[cfg(all(feature = "geometry", feature = "chemistry"))]
mod side_chains;

pub use collections::{
    Atoms, AtomsIter, Chains, ChainsIter, Models, ModelsIter, Residues, ResiduesIter,
};
pub use difference::structure_difference;
pub use editor::StructureEditor;
pub use handle::Structure;
pub use selection::Selection;

#[cfg(all(feature = "geometry", feature = "chemistry"))]
pub use backbone::{
    BackboneTorsionRecord, ProteinAlphaTrace, structure_backbone_torsions,
    structure_backbone_torsions_model, structure_protein_alpha_traces,
};
#[cfg(all(feature = "chemistry", feature = "spatial"))]
pub use bonds::{
    BondInference, BondInferenceReport, DEFAULT_BOND_RADIUS_SCALE, DEFAULT_MINIMUM_BOND_DISTANCE,
    infer_bonds,
};
#[cfg(feature = "query")]
pub use query::QueryStructure;
#[cfg(all(feature = "geometry", feature = "chemistry"))]
pub use side_chains::{
    SideChainTorsionRecord, SideChainTorsionReport, structure_side_chain_torsions,
};
