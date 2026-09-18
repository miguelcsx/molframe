//! Structure-scoped operations layered on the domain crates.
//!
//! Each of these takes a [`Structure`](crate::Structure) and projects it into
//! chemistry, geometry or query kernels. They sit apart from `crate::operations`
//! because they are not declarative requests executed by a plan: each is a
//! direct computation, gated on the feature combination it needs and nothing
//! more.

#[cfg(all(feature = "geom", feature = "chem"))]
mod backbone;
#[cfg(all(feature = "chem", feature = "spatial"))]
mod bonds;
#[cfg(feature = "query")]
mod query;
#[cfg(all(feature = "geom", feature = "chem"))]
mod roles;
#[cfg(all(feature = "geom", feature = "chem"))]
mod side_chains;

#[cfg(all(feature = "geom", feature = "chem"))]
pub use backbone::{
    BackboneTorsionRecord, ProteinAlphaTrace, structure_backbone_torsions,
    structure_backbone_torsions_model, structure_protein_alpha_traces,
};
#[cfg(all(feature = "chem", feature = "spatial"))]
pub use bonds::{
    BondInference, BondInferenceReport, DEFAULT_BOND_RADIUS_SCALE, DEFAULT_MINIMUM_BOND_DISTANCE,
    infer_bonds,
};
#[cfg(feature = "query")]
pub use query::QueryStructure;
#[cfg(all(feature = "geom", feature = "chem"))]
pub use side_chains::{
    SideChainTorsionRecord, SideChainTorsionReport, structure_side_chain_torsions,
};
