//! Residue-aware lDDT between two structures.
//!
//! The coordinate kernel in [`crate::distance_difference`] scores two
//! equally long point lists and knows nothing of residues, so it cannot say that
//! two oxygens of a carboxylate are interchangeable, that a residue missing from
//! the model is a failure, or that atoms of one residue do not test each other.
//! This module supplies those definitions on top of structures.

mod options;
mod resolve;
mod score;
mod symmetry;
mod table;

pub use options::{ResidueLddt, ResidueLddtOptions, ResidueScore, SymmetryPolicy};
pub use score::residue_lddt;
pub use table::ResidueCorrespondence;
