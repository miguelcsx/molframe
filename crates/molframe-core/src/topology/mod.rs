//! The hierarchy, stored as offset tables rather than nested collections.
//!
//! Models contain chains, chains contain residues, residues contain atoms. Each
//! table records where its first child sits and how many children follow, and
//! because children are contiguous, a parent's children are a range. Finding the
//! atoms of a residue is then two loads and an addition — no pointer chasing, no
//! search, and no allocation per residue.

mod chain;
mod csr;
mod entity;
mod hierarchy;
mod model;
mod residue;

pub use chain::{ChainRecord, ChainTable, PolymerKind};
pub use csr::{Csr, CsrBuilder};
pub use entity::{EntityKind, EntityTable};
pub use hierarchy::Topology;
pub use model::ModelTable;
pub use residue::{ResidueRecord, ResidueTable};
