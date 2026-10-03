//! Ordered PDB metadata records.

mod records;
mod remark;
mod typed;

pub use records::*;
pub use remark::{Biomolecule, Biomt, BiomtGroup, MissingResidue};
pub use typed::{Link, ResidueId, SeqRes, SsBond};
