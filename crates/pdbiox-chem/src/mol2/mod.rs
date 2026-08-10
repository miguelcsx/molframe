//! Tripos MOL2 records.

mod model;
mod parse;
mod write;

pub use model::{Mol2AtomMetadata, Mol2BondMetadata, Mol2Error, Mol2Record, Mol2Section};
pub use parse::parse_mol2_record;
pub use write::write_mol2;

#[cfg(test)]
#[path = "mol2_tests.rs"]
mod tests;
