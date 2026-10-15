//! GROMACS portable run-input topology reader.

mod constants;
mod error;
mod model;
mod parameters;
mod parser;
mod topology;
mod xdr;

pub use error::TprError;
pub use model::{TprAtom, TprBond, TprHeader, TprResidue, TprTopology};
pub use parser::parse_tpr;

#[cfg(test)]
#[path = "tpr_tests.rs"]
mod tests;
