//! MDL MOL and multi-record SDF input and output.

mod model;
mod parse;
mod write;

pub use model::{
    MolAtom, MolAtomMetadata, MolBond, MolBondMetadata, MolError, MolRecord, MolVersion, Molecule,
    SdfProperty,
};
pub use parse::{parse_mol_record, parse_sdf_records};
pub use write::{write_mol, write_sdf};

#[cfg(test)]
#[path = "mol_tests.rs"]
mod tests;
