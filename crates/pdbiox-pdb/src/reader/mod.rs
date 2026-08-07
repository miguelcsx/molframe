//! Reading the legacy fixed-column format.
//!
//! Records are read by column position, because that is what the format is: a
//! value's meaning comes from where it sits, not from any delimiter. A line that
//! stops early is padded, a field that will not parse becomes a finding rather
//! than a guess, and the read continues so the caller gets the structure *and*
//! the list of what was wrong with it.
//!
//! A new residue begins whenever the chain, the residue number or the insertion
//! code changes — the insertion code included, because 163, 163A and 163B are
//! three residues and an antibody numbering scheme depends on it.

mod ensemble;
mod entry;
mod lines;
mod state;

pub use entry::{PdbReader, read, read_pdbqt, read_pqr};
