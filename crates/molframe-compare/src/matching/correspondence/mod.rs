//! Sequence extraction, globally optimal chain assignment and query mapping.

mod assignment;
mod atom;
mod chain;
mod identifier;
mod residue;
mod sequence;

pub use assignment::{ChainAlternative, ChainAssignment, ChainMapping, assign_chains, map_chains};
pub use atom::{AtomCorrespondence, AtomMappingAlternative, map_atoms};
pub use chain::{ChainSequence, chain_sequences, pair_chains};
pub use identifier::pair_residues_by_identifier;
pub use residue::{ResiduePair, map_residues};
pub use sequence::{ResidueMatch, map_sequence_to_structure};

#[cfg(test)]
#[path = "../correspondence_tests.rs"]
mod tests;
