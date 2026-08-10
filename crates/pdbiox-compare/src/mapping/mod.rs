//! Sequence extraction, globally optimal chain assignment and query mapping.

mod assignment;
mod chain;
mod sequence;

pub use assignment::{ChainAlternative, ChainAssignment, ChainMapping, assign_chains, map_chains};
pub use chain::{ChainSequence, chain_sequences};
pub use sequence::{ResidueMatch, map_sequence_to_structure};

#[cfg(test)]
#[path = "../mapping_tests.rs"]
mod tests;
