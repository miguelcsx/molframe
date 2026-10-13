//! Small-molecule library chemistry: the properties a screening campaign filters
//! on, a Morgan-style fingerprint, and similarity clustering.
//!
//! Everything here is a graph computation over a [`MolRecord`](crate::mol::MolRecord); no geometry is
//! consulted. Each result says what it did not establish: hydrogens are taken
//! from the file when it has explicit ones and otherwise from standard valences
//! (which is not possible for aromatic-typed bonds without explicit hydrogens,
//! so those molecules get no mass or donor count rather than a guess), the
//! aromatic ring count is the conservative one of `molframe-chem`'s own
//! aromaticity perception, and no octanol-water partition is estimated.

mod descriptors;
mod graph;
mod morgan;

pub use descriptors::{LigandDescriptors, RuleOfFive, describe};
pub use morgan::{FINGERPRINT_BITS, Fingerprint, cluster, fingerprint, tanimoto};

#[cfg(test)]
#[path = "library_tests.rs"]
mod tests;
