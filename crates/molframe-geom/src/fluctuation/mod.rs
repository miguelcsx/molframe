//! Root-mean-square fluctuation of atoms across a set of frames.
//!
//! Complexity: O(frames x atoms) time; O(atoms) space for the accumulator.

mod rmsf;

pub use rmsf::*;
