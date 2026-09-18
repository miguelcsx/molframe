//! Adaptive atom sets represented as runs, sparse positions, or dense bits.

mod merge;
mod operations;
mod representation;
mod runs;

pub use representation::{AtomSelection, SelectionIter};

#[cfg(test)]
#[path = "query_tests.rs"]
mod tests;
