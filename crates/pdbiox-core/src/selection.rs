//! Adaptive atom sets represented as runs, sparse positions, or dense bits.

mod merge;
mod operations;
mod representation;
mod runs;

pub use representation::AtomSelection;

#[cfg(test)]
#[path = "selection_tests.rs"]
mod tests;
