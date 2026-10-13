//! Optimal rigid superposition and deviation measures between corresponding point sets.
//!
//! Complexity: O(n) time and O(1) extra space for n points, plus a fixed-size
//! four-by-four solve; the fit never searches for correspondence.

mod deviation;
mod fit;
mod qcp;

pub use deviation::{rmsd, rmsd_flat};
pub use fit::*;
