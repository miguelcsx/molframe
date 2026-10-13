//! Eigen-decomposition of small symmetric matrices of compile-time size N.
//!
//! Complexity: O(N^3) time per sweep and O(N^2) space; N is 3 or 4 in practice.

mod solver;

pub use solver::*;
