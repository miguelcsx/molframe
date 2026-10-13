//! Backbone frames and helix geometry from ordered Cα traces.
//!
//! Complexity: O(n) time and space for n residues; helix fitting adds a fixed-size
//! eigen solve.

mod geometry;

pub use geometry::*;
