//! Rooted phylogenetic trees and deterministic Newick I/O.

mod model;

pub use model::{LadderDirection, NewickError, RerootError, TraversalOrder, Tree};
