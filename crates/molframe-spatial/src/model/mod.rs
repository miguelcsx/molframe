//! Public spatial value model.

mod backend_names;
mod columns;
mod types;

pub use columns::PairColumns;

pub use types::{NeighborPair, SpatialBackend, SpatialError, SpatialOption};

#[cfg(test)]
mod columns_tests;
