//! Deterministic spatial indices and neighbour search.

#![forbid(unsafe_code)]

mod brute;
mod cell;
mod kd;
mod neighbor;
mod periodic;
mod resolver;
mod search;
mod types;

pub use cell::CellList;
pub use kd::KdTree;
pub use neighbor::NeighborList;
pub use periodic::PeriodicBox;
pub use resolver::StructureSpatial;
pub use search::{pairs_within, within};
pub use types::{NeighborPair, SpatialBackend, SpatialError};
