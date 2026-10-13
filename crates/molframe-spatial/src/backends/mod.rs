//! Spatial index implementations.

pub(crate) mod brute;
mod brute_simd;
pub(crate) mod cell;
pub(crate) mod kd;
pub(crate) mod neighbor;

pub use cell::CellList;
pub use kd::{KdTree, KnnOptions, KnnResult, NO_NEIGHBOR};
pub use neighbor::NeighborList;
