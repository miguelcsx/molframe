//! Deterministic spatial indices and neighbour search.

#![forbid(unsafe_code)]

mod backends;
mod geometry;
mod model;
mod planning;
mod queries;

pub(crate) use backends::brute;
pub(crate) use geometry::numeric;

pub use backends::{CellList, KdTree, NeighborList};
pub use geometry::{PeriodicBox, PeriodicImage};
pub use model::{NeighborPair, SpatialBackend, SpatialError, SpatialOption};
pub use planning::{
    AutoBackendProfile, CellGridOptions, KdPeriodicOptions, NeighborListOptions,
    NeighborSkinProfile, SpatialPlan, SpatialSearchOptions, StructureSpatial,
};
pub use queries::{
    PairQuery, count_pairs_within, for_each_pairs_within_unsorted, pairs_within,
    pairs_within_unsorted, pairs_within_unsorted_with_options, pairs_within_with_options,
    reduce_pairs_within_unsorted, within, within_with_options,
};
