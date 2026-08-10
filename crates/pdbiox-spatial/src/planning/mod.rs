//! Spatial planning and structure-bound resolution.

pub(crate) mod options;
pub(crate) mod resolver;

pub use options::{
    AutoBackendProfile, CellGridOptions, KdPeriodicOptions, NeighborListOptions,
    NeighborSkinProfile, SpatialPlan, SpatialSearchOptions,
};
pub use resolver::StructureSpatial;
