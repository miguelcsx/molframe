//! Declarative path-based trajectory I/O.

mod accounted;
mod error;
mod model;
mod read;
mod stream;
mod write;

pub use accounted::read_trajectory_in;
pub use error::TrajectoryIoError;
pub use model::{
    AmberAsciiReadOptions, FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryMetadata,
    TrajectoryReadOptions, TrajectoryReaderOptions, TrajectoryWriteOptions, TrzWriteOptions,
};
pub use read::read_trajectory_materialized;
pub use stream::read_trajectory;
pub use write::write_trajectory;

#[cfg(test)]
#[path = "dispatch_tests.rs"]
mod tests;
