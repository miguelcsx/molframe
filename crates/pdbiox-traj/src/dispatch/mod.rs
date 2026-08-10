//! Declarative path-based trajectory I/O.

mod error;
mod model;
mod read;
mod write;

pub use error::TrajectoryIoError;
pub use model::{
    AmberAsciiReadOptions, FormatMetadata, TrajectoryData, TrajectoryFormat, TrajectoryMetadata,
    TrajectoryReadOptions, TrajectoryWriteOptions, TrzWriteOptions,
};
pub use read::read_trajectory;
pub use write::write_trajectory;

#[cfg(test)]
#[path = "dispatch_tests.rs"]
mod tests;
