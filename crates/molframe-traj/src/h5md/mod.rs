//! H5MD trajectory I/O.

mod error;
mod options;
mod read;
mod schema;
mod write;

pub use error::H5mdError;
pub use options::{H5mdOptions, H5mdUnitSystem};
pub use read::{
    H5mdMetadata, H5mdTrajectory, parse_h5md, parse_h5md_record_with_options,
    parse_h5md_with_options,
};
pub use write::{write_h5md, write_h5md_with_metadata};

#[cfg(test)]
#[path = "h5md_tests.rs"]
mod tests;
