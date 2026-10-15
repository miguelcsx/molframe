//! AMBER `NetCDF` trajectory I/O.

mod error;
mod read;
mod schema;
mod write;

pub use error::AmberNetcdfError;
pub use read::{
    AmberNetcdfMetadata, AmberNetcdfTrajectory, parse_amber_netcdf, parse_amber_netcdf_record,
};
pub use write::{AmberNetcdfPrecision, AmberNetcdfWriteOptions, write_amber_netcdf};

#[cfg(test)]
#[path = "amber_netcdf_tests.rs"]
mod tests;
