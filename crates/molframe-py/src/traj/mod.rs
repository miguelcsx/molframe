//! Declarative trajectory bindings over native Rust format dispatch.

mod amber_formats;
mod amber_netcdf_formats;
mod arrays;
mod binary_formats;
mod container_formats;
mod data;
mod dlpoly_formats;
pub(crate) mod dms;
pub(crate) mod dms_models;
mod formats;
mod gamess;
mod generic;
mod gromacs_formats;
mod h5md_formats;
mod imd;
mod interpolation;
mod io;
mod kernels;
mod model;
mod namd_formats;
mod neighbors;
mod operations;
mod reader_types;
mod registration;
mod selection;
mod stream_frame;
mod streaming;
mod text_extended;
mod text_formats;
mod topology_extended;
mod topology_formats;
mod transforms;
mod types;
mod xvg;

pub(crate) use amber_formats::PyAmberRestartLayout;
pub(crate) use model::PyTrajectory;
pub(crate) use types::{PyTrajectoryFormat, PyTrajectoryUnits, PyTrajectoryWriteOptions};

pub(crate) use interpolation::register as register_interpolation;
pub(crate) use registration::register;

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;
