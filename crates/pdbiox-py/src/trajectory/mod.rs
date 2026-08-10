//! Declarative trajectory bindings over native Rust format dispatch.

mod arrays;
mod io;
mod model;
mod types;

pub(crate) use model::PyTrajectory;
pub(crate) use types::{PyTrajectoryFormat, PyTrajectoryUnits, PyTrajectoryWriteOptions};

#[cfg(test)]
#[path = "trajectory_tests.rs"]
mod tests;
