//! Spatial backend names accepted by the Python contract.

use pyo3::{exceptions::PyValueError, prelude::*};

pub(crate) fn parse(value: &str) -> PyResult<molframe::spatial::SpatialBackend> {
    match value {
        "auto" => Ok(molframe::spatial::SpatialBackend::Auto),
        "cell" => Ok(molframe::spatial::SpatialBackend::CellList),
        "kd_tree" => Ok(molframe::spatial::SpatialBackend::KdTree),
        "brute_force" => Ok(molframe::spatial::SpatialBackend::BruteForce),
        _ => Err(PyValueError::new_err(
            "backend must be 'auto', 'cell', 'kd_tree', or 'brute_force'",
        )),
    }
}
