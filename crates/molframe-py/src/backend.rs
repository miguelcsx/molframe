//! Spatial backend names accepted by the Python contract.

use pyo3::prelude::*;

/// The backend a name selects, through the spelling Rust owns.
pub(crate) fn parse(value: &str) -> PyResult<molframe::spatial::SpatialBackend> {
    value.parse().map_err(crate::error::kernel)
}
