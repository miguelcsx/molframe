//! Validation and complexity metadata for coordinate plan operations.

use crate::geometry::borrowed_coordinates;
use crate::graph::PySpatialBackend;
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

pub(super) fn validate_contacts(left: &str, right: &str, cutoff: f32) -> PyResult<()> {
    if left.is_empty() || right.is_empty() {
        return Err(PyValueError::new_err(
            "left and right selections are required",
        ));
    }
    if !cutoff.is_finite() || cutoff <= 0.0 {
        return Err(PyValueError::new_err("cutoff must be finite and positive"));
    }
    Ok(())
}

/// Validates two borrowed coordinate buffers without materialising either one.
pub(super) fn validate_coordinate_pair(
    py: Python<'_>,
    mobile: &Py<PyAny>,
    reference: &Py<PyAny>,
) -> PyResult<()> {
    let mobile = mobile.bind(py).extract::<PyReadonlyArray2<'_, f32>>()?;
    let mobile_count = borrowed_coordinates(&mobile)?.len();
    let reference = reference.bind(py).extract::<PyReadonlyArray2<'_, f32>>()?;
    let reference_count = borrowed_coordinates(&reference)?.len();
    if mobile_count == reference_count {
        Ok(())
    } else {
        Err(PyValueError::new_err(
            "mobile and reference must contain the same number of coordinates",
        ))
    }
}

pub(super) fn contacts_complexity(backend: PySpatialBackend) -> &'static str {
    match backend {
        PySpatialBackend::BruteForce => "O(n_left * n_right + k)",
        PySpatialBackend::CellList => "O(n_left + n_right + k) expected",
        PySpatialBackend::KdTree => "O(n_right log n_right + n_left log n_right + k) expected",
        PySpatialBackend::NeighborList => "O(n_left + k) amortized after index build",
        PySpatialBackend::Auto => "planner-dependent",
    }
}

pub(super) fn validate_inclusion_radius(value: f64) -> PyResult<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(PyValueError::new_err(
            "inclusion_radius must be finite and positive",
        ))
    }
}
