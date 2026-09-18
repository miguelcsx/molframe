//! Explanations and reproducible configuration for geometry nodes.

use super::model::PyGeometryOperation;
use crate::geometry::PyEigenOptions;
use pyo3::prelude::*;
use pyo3::types::PyDict;

pub(crate) fn explain_operation<'py>(
    py: Python<'py>,
    operation: &PyGeometryOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation_tag(operation))?;
    result.set_item("execution", "native")?;
    result.set_item("requires_c_contiguous", true)?;
    result.set_item(
        "materializes",
        matches!(
            operation,
            PyGeometryOperation::DistanceMatrix { .. }
                | PyGeometryOperation::DistanceMatrixBetween { .. }
                | PyGeometryOperation::Rmsf { .. }
        ),
    )?;
    result.set_item(
        "complexity",
        match operation {
            PyGeometryOperation::DistanceMatrix { .. }
            | PyGeometryOperation::DistanceMatrixBetween { .. } => "O(n²) time and output",
            PyGeometryOperation::Rmsf { .. } => "O(frames × atoms)",
            _ => "O(n)",
        },
    )?;
    Ok(result)
}

pub(crate) fn operation_to_dict<'py>(
    py: Python<'py>,
    operation: &PyGeometryOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation_tag(operation))?;
    match operation {
        PyGeometryOperation::Centroid { positions }
        | PyGeometryOperation::DistanceMatrix { positions } => {
            result.set_item("positions", positions.bind(py))?;
        }
        PyGeometryOperation::CentreOfMass { positions, masses }
        | PyGeometryOperation::RadiusOfGyration { positions, masses }
        | PyGeometryOperation::InertiaTensor { positions, masses } => {
            result.set_item("positions", positions.bind(py))?;
            result.set_item(
                "masses",
                masses
                    .as_ref()
                    .map_or_else(|| py.None(), |value| value.clone_ref(py)),
            )?;
        }
        PyGeometryOperation::PrincipalAxes {
            positions,
            masses,
            options,
        } => {
            result.set_item("positions", positions.bind(py))?;
            result.set_item(
                "masses",
                masses
                    .as_ref()
                    .map_or_else(|| py.None(), |value| value.clone_ref(py)),
            )?;
            result.set_item("options", Py::new(py, PyEigenOptions { inner: *options })?)?;
        }
        PyGeometryOperation::Asphericity { positions, options }
        | PyGeometryOperation::GyrationAxes { positions, options } => {
            result.set_item("positions", positions.bind(py))?;
            result.set_item("options", Py::new(py, PyEigenOptions { inner: *options })?)?;
        }
        PyGeometryOperation::DistanceMatrixBetween { left, right } => {
            result.set_item("left", left.bind(py))?;
            result.set_item("right", right.bind(py))?;
        }
        PyGeometryOperation::Rmsf { frames } => {
            result.set_item("frames", frames.bind(py))?;
        }
    }
    Ok(result)
}

fn operation_tag(operation: &PyGeometryOperation) -> &'static str {
    match operation {
        PyGeometryOperation::Centroid { .. } => "centroid",
        PyGeometryOperation::CentreOfMass { .. } => "centre_of_mass",
        PyGeometryOperation::RadiusOfGyration { .. } => "radius_of_gyration",
        PyGeometryOperation::InertiaTensor { .. } => "inertia_tensor",
        PyGeometryOperation::PrincipalAxes { .. } => "principal_axes",
        PyGeometryOperation::Asphericity { .. } => "asphericity",
        PyGeometryOperation::GyrationAxes { .. } => "gyration_axes",
        PyGeometryOperation::DistanceMatrix { .. } => "distance_matrix",
        PyGeometryOperation::DistanceMatrixBetween { .. } => "distance_matrix_between",
        PyGeometryOperation::Rmsf { .. } => "rmsf",
    }
}
