//! Dense `NumPy` geometry calls delegated to shared Rust kernels.

use super::{PyRigid, PySuperposition};
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn distance_matrix<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let positions = coordinates(positions)?;
    let matrix = py
        .detach(move || pdbiox::distance_matrix(&positions))
        .map_err(matrix_error)?;
    Array2::from_shape_vec(
        (matrix.rows(), matrix.columns()),
        matrix.as_slice().to_vec(),
    )
    .map_err(|error| PyValueError::new_err(error.to_string()))
    .map(|array| array.into_pyarray(py))
}

pub(super) fn matrix_error(error: pdbiox::MatrixError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[pyfunction]
pub(crate) fn rmsd(
    py: Python<'_>,
    mobile: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
) -> PyResult<f64> {
    let mobile = coordinates(mobile)?;
    let reference = coordinates(reference)?;
    py.detach(|| pdbiox::rmsd(&mobile, &reference))
        .map_err(superpose_error)
}

#[pyfunction]
pub(crate) fn superpose(
    py: Python<'_>,
    mobile: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
) -> PyResult<PySuperposition> {
    let mobile = coordinates(mobile)?;
    let reference = coordinates(reference)?;
    py.detach(|| pdbiox::superpose(&mobile, &reference))
        .map(PySuperposition::from)
        .map_err(superpose_error)
}

pub(crate) fn coordinates(array: PyReadonlyArray2<'_, f32>) -> PyResult<Vec<[f32; 3]>> {
    if array.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("coordinates must have shape (n, 3)"));
    }
    let coordinates = array
        .as_array()
        .rows()
        .into_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    drop(array);
    Ok(coordinates)
}

fn superpose_error(error: pdbiox::SuperposeError) -> PyErr {
    let message = match error {
        pdbiox::SuperposeError::LengthMismatch => "coordinate sets differ in length",
        pdbiox::SuperposeError::TooFewPoints => "superposition requires at least three points",
        pdbiox::SuperposeError::Degenerate => "superposition is undefined for collinear points",
        pdbiox::SuperposeError::InvalidOptions => "superposition options are invalid",
        pdbiox::SuperposeError::TooManyPoints => "coordinate set exceeds the supported size",
        pdbiox::SuperposeError::Eigen(_) => "superposition eigendecomposition failed",
    };
    PyValueError::new_err(message)
}

impl From<pdbiox::Superposition> for PySuperposition {
    fn from(value: pdbiox::Superposition) -> Self {
        Self {
            transform: PyRigid(value.transform),
            rmsd: value.rmsd,
        }
    }
}
