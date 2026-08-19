//! Dense `NumPy` geometry calls delegated to shared Rust kernels.

use super::{PyRigid, PySuperposeOptions, PySuperposition};
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyfunction]
pub(crate) fn distance_matrix<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<'_, f32>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let positions = borrowed_coordinates(&positions)?;
    let matrix = py
        .detach(move || pdbiox::distance_matrix(positions))
        .map_err(matrix_error)?;
    distance_matrix_value(py, matrix)
}

/// Converts one native dense matrix into one contiguous `NumPy` allocation.
pub(crate) fn distance_matrix_value(
    py: Python<'_>,
    matrix: pdbiox::DistanceMatrix,
) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let shape = (matrix.rows(), matrix.columns());
    Array2::from_shape_vec(shape, matrix.into_values())
        .map_err(|error| PyValueError::new_err(error.to_string()))
        .map(|array| array.into_pyarray(py))
}

pub(super) fn matrix_error(error: pdbiox::MatrixError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[pyfunction]
pub(crate) fn rmsd(
    py: Python<'_>,
    mobile: &Bound<'_, PyAny>,
    reference: &Bound<'_, PyAny>,
) -> PyResult<f64> {
    rmsd_kernel(py, mobile, reference)
}

/// RMSD over C-contiguous `(n, 3)` coordinate arrays without materialising
/// an intermediate coordinate buffer.
#[pyfunction(name = "rmsd_flat")]
pub(crate) fn rmsd_flat(
    py: Python<'_>,
    mobile: &Bound<'_, PyAny>,
    reference: &Bound<'_, PyAny>,
) -> PyResult<f64> {
    rmsd_kernel(py, mobile, reference)
}

fn rmsd_kernel(
    py: Python<'_>,
    mobile: &Bound<'_, PyAny>,
    reference: &Bound<'_, PyAny>,
) -> PyResult<f64> {
    let mobile = mobile.extract::<PyReadonlyArray2<'_, f32>>()?;
    let reference = reference.extract::<PyReadonlyArray2<'_, f32>>()?;
    let mobile = flat_coordinates(&mobile)?;
    let reference = flat_coordinates(&reference)?;
    py.detach(|| pdbiox::rmsd_flat(mobile, reference))
        .map_err(superpose_error)
}

fn flat_coordinates<'a>(array: &'a PyReadonlyArray2<'_, f32>) -> PyResult<&'a [f32]> {
    if array.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("coordinates must have shape (n, 3)"));
    }
    array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "coordinates must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })
}

#[pyfunction]
pub(crate) fn superpose(
    py: Python<'_>,
    mobile: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
) -> PyResult<PySuperposition> {
    let mobile = borrowed_coordinates(&mobile)?;
    let reference = borrowed_coordinates(&reference)?;
    py.detach(|| pdbiox::superpose(mobile, reference))
        .map(PySuperposition::from)
        .map_err(superpose_error)
}

#[pyfunction]
pub(crate) fn superpose_with_options(
    py: Python<'_>,
    mobile: PyReadonlyArray2<'_, f32>,
    reference: PyReadonlyArray2<'_, f32>,
    options: &PySuperposeOptions,
) -> PyResult<PySuperposition> {
    let mobile = borrowed_coordinates(&mobile)?;
    let reference = borrowed_coordinates(&reference)?;
    py.detach(|| pdbiox::superpose_with_options(mobile, reference, options.0))
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

/// Borrows a C-contiguous `(n, 3)` `NumPy` coordinate array without copying it.
///
/// Callers must retain `array` for the whole native call.  This is suitable for
/// one-shot kernels; persistent indexes require an explicitly owned snapshot.
pub(crate) fn borrowed_coordinates<'a>(
    array: &'a PyReadonlyArray2<'_, f32>,
) -> PyResult<&'a [[f32; 3]]> {
    if array.shape().get(1).copied() != Some(3) {
        return Err(PyValueError::new_err("coordinates must have shape (n, 3)"));
    }
    let values = array.as_slice().map_err(|_| {
        PyValueError::new_err(
            "coordinates must be C-contiguous; pass copy=True explicitly to materialise them",
        )
    })?;
    let (coordinates, remainder) = values.as_chunks::<3>();
    if remainder.is_empty() {
        Ok(coordinates)
    } else {
        Err(PyValueError::new_err("coordinates must have shape (n, 3)"))
    }
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
