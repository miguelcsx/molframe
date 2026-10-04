//! Stateless geometry over `NumPy` coordinate arrays.

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::prelude::*;

pub(crate) fn coordinates<'a>(array: &'a PyReadonlyArray2<'_, f32>) -> PyResult<&'a [[f32; 3]]> {
    let shape = array.shape();
    if shape.len() != 2 || shape[1] != 3 {
        return Err(crate::error::value("coordinates must have shape (n, 3)"));
    }
    let contiguous = array.as_slice().map_err(|_| {
        crate::error::value("coordinates must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    Ok(contiguous.as_chunks::<3>().0)
}

#[pyfunction]
pub(crate) fn centroid(array: &Bound<'_, PyArray2<f32>>) -> PyResult<Option<[f64; 3]>> {
    let array = array.readonly();
    Ok(molframe::geometry::centroid(coordinates(&array)?))
}

#[pyfunction]
pub(crate) fn rmsd(
    mobile: &Bound<'_, PyArray2<f32>>,
    reference: &Bound<'_, PyArray2<f32>>,
) -> PyResult<f64> {
    let mobile = mobile.readonly();
    let reference = reference.readonly();
    molframe::geometry::rmsd(coordinates(&mobile)?, coordinates(&reference)?)
        .map_err(crate::error::kernel)
}

#[pyfunction]
pub(crate) fn distance_matrix<'py>(
    py: Python<'py>,
    array: &Bound<'_, PyArray2<f32>>,
) -> PyResult<Bound<'py, PyArray2<f64>>> {
    let array = array.readonly();
    let positions = coordinates(&array)?;
    let matrix = py
        .detach(|| molframe::geometry::distance_matrix(positions))
        .map_err(crate::error::kernel)?;
    let rows = matrix.rows();
    matrix.into_values().into_pyarray(py).reshape((rows, rows))
}
