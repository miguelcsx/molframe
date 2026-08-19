use super::{PyDiffusionMap, PyPcaResult};
use numpy::ndarray::{Array1, Array2};
use numpy::{
    IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3,
    PyUntypedArrayMethods,
};
use pdbiox::traj::{DiffusionMap, PcaResult};
use pyo3::Py;
use pyo3::prelude::*;

pub(crate) fn triples_owned(values: PyReadonlyArray2<'_, f32>) -> PyResult<Box<[[f32; 3]]>> {
    if values.shape().get(1) != Some(&3) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "coordinates must have shape (atoms, 3)",
        ));
    }
    let result = values
        .as_array()
        .outer_iter()
        .map(|row| [row[0], row[1], row[2]])
        .collect();
    drop(values);
    Ok(result)
}

/// Borrows a validated C-contiguous `(frames, atoms, 3)` `NumPy` array.
///
/// The resulting native view retains no Python object itself; callers must keep
/// the `NumPy` owner alive until the Rust kernel returns.
pub(crate) fn borrowed_frame_input<'py>(
    values: &'py PyReadonlyArray3<'py, f32>,
) -> PyResult<pdbiox::FrameInput<'py>> {
    let shape = values.shape();
    if shape.get(2) != Some(&3) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "trajectory coordinates must have shape (frames, atoms, 3)",
        ));
    }
    let values = values.as_slice().map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "trajectory coordinates must be C-contiguous; call numpy.ascontiguousarray explicitly",
        )
    })?;
    let (positions, remainder) = values.as_chunks::<3>();
    if !remainder.is_empty() {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "trajectory coordinates must have shape (frames, atoms, 3)",
        ));
    }
    Ok(pdbiox::FrameInput {
        positions,
        frame_count: shape[0],
        atom_count: shape[1],
    })
}

/// Borrows a validated C-contiguous `(frames, atoms, 3)` `NumPy` array.
///
/// The resulting native view retains no Python object itself; callers must keep
/// the `NumPy` owner alive until the Rust kernel returns.
pub(crate) fn borrowed_frames<'py>(
    values: &'py PyReadonlyArray3<'py, f32>,
) -> PyResult<pdbiox::traj::FrameView<'py>> {
    let input = borrowed_frame_input(values)?;
    pdbiox::traj::FrameView::new(input.positions, input.frame_count, input.atom_count)
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))
}

pub(crate) fn scalar_values(values: PyReadonlyArray1<'_, f32>) -> Box<[f32]> {
    let result = values.as_array().iter().copied().collect();
    drop(values);
    result
}

pub(crate) fn pca_result(py: Python<'_>, value: PcaResult) -> Result<PyPcaResult, &'static str> {
    Ok(PyPcaResult {
        mean: Array1::from_vec(value.mean.into_vec())
            .into_pyarray(py)
            .unbind(),
        eigenvalues: Array1::from_vec(value.eigenvalues.into_vec())
            .into_pyarray(py)
            .unbind(),
        components: nested_array2(py, value.components.into_vec())?,
        projections: nested_array2(py, value.projections.into_vec())?,
    })
}

pub(crate) fn diffusion_result(
    py: Python<'_>,
    value: DiffusionMap,
) -> Result<PyDiffusionMap, &'static str> {
    Ok(PyDiffusionMap {
        eigenvalues: Array1::from_vec(value.eigenvalues.into_vec())
            .into_pyarray(py)
            .unbind(),
        coordinates: nested_array2(py, value.coordinates.into_vec())?,
        graph_components: Array1::from_vec(value.graph_components.into_vec())
            .into_pyarray(py)
            .unbind(),
    })
}

fn nested_array2(
    py: Python<'_>,
    values: Vec<Box<[f64]>>,
) -> Result<Py<PyArray2<f64>>, &'static str> {
    let rows = values.len();
    let columns = values.first().map_or(0, |row| row.len());
    if values.iter().any(|row| row.len() != columns) {
        return Err("native result is not rectangular");
    }
    Array2::from_shape_vec((rows, columns), values.into_iter().flatten().collect())
        .map(|array| array.into_pyarray(py).unbind())
        .map_err(|_| "native result has an invalid shape")
}

pub(crate) fn array2<T: numpy::Element + Copy>(
    py: Python<'_>,
    values: &[[T; 3]],
) -> Result<Py<PyArray2<T>>, &'static str> {
    Array2::from_shape_vec(
        (values.len(), 3),
        values.iter().flat_map(|row| row.iter().copied()).collect(),
    )
    .map(|array| array.into_pyarray(py).unbind())
    .map_err(|_| "native result has an invalid shape")
}
