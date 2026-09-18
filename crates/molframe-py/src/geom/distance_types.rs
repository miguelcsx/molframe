//! Dense distance-matrix ownership and zero-copy `NumPy` projection.

use super::arrays::{borrowed_coordinates, matrix_error};
use ::molframe;
use numpy::ndarray::ArrayView2;
use numpy::{PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

/// Retains the immutable Rust allocation used by one exported array.
#[pyclass(frozen)]
struct DistanceMatrixOwner {
    matrix: Arc<molframe::DistanceMatrix>,
}

/// An immutable, row-major dense matrix of pairwise distances.
#[pyclass(name = "DistanceMatrix", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDistanceMatrix(pub(crate) Arc<molframe::DistanceMatrix>);

impl PyDistanceMatrix {
    pub(crate) fn from_native(matrix: molframe::DistanceMatrix) -> Self {
        Self(Arc::new(matrix))
    }
}

#[pymethods]
impl PyDistanceMatrix {
    #[staticmethod]
    fn from_positions(py: Python<'_>, positions: PyReadonlyArray2<'_, f32>) -> PyResult<Self> {
        let positions = borrowed_coordinates(&positions)?;
        py.detach(move || molframe::distance_matrix(positions))
            .map(Self::from_native)
            .map_err(matrix_error)
    }

    #[staticmethod]
    fn from_sets(
        py: Python<'_>,
        left: PyReadonlyArray2<'_, f32>,
        right: PyReadonlyArray2<'_, f32>,
    ) -> PyResult<Self> {
        let left = borrowed_coordinates(&left)?;
        let right = borrowed_coordinates(&right)?;
        py.detach(move || molframe::distance_matrix_between(left, right))
            .map(Self::from_native)
            .map_err(matrix_error)
    }

    #[getter]
    fn rows(&self) -> usize {
        self.0.rows()
    }

    #[getter]
    fn columns(&self) -> usize {
        self.0.columns()
    }

    fn get(&self, row: usize, column: usize) -> Option<f64> {
        self.0.get(row, column)
    }

    /// Materialises a Python list for compatibility with scalar-oriented code.
    fn as_slice(&self) -> Vec<f64> {
        self.0.as_slice().to_vec()
    }

    /// Returns a read-only zero-copy `float64[rows, columns]` `NumPy` view.
    fn as_array<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        numpy_matrix(py, Arc::clone(&self.0))
    }
}

fn numpy_matrix(
    py: Python<'_>,
    matrix: Arc<molframe::DistanceMatrix>,
) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let owner = Bound::new(py, DistanceMatrixOwner { matrix })?;
    let (rows, columns, length, pointer) = {
        let owner = owner.borrow();
        let matrix = owner.matrix.as_ref();
        (
            matrix.rows(),
            matrix.columns(),
            matrix.as_slice().len(),
            matrix.as_slice().as_ptr(),
        )
    };
    if rows.checked_mul(columns) != Some(length) {
        return Err(PyValueError::new_err(
            "distance-matrix buffer length does not match its declared shape",
        ));
    }
    // SAFETY: the exact row-major element count was validated above; the Arc
    // allocation is aligned for `f64` and remains alive through `owner`.
    let view = unsafe { ArrayView2::from_shape_ptr((rows, columns), pointer) };
    // SAFETY: the owner becomes NumPy's base object, retaining the immutable
    // Arc for at least as long as the view. The result is read-only below.
    let result = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
    let _readonly = result.readwrite().make_nonwriteable();
    Ok(result)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDistanceMatrix>()
}
