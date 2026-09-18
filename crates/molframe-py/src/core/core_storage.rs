//! Native coordinate blocks with read-only `NumPy` views.

use crate::core_values::PyAabb;
use numpy::ndarray::ArrayView2;
use numpy::{PyArray2, PyArrayMethods};
use pyo3::prelude::*;

#[pyclass(name = "CoordinateBlock", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCoordinateBlock(pub(crate) molframe::CoordinateBlock);

impl PyCoordinateBlock {
    pub(crate) fn view<'py>(
        &self,
        py: Python<'py>,
        range: Option<(u32, u32)>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let owner = Bound::new(py, self.clone())?;
        let (length, pointer) = {
            let snapshot = owner.borrow();
            let positions = match range {
                Some((start, end)) => snapshot.0.range(start..end).ok_or_else(|| {
                    pyo3::exceptions::PyIndexError::new_err("coordinate range is outside the block")
                })?,
                None => snapshot.0.as_slice(),
            };
            (positions.len(), positions.as_ptr().cast::<f32>())
        };
        let shape = (length, 3);
        // SAFETY: CoordinateBlock stores contiguous `[f32; 3]` values.  The
        // owner is retained as NumPy's base object and the view is read-only.
        let view = unsafe { ArrayView2::from_shape_ptr(shape, pointer) };
        // SAFETY: `view` points into the cloned block retained by `owner`.
        let array = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
        array.readwrite().make_nonwriteable();
        Ok(array)
    }
}

#[pymethods]
impl PyCoordinateBlock {
    #[new]
    fn new() -> Self {
        Self(molframe::CoordinateBlock::new())
    }

    #[staticmethod]
    fn with_capacity(positions: usize) -> Self {
        Self(molframe::CoordinateBlock::with_capacity(positions))
    }

    #[getter]
    fn len(&self) -> u32 {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn push(&mut self, position: [f32; 3]) {
        self.0.push(position);
    }

    #[getter]
    fn positions<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.view(py, None)
    }

    fn as_slice<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.view(py, None)
    }

    fn range<'py>(
        &self,
        py: Python<'py>,
        start: u32,
        end: u32,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        self.view(py, Some((start, end)))
    }

    fn bounds(&self, start: u32, end: u32) -> PyAabb {
        PyAabb(self.0.bounds(start..end))
    }

    fn allocated_bytes(&self) -> usize {
        self.0.allocated_bytes()
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCoordinateBlock>()?;
    Ok(())
}
