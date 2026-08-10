//! Python handles and zero-copy array conversion for a structure snapshot.

use crate::errors::read_error;
use numpy::ndarray::ArrayView2;
use numpy::{PyArray2, PyArrayMethods};
use pdbiox::Structure;
use pyo3::prelude::*;
use std::path::PathBuf;

#[pyclass(name = "Structure", skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyStructure {
    inner: Structure,
}

#[pyclass(frozen)]
struct SnapshotOwner {
    inner: Structure,
}

impl PyStructure {
    pub(crate) const fn new(inner: Structure) -> Self {
        Self { inner }
    }

    pub(crate) const fn structure(&self) -> &Structure {
        &self.inner
    }

    pub(crate) fn replace(&mut self, inner: Structure) {
        self.inner = inner;
    }
}

#[pymethods]
impl PyStructure {
    fn __repr__(&self) -> String {
        self.inner.to_string()
    }

    fn write(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        pdbiox::write(path, &self.inner).map_err(|findings| read_error(py, &findings))
    }

    #[getter]
    fn atom_count(&self) -> u32 {
        self.inner.atom_count()
    }

    #[getter]
    fn model_count(&self) -> usize {
        self.inner.model_count()
    }

    #[getter]
    fn chain_count(&self) -> usize {
        self.inner.chain_count()
    }

    #[getter]
    fn residue_count(&self) -> usize {
        self.inner.residue_count()
    }

    #[getter]
    fn xyz<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let owner = Bound::new(
            py,
            SnapshotOwner {
                inner: self.inner.clone(),
            },
        )?;
        let (position_count, pointer) = {
            let snapshot = owner.borrow();
            let positions = snapshot.inner.positions();
            (positions.len(), positions.as_ptr().cast::<f32>())
        };
        position_count.checked_mul(3).ok_or_else(|| {
            pyo3::exceptions::PyOverflowError::new_err("coordinate array shape overflows usize")
        })?;
        let shape = (position_count, 3);
        // SAFETY: `CoordinateBlock` is an immutable contiguous array of f32
        // triples with no padding; the checked product above proves the shape
        // is representable. `owner` holds the Structure snapshot and becomes the NumPy
        // base object, so the allocation cannot move or be freed while Python
        // can reach this view. The dimensions cover exactly len * 3 elements.
        let view = unsafe { ArrayView2::from_shape_ptr(shape, pointer) };
        // SAFETY: `view` points into the snapshot retained by `owner`, and that
        // snapshot is immutable and never reallocates its coordinate block;
        // the array is made non-writeable before it escapes this function.
        let array = unsafe { PyArray2::borrow_from_array(&view, owner.into_any()) };
        let _readonly = array.readwrite().make_nonwriteable();
        Ok(array)
    }
}

#[cfg(test)]
#[path = "binding_tests.rs"]
mod tests;
