//! Scoped `NumPy` coordinate mutation backed by the Rust transaction API.

use crate::errors::read_error;
use crate::structure::PyStructure;
use numpy::ndarray::ArrayView2;
use numpy::{PyArray2, PyArrayMethods};
use pdbiox::{CoordinateEditor, ModelIndex};
use pyo3::prelude::*;

#[pyclass(name = "CoordinateEdit", skip_from_py_object)]
pub(crate) struct PyCoordinateEdit {
    owner: Py<PyStructure>,
    editor: CoordinateEditor,
    array: Option<Py<PyArray2<f32>>>,
}

impl PyCoordinateEdit {
    fn new(owner: Py<PyStructure>, editor: CoordinateEditor) -> Self {
        Self {
            owner,
            editor,
            array: None,
        }
    }
}

#[pymethods]
impl PyStructure {
    fn edit_coordinates(slf: &Bound<'_, Self>) -> PyCoordinateEdit {
        PyCoordinateEdit::new(
            slf.clone().unbind(),
            slf.borrow().structure().edit_coordinates(),
        )
    }
}

#[pymethods]
impl PyCoordinateEdit {
    fn __enter__<'py>(slf: &Bound<'py, Self>) -> PyResult<Bound<'py, PyArray2<f32>>> {
        if let Some(existing) = &slf.borrow().array {
            return Ok(existing.bind(slf.py()).clone());
        }
        let (position_count, pointer) = {
            let mut context = slf.borrow_mut();
            let positions = context
                .editor
                .try_positions_mut(ModelIndex::new(0))
                .map_err(|finding| read_error(slf.py(), std::slice::from_ref(&finding)))?;
            (positions.len(), positions.as_ptr().cast::<f32>())
        };
        let shape = (position_count, 3);
        // SAFETY: the editor owns a fixed-size contiguous f32 triple buffer and
        // `slf` becomes the NumPy base object. The editor is retained after the
        // scope, so the pointer remains valid even if the array escapes.
        let view = unsafe { ArrayView2::from_shape_ptr(shape, pointer) };
        // SAFETY: the view points into the editor retained by the base object;
        // its allocation cannot move while the NumPy array exists.
        let array = unsafe { PyArray2::borrow_from_array(&view, slf.clone().into_any()) };
        slf.borrow_mut().array = Some(array.clone().unbind());
        Ok(array)
    }

    #[pyo3(signature = (exc_type, _exc_value, _traceback))]
    fn __exit__(
        &mut self,
        py: Python<'_>,
        exc_type: Option<&Bound<'_, PyAny>>,
        _exc_value: Option<&Bound<'_, PyAny>>,
        _traceback: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<bool> {
        if let Some(array) = &self.array {
            let _readonly = array.bind(py).readwrite().make_nonwriteable();
        }
        if exc_type.is_some() {
            return Ok(false);
        }
        let edited = self
            .editor
            .snapshot()
            .map_err(|findings| read_error(py, &findings))?;
        self.owner.borrow_mut(py).replace(edited);
        Ok(false)
    }
}
