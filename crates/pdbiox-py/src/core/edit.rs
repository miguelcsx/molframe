//! Scoped `NumPy` coordinate mutation backed by the Rust transaction API.

use crate::errors::read_error;
use crate::structure::PyStructure;
use numpy::ndarray::ArrayViewMut2;
use numpy::{PyArray2, PyArrayMethods};
use pdbiox::{CoordinateEditor, ModelIndex};
use pyo3::exceptions::{PyMemoryError, PyValueError};
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
    fn edit_coordinates(
        slf: &Bound<'_, Self>,
        context: &crate::core::execution::PyExecutionContext,
    ) -> PyResult<PyCoordinateEdit> {
        slf.borrow()
            .structure()
            .edit_coordinates(&context.native())
            .map(|editor| PyCoordinateEdit::new(slf.clone().unbind(), editor))
            .map_err(|error| PyMemoryError::new_err(error.to_string()))
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
            let Some(scalar_count) = positions.len().checked_mul(3) else {
                return Err(PyValueError::new_err(
                    "coordinate array shape overflows usize",
                ));
            };
            let _ = scalar_count;
            (positions.len(), positions.as_mut_ptr().cast::<f32>())
        };
        let shape = (position_count, 3);
        // SAFETY: the editor owns a stable transactional coordinate buffer for
        // the whole context. The returned NumPy array retains `slf` as its base,
        // so the editor and its buffer cannot be dropped while the view exists.
        let view = unsafe { ArrayViewMut2::from_shape_ptr(shape, pointer) };
        // SAFETY: `view` covers exactly the editor's dense triples and the
        // editor never reallocates them during the scoped transaction.
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
            self.array = None;
            return Ok(false);
        }
        let edited = self
            .editor
            .snapshot()
            .map_err(|findings| read_error(py, &findings))?;
        self.array = None;
        self.owner.borrow_mut(py).replace(edited);
        Ok(false)
    }
}
