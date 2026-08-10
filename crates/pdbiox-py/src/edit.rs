//! Scoped `NumPy` coordinate mutation backed by the Rust transaction API.

use crate::errors::read_error;
use crate::structure::PyStructure;
use numpy::ndarray::Array2;
use numpy::{PyArray2, PyArrayMethods};
use pdbiox::{CoordinateEditor, ModelIndex};
use pyo3::exceptions::PyValueError;
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
        let (position_count, flattened) = {
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
            let mut flattened = Vec::with_capacity(scalar_count);
            for position in positions.iter() {
                flattened.extend_from_slice(position);
            }
            (positions.len(), flattened)
        };
        let owned = Array2::from_shape_vec((position_count, 3), flattened)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let array = PyArray2::from_owned_array(slf.py(), owned);
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
        if let Some(array) = &self.array {
            let array = array.bind(py).readonly();
            let values = array
                .as_slice()
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
            let positions = self
                .editor
                .try_positions_mut(ModelIndex::new(0))
                .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))?;
            for (target, source) in positions.iter_mut().zip(values.chunks_exact(3)) {
                target.copy_from_slice(source);
            }
        }
        let edited = self
            .editor
            .snapshot()
            .map_err(|findings| read_error(py, &findings))?;
        self.owner.borrow_mut(py).replace(edited);
        Ok(false)
    }
}
