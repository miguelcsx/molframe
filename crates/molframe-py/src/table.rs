//! A read-only columnar result: named `NumPy` columns of equal length.

use numpy::{PyArray1, PyArrayMethods, ToPyArray};
use pyo3::{prelude::*, types::PyDict};

/// Named columns of equal length, as read-only arrays.
#[derive(Debug)]
#[pyclass(name = "Table", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyTable {
    length: usize,
    columns: Vec<(&'static str, Py<PyAny>)>,
}

/// Builds a table one typed column at a time.
pub(crate) struct TableBuilder<'py> {
    py: Python<'py>,
    length: usize,
    columns: Vec<(&'static str, Py<PyAny>)>,
}

impl<'py> TableBuilder<'py> {
    pub(crate) const fn new(py: Python<'py>, length: usize) -> Self {
        Self {
            py,
            length,
            columns: Vec::new(),
        }
    }

    /// Adds a `u32` atom-index column.
    pub(crate) fn indices(self, name: &'static str, values: &[u32]) -> Self {
        self.column(name, values)
    }

    /// Adds an `f32` column.
    pub(crate) fn single(self, name: &'static str, values: &[f32]) -> Self {
        self.column(name, values)
    }

    /// Adds an `f64` column.
    pub(crate) fn double(self, name: &'static str, values: &[f64]) -> Self {
        self.column(name, values)
    }

    fn column<T: numpy::Element>(mut self, name: &'static str, values: &[T]) -> Self {
        let array: Bound<'py, PyArray1<T>> = values.to_pyarray(self.py);
        array.readwrite().make_nonwriteable();
        self.columns.push((name, array.into_any().unbind()));
        self
    }

    pub(crate) fn finish(self) -> PyTable {
        PyTable {
            length: self.length,
            columns: self.columns,
        }
    }
}

#[pymethods]
impl PyTable {
    fn __len__(&self) -> usize {
        self.length
    }

    /// Column names in order.
    #[getter]
    fn names(&self) -> Vec<&'static str> {
        self.columns.iter().map(|(name, _)| *name).collect()
    }

    fn __getitem__(&self, py: Python<'_>, name: &str) -> PyResult<Py<PyAny>> {
        self.columns
            .iter()
            .find(|(candidate, _)| *candidate == name)
            .map(|(_, column)| column.clone_ref(py))
            .ok_or_else(|| crate::error::key(name))
    }

    fn __contains__(&self, name: &str) -> bool {
        self.columns.iter().any(|(candidate, _)| *candidate == name)
    }

    /// The columns as a dictionary of arrays.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dictionary = PyDict::new(py);
        for (name, column) in &self.columns {
            dictionary.set_item(name, column.clone_ref(py))?;
        }
        Ok(dictionary)
    }

    /// The columns as an Arrow stream, copied once.
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, pyo3::types::PyCapsule>> {
        crate::interop::columns_stream(py, &self.columns)
    }

    fn __repr__(&self) -> String {
        let names: Vec<&str> = self.columns.iter().map(|(name, _)| *name).collect();
        format!("Table(rows={}, columns={})", self.length, names.join(", "))
    }
}
