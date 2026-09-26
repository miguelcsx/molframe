//! Python view of named query definitions.

use crate::bindings::PyQuery;
use pyo3::prelude::*;

/// Named query definitions resolved to closed queries.
///
/// Each method converts its arguments, makes one call into
/// `molframe::QueryAliases`, and converts the result.
#[derive(Clone, Debug, Default)]
#[pyclass(name = "QueryAliases", skip_from_py_object)]
pub(crate) struct PyQueryAliases {
    inner: molframe::QueryAliases,
}

#[pymethods]
impl PyQueryAliases {
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Defines or redefines `name`, returning the definition it replaced.
    fn define(&mut self, name: &str, query: &PyQuery) -> PyResult<Option<PyQuery>> {
        self.inner
            .define(name, query.native().clone())
            .map(|previous| previous.map(PyQuery::from_native))
            .map_err(|finding| {
                crate::bindings::findings_error(&molframe::Findings::from(vec![finding]))
            })
    }

    /// Removes `name`, returning its definition.
    fn remove(&mut self, name: &str) -> Option<PyQuery> {
        self.inner.remove(name).map(PyQuery::from_native)
    }

    /// The definition of `name`, as written.
    fn get(&self, name: &str) -> Option<PyQuery> {
        self.inner.get(name).cloned().map(PyQuery::from_native)
    }

    /// Defined names in ascending order.
    #[getter]
    fn names(&self) -> Vec<String> {
        self.inner.names().map(str::to_owned).collect()
    }

    /// Replaces every `$name` reference with its definition.
    fn resolve(&self, query: &PyQuery) -> PyResult<PyQuery> {
        self.inner
            .resolve(query.native())
            .map(PyQuery::from_native)
            .map_err(|findings| {
                crate::bindings::findings_error(&molframe::Findings::from(findings))
            })
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    fn __contains__(&self, name: &str) -> bool {
        self.inner.contains(name)
    }
}
