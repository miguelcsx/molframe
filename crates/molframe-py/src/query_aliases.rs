//! Python view of named query definitions.

use crate::bindings::{PyQuery, PyStructure};
use pyo3::prelude::*;

use molframe::CompletionKind;

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

/// Complete a query at a UTF-8 byte cursor using the native language registry.
///
/// `structure` is optional; when supplied, identifier-valued columns such as
/// `chain`, `resname`, and `name` are completed from the structure's own
/// native indexes rather than any client-maintained vocabulary.
#[pyfunction]
#[pyo3(signature = (source, cursor, aliases = None, structure = None))]
pub(crate) fn complete(
    source: &str,
    cursor: usize,
    aliases: Option<&PyQueryAliases>,
    structure: Option<&PyStructure>,
) -> (usize, usize, Vec<(String, String)>) {
    let empty = molframe::QueryAliases::default();
    let aliases = aliases.map_or(&empty, |value| &value.inner);
    let values = structure.map(|structure| &structure.inner as &dyn molframe::StructureValues);
    let result = molframe::complete(source, cursor, aliases, values);
    let items = result
        .items
        .into_iter()
        .map(|item| {
            let kind = match item.kind {
                CompletionKind::Keyword => "keyword",
                CompletionKind::Column => "column",
                CompletionKind::Macro => "macro",
                CompletionKind::Alias => "alias",
                CompletionKind::Value => "value",
            };
            (item.label, kind.to_owned())
        })
        .collect();
    (result.replacement_start, result.replacement_end, items)
}
