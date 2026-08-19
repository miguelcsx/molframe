//! Explicit document-to-structure lowering.

use crate::cif_document::PyCifDocument;
use crate::errors::read_error;
use crate::io::{PyReadOptions, PyReadReport};
use pyo3::prelude::*;

#[pyfunction(name = "cif_lower")]
#[pyo3(signature = (document, options=None))]
pub(crate) fn lower(
    py: Python<'_>,
    document: &PyCifDocument,
    options: Option<&PyReadOptions>,
) -> PyResult<PyReadReport> {
    let document = document.inner.clone();
    let options = options.map_or_else(pdbiox::ReadOptions::new, |value| value.0.clone());
    py.detach(move || pdbiox::cif::lower(&document, &options))
        .map(|(structure, findings)| PyReadReport {
            structure: crate::structure::PyStructure::new(structure),
            findings: findings
                .into_iter()
                .map(|finding| finding.to_string())
                .collect(),
        })
        .map_err(|findings| read_error(py, &findings))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(lower, module)?)?;
    Ok(())
}
