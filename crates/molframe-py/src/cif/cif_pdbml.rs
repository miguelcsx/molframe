//! Native PDBML/XML conversion over the shared CIF document model.

use crate::cif_document::PyCifDocument;
use crate::errors::read_error;
use crate::io::{PyReadOptions, PyReadReport};
use pyo3::create_exception;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

create_exception!(_native, PdbmlError, PyValueError);
create_exception!(_native, PdbmlReadError, PyValueError);

#[pyfunction]
pub(crate) fn parse_pdbml_document(py: Python<'_>, bytes: Vec<u8>) -> PyResult<PyCifDocument> {
    py.detach(move || -> PyResult<PyCifDocument> {
        molframe::cif::parse_pdbml_document(&bytes)
            .map(Into::into)
            .map_err(|error| PdbmlError::new_err(error.to_string()))
    })
}

#[pyfunction]
#[pyo3(signature = (bytes, options=None))]
pub(crate) fn read_pdbml(
    py: Python<'_>,
    bytes: Vec<u8>,
    options: Option<&PyReadOptions>,
) -> PyResult<(PyCifDocument, PyReadReport)> {
    let options = options.map_or_else(molframe::ReadOptions::new, |value| value.0.clone());
    py.detach(move || molframe::cif::read_pdbml(&bytes, &options))
        .map(|(document, structure, findings)| {
            (
                PyCifDocument::from(document),
                PyReadReport {
                    structure: crate::structure::PyStructure::new(structure),
                    findings: findings.into_iter().map(Into::into).collect(),
                },
            )
        })
        .map_err(|error| match error {
            molframe::cif::PdbmlReadError::Findings(findings) => read_error(py, &findings),
            molframe::cif::PdbmlReadError::Pdbml(error) => {
                PdbmlReadError::new_err(error.to_string())
            }
            _ => PdbmlReadError::new_err("unknown PDBML read failure"),
        })
}

#[pyfunction]
pub(crate) fn write_pdbml(py: Python<'_>, document: &PyCifDocument) -> PyResult<String> {
    py.detach(move || -> PyResult<String> {
        molframe::cif::write_pdbml(&document.inner)
            .map_err(|error| PdbmlError::new_err(error.to_string()))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("PdbmlError", module.py().get_type::<PdbmlError>())?;
    module.add("PdbmlReadError", module.py().get_type::<PdbmlReadError>())?;
    module.add_function(wrap_pyfunction!(parse_pdbml_document, module)?)?;
    module.add_function(wrap_pyfunction!(read_pdbml, module)?)?;
    module.add_function(wrap_pyfunction!(write_pdbml, module)?)?;
    Ok(())
}
