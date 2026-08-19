//! Native CIF parser and reader entry points.

use crate::cif_document::PyCifDocument;
use crate::errors::read_error;
use crate::io::{PyReadOptions, PyReadReport};
use pyo3::prelude::*;

#[pyclass(name = "CifReader", frozen, from_py_object)]
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct PyCifReader;

#[pymethods]
impl PyCifReader {
    #[new]
    const fn new() -> Self {
        Self
    }

    #[staticmethod]
    fn parse(py: Python<'_>, data: Vec<u8>) -> PyResult<(PyCifDocument, Vec<String>)> {
        parse_native(py, data)
    }

    #[staticmethod]
    #[pyo3(signature = (data, options=None))]
    fn read(
        py: Python<'_>,
        data: Vec<u8>,
        options: Option<&PyReadOptions>,
    ) -> PyResult<PyReadReport> {
        read_native(py, data, options)
    }

    #[staticmethod]
    #[pyo3(signature = (data, options=None))]
    fn read_with_document(
        py: Python<'_>,
        data: Vec<u8>,
        options: Option<&PyReadOptions>,
    ) -> PyResult<(PyCifDocument, PyReadReport)> {
        read_with_document_native(py, data, options)
    }
}

#[pyfunction(name = "parse")]
pub(crate) fn parse_cif(py: Python<'_>, data: Vec<u8>) -> PyResult<(PyCifDocument, Vec<String>)> {
    parse_native(py, data)
}

#[pyfunction(name = "cif_read")]
#[pyo3(signature = (data, options=None))]
pub(crate) fn read_cif(
    py: Python<'_>,
    data: Vec<u8>,
    options: Option<&PyReadOptions>,
) -> PyResult<PyReadReport> {
    read_native(py, data, options)
}

#[pyfunction(name = "read_with_document")]
#[pyo3(signature = (data, options=None))]
pub(crate) fn read_cif_with_document(
    py: Python<'_>,
    data: Vec<u8>,
    options: Option<&PyReadOptions>,
) -> PyResult<(PyCifDocument, PyReadReport)> {
    read_with_document_native(py, data, options)
}

#[pyfunction]
pub(crate) fn split_tag(tag: &str) -> (&str, &str) {
    pdbiox::cif::split_tag(tag)
}

fn parse_native(py: Python<'_>, data: Vec<u8>) -> PyResult<(PyCifDocument, Vec<String>)> {
    py.detach(move || pdbiox::cif::parse(&pdbiox::InputBuffer::from_bytes(data)))
        .map(|(document, findings)| {
            (
                PyCifDocument::from(document),
                findings
                    .into_iter()
                    .map(|finding| finding.to_string())
                    .collect(),
            )
        })
        .map_err(|findings| read_error(py, &findings))
}

fn read_native(
    py: Python<'_>,
    data: Vec<u8>,
    options: Option<&PyReadOptions>,
) -> PyResult<PyReadReport> {
    let options = options.map_or_else(pdbiox::ReadOptions::new, |value| value.0.clone());
    py.detach(move || pdbiox::cif::read(&pdbiox::InputBuffer::from_bytes(data), &options))
        .map(|(structure, findings)| PyReadReport {
            structure: crate::structure::PyStructure::new(structure),
            findings: findings
                .into_iter()
                .map(|finding| finding.to_string())
                .collect(),
        })
        .map_err(|findings| read_error(py, &findings))
}

fn read_with_document_native(
    py: Python<'_>,
    data: Vec<u8>,
    options: Option<&PyReadOptions>,
) -> PyResult<(PyCifDocument, PyReadReport)> {
    let options = options.map_or_else(pdbiox::ReadOptions::new, |value| value.0.clone());
    py.detach(move || {
        pdbiox::cif::read_with_document(&pdbiox::InputBuffer::from_bytes(data), &options)
    })
    .map(|(document, structure, findings)| {
        (
            PyCifDocument::from(document),
            PyReadReport {
                structure: crate::structure::PyStructure::new(structure),
                findings: findings
                    .into_iter()
                    .map(|finding| finding.to_string())
                    .collect(),
            },
        )
    })
    .map_err(|findings| read_error(py, &findings))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCifReader>()?;
    module.add_function(wrap_pyfunction!(parse_cif, module)?)?;
    module.add_function(wrap_pyfunction!(read_cif, module)?)?;
    module.add_function(wrap_pyfunction!(read_cif_with_document, module)?)?;
    module.add_function(wrap_pyfunction!(split_tag, module)?)?;
    Ok(())
}
