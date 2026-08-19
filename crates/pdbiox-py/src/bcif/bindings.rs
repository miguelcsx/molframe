//! Lazy `BinaryCIF` container access and structure decoding.

use crate::cif_document::{PyCifCategory, PyCifDocument};
use crate::errors::read_error;
use crate::io::{PyLimits, PyReadOptions, PyReadReport};
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[pyclass(name = "BinaryDocument", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBinaryDocument(pub(crate) pdbiox::BinaryDocument);

#[pymethods]
impl PyBinaryDocument {
    #[new]
    fn new(py: Python<'_>, data: &Bound<'_, PyBytes>, limits: &PyLimits) -> PyResult<Self> {
        let input = pdbiox::InputBuffer::from_bytes(data.as_bytes().to_vec());
        pdbiox::bcif::read_document(&input, limits.0)
            .map(Self)
            .map_err(|finding| read_error(py, &[finding]))
    }

    #[getter]
    fn version(&self) -> &str {
        self.0.version()
    }

    #[getter]
    fn encoder(&self) -> &str {
        self.0.encoder()
    }

    #[getter]
    fn block_count(&self) -> usize {
        self.0.block_count()
    }

    fn category(
        &self,
        py: Python<'_>,
        block: usize,
        name: &str,
    ) -> PyResult<Option<PyCifCategory>> {
        self.0
            .category(block, name)
            .map(|value| value.map(|inner| PyCifCategory { inner }))
            .map_err(|finding| read_error(py, &[finding]))
    }

    fn to_document(&self, py: Python<'_>) -> PyResult<PyCifDocument> {
        self.0
            .to_document()
            .map(|inner| PyCifDocument { inner })
            .map_err(|finding| read_error(py, &[finding]))
    }
}

#[pyclass(name = "BcifReader", frozen, skip_from_py_object)]
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct PyBcifReader;

#[pymethods]
impl PyBcifReader {
    #[staticmethod]
    fn read_document(
        py: Python<'_>,
        data: &Bound<'_, PyBytes>,
        limits: &PyLimits,
    ) -> PyResult<PyBinaryDocument> {
        PyBinaryDocument::new(py, data, limits)
    }

    #[staticmethod]
    fn read(
        py: Python<'_>,
        data: &Bound<'_, PyBytes>,
        options: &PyReadOptions,
    ) -> PyResult<PyReadReport> {
        let input = pdbiox::InputBuffer::from_bytes(data.as_bytes().to_vec());
        pdbiox::bcif::read(&input, &options.0)
            .map(|(structure, findings)| PyReadReport {
                structure: crate::structure::PyStructure::new(structure),
                findings: findings
                    .into_iter()
                    .map(|value| value.to_string())
                    .collect(),
            })
            .map_err(|findings| read_error(py, &findings))
    }
}

#[pyfunction]
pub(crate) fn read_bcif_document(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    limits: &PyLimits,
) -> PyResult<PyBinaryDocument> {
    PyBinaryDocument::new(py, data, limits)
}

#[pyfunction]
pub(crate) fn read_bcif(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    options: &PyReadOptions,
) -> PyResult<PyReadReport> {
    PyBcifReader::read(py, data, options)
}

#[pyfunction]
pub(crate) fn read_bcif_with_document(
    py: Python<'_>,
    data: &Bound<'_, PyBytes>,
    options: &PyReadOptions,
) -> PyResult<(PyCifDocument, PyReadReport)> {
    let input = pdbiox::InputBuffer::from_bytes(data.as_bytes().to_vec());
    let options = options.0.clone();
    py.detach(move || pdbiox::bcif::read_with_document(&input, &options))
        .map(|(document, structure, findings)| {
            (
                PyCifDocument { inner: document },
                PyReadReport {
                    structure: crate::structure::PyStructure::new(structure),
                    findings: findings
                        .into_iter()
                        .map(|value| value.to_string())
                        .collect(),
                },
            )
        })
        .map_err(|findings| read_error(py, &findings))
}

#[pyfunction]
pub(crate) fn write_bcif_document<'py>(
    py: Python<'py>,
    document: &PyCifDocument,
) -> PyResult<Bound<'py, PyBytes>> {
    let document = document.inner.clone();
    py.detach(move || pdbiox::bcif::write_document(&document))
        .map(|data| PyBytes::new(py, &data))
        .map_err(|finding| read_error(py, &[finding]))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::bcif_codec::register(module)?;
    module.add_class::<PyBinaryDocument>()?;
    module.add_class::<PyBcifReader>()?;
    module.add_function(wrap_pyfunction!(read_bcif_document, module)?)?;
    module.add_function(wrap_pyfunction!(read_bcif, module)?)?;
    module.add_function(wrap_pyfunction!(read_bcif_with_document, module)?)?;
    module.add_function(wrap_pyfunction!(write_bcif_document, module)?)?;
    Ok(())
}
