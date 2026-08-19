//! Python adaptations of the facade's input and reader contracts.

use crate::errors::read_error;
use crate::io::{PyLimits, PyReadOptions, PyReadReport};
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use std::io::Cursor;
use std::path::PathBuf;

#[pyclass(name = "Compression", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCompression {
    NoCompression,
    Gzip,
    Zstd,
}

#[pyclass(name = "InputKind", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyInputKind {
    Owned,
    Mapped,
}

#[pyclass(name = "InputBuffer", from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyInputBuffer(pub(crate) pdbiox::InputBuffer);

#[pyclass(name = "Select", frozen, from_py_object)]
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct PySelect;

#[pyclass(name = "Reader", frozen)]
pub(crate) struct PyReader;

impl From<pdbiox::Compression> for PyCompression {
    fn from(value: pdbiox::Compression) -> Self {
        match value {
            pdbiox::Compression::Gzip => Self::Gzip,
            pdbiox::Compression::Zstd => Self::Zstd,
            _ => Self::NoCompression,
        }
    }
}

impl From<pdbiox::InputKind> for PyInputKind {
    fn from(value: pdbiox::InputKind) -> Self {
        match value {
            pdbiox::InputKind::Owned => Self::Owned,
            pdbiox::InputKind::Mapped => Self::Mapped,
        }
    }
}

#[pymethods]
impl PyCompression {
    #[staticmethod]
    fn sniff(data: &Bound<'_, PyBytes>) -> Self {
        pdbiox::Compression::sniff(data.as_bytes()).into()
    }
}

#[pymethods]
impl PyInputBuffer {
    #[staticmethod]
    fn from_bytes(data: &Bound<'_, PyBytes>) -> Self {
        Self(pdbiox::InputBuffer::from_bytes(data.as_bytes().to_vec()))
    }

    #[staticmethod]
    #[pyo3(signature = (path, limits=None))]
    fn open(py: Python<'_>, path: PathBuf, limits: Option<&PyLimits>) -> PyResult<Self> {
        let limits = limits.map_or_else(pdbiox::Limits::default, |value| value.0);
        py.detach(move || pdbiox::InputBuffer::open(path, limits))
            .map(Self)
            .map_err(|finding| read_error(py, &[finding]))
    }

    #[staticmethod]
    #[pyo3(signature = (data, limits=None))]
    fn from_reader(data: &Bound<'_, PyBytes>, limits: Option<&PyLimits>) -> PyResult<Self> {
        let limits = limits.map_or_else(pdbiox::Limits::default, |value| value.0);
        pdbiox::InputBuffer::from_reader(Cursor::new(data.as_bytes()), limits)
            .map(Self)
            .map_err(|finding| read_error(data.py(), &[finding]))
    }

    #[getter]
    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    #[getter]
    fn kind(&self) -> PyInputKind {
        self.0.kind().into()
    }

    #[getter]
    fn origin(&self) -> Option<String> {
        self.0.origin().map(str::to_owned)
    }

    #[getter]
    fn bytes(&self, py: Python<'_>) -> Py<PyBytes> {
        PyBytes::new(py, self.0.as_bytes()).unbind()
    }

    fn with_origin(&self, origin: &str) -> Self {
        Self(self.0.clone().with_origin(origin))
    }
}

#[pymethods]
impl PySelect {
    #[new]
    fn new() -> Self {
        Self
    }

    fn accept_model(&self, model: u32) -> bool {
        let _ = model;
        true
    }

    fn accept_chain(&self, chain: u32) -> bool {
        let _ = chain;
        true
    }

    fn accept_residue(&self, residue: u32) -> bool {
        let _ = residue;
        true
    }

    fn accept_atom(&self, atom: u32) -> bool {
        let _ = atom;
        true
    }
}

#[pymethods]
impl PyReader {
    #[staticmethod]
    fn read(
        py: Python<'_>,
        input: &PyInputBuffer,
        options: &PyReadOptions,
    ) -> PyResult<PyReadReport> {
        let data = input.0.as_bytes().to_vec();
        let options = options.0.clone();
        py.detach(move || pdbiox::read_bytes(data, None, &options))
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
fn write_output(py: Python<'_>, path: PathBuf, data: &Bound<'_, PyBytes>) -> PyResult<()> {
    let bytes = data.as_bytes().to_vec();
    py.detach(move || pdbiox::core::write_output(path, &bytes))
        .map_err(|finding| read_error(py, std::slice::from_ref(&finding)))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyCompression>()?;
    module.getattr("Compression")?.setattr(
        "None",
        module.getattr("Compression")?.getattr("NoCompression")?,
    )?;
    module.add_class::<PyInputKind>()?;
    module.add_class::<PyInputBuffer>()?;
    module.add_class::<PySelect>()?;
    module.add_class::<PyReader>()?;
    module.add_function(wrap_pyfunction!(write_output, module)?)?;
    module.add("SelectAll", module.getattr("Select")?)?;
    module.add("ReadResult", module.getattr("ReadReport")?)?;
    Ok(())
}
