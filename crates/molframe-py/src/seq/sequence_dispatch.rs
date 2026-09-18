//! Common typed sequence-format dispatch delegated to `molframe-seq`.

use super::sequence_errors::SequenceFormatError;
use super::{PyFastaRecord, PyFastqRecord, PyTree};
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;

#[pyclass(name = "SequenceFormat", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySequenceFormat {
    Fasta,
    Fastq,
    Stockholm,
    Clustal,
    Phylip,
    A2m,
    A3m,
    Newick,
}

#[pyclass(name = "SequenceDocumentKind", frozen, eq, eq_int, skip_from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PySequenceDocumentKind {
    Records,
    Fastq,
    Tree,
}

#[pyclass(name = "SequenceDocument", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySequenceDocument(molframe::seq::SequenceDocument);

#[pymethods]
impl PySequenceDocument {
    #[staticmethod]
    fn records(records: Vec<PyFastaRecord>) -> Self {
        Self(molframe::seq::SequenceDocument::Records(
            records.into_iter().map(Into::into).collect(),
        ))
    }

    #[staticmethod]
    fn fastq(records: Vec<PyFastqRecord>) -> Self {
        Self(molframe::seq::SequenceDocument::Fastq(
            records.into_iter().map(Into::into).collect(),
        ))
    }

    #[staticmethod]
    fn tree(tree: &PyTree) -> Self {
        Self(molframe::seq::SequenceDocument::Tree(tree.0.clone()))
    }

    #[getter]
    const fn kind(&self) -> PySequenceDocumentKind {
        match &self.0 {
            molframe::seq::SequenceDocument::Records(_) => PySequenceDocumentKind::Records,
            molframe::seq::SequenceDocument::Fastq(_) => PySequenceDocumentKind::Fastq,
            molframe::seq::SequenceDocument::Tree(_) => PySequenceDocumentKind::Tree,
        }
    }

    fn record_values(&self) -> PyResult<Vec<PyFastaRecord>> {
        match &self.0 {
            molframe::seq::SequenceDocument::Records(records) => {
                Ok(records.iter().cloned().map(Into::into).collect())
            }
            _ => Err(wrong_document("records")),
        }
    }

    fn fastq_values(&self) -> PyResult<Vec<PyFastqRecord>> {
        match &self.0 {
            molframe::seq::SequenceDocument::Fastq(records) => {
                Ok(records.iter().cloned().map(Into::into).collect())
            }
            _ => Err(wrong_document("FASTQ records")),
        }
    }

    fn tree_value(&self) -> PyResult<PyTree> {
        match &self.0 {
            molframe::seq::SequenceDocument::Tree(tree) => Ok(PyTree(tree.clone())),
            _ => Err(wrong_document("tree")),
        }
    }
}

#[pyfunction]
pub(crate) fn read_sequence(
    py: Python<'_>,
    text: String,
    format: PySequenceFormat,
) -> PyResult<PySequenceDocument> {
    py.detach(move || molframe::seq::read_sequence(&text, format.into()))
        .map(PySequenceDocument)
        .map_err(|error| SequenceFormatError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_sequence(
    py: Python<'_>,
    document: &PySequenceDocument,
    format: PySequenceFormat,
) -> PyResult<String> {
    let document = document.0.clone();
    py.detach(move || molframe::seq::write_sequence(&document, format.into()))
        .map_err(|error| SequenceFormatError::new_err(error.to_string()))
}

impl From<PySequenceFormat> for molframe::seq::SequenceFormat {
    fn from(value: PySequenceFormat) -> Self {
        match value {
            PySequenceFormat::Fasta => Self::Fasta,
            PySequenceFormat::Fastq => Self::Fastq,
            PySequenceFormat::Stockholm => Self::Stockholm,
            PySequenceFormat::Clustal => Self::Clustal,
            PySequenceFormat::Phylip => Self::Phylip,
            PySequenceFormat::A2m => Self::A2m,
            PySequenceFormat::A3m => Self::A3m,
            PySequenceFormat::Newick => Self::Newick,
        }
    }
}

fn wrong_document(expected: &str) -> PyErr {
    PyTypeError::new_err(format!("sequence document does not contain {expected}"))
}
