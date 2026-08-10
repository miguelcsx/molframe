//! Strict sequence formats and bounded sequence-search policies.

use super::{PyAlignment, PyScoring, PySubstitutionMatrix};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[pyclass(name = "FastqRecord", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFastqRecord {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    description: String,
    sequence: Vec<u8>,
    quality: Vec<u8>,
}

#[pymethods]
impl PyFastqRecord {
    #[new]
    fn new(id: String, description: String, sequence: Vec<u8>, quality: Vec<u8>) -> Self {
        Self {
            id,
            description,
            sequence,
            quality,
        }
    }

    #[getter]
    fn sequence<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.sequence)
    }

    #[getter]
    fn quality<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.quality)
    }
}

#[pyfunction]
pub(crate) fn parse_fastq(py: Python<'_>, text: String) -> PyResult<Vec<PyFastqRecord>> {
    py.detach(move || pdbiox::seq::parse_fastq(&text))
        .map(|records| records.into_iter().map(Into::into).collect())
        .map_err(value_error)
}

#[pyfunction]
pub(crate) fn write_fastq(py: Python<'_>, records: Vec<PyFastqRecord>) -> PyResult<String> {
    let records = records.into_iter().map(Into::into).collect::<Vec<_>>();
    py.detach(|| pdbiox::seq::write_fastq(&records))
        .map_err(value_error)
}

macro_rules! strict_fasta_format {
    ($parse:ident, $write:ident, $rust_parse:path, $rust_write:path) => {
        #[pyfunction]
        pub(crate) fn $parse(py: Python<'_>, text: String) -> PyResult<Vec<super::PyFastaRecord>> {
            py.detach(|| $rust_parse(&text))
                .map(|records| records.into_iter().map(Into::into).collect())
                .map_err(value_error)
        }

        #[pyfunction]
        pub(crate) fn $write(
            py: Python<'_>,
            records: Vec<super::PyFastaRecord>,
        ) -> PyResult<String> {
            let records = records.into_iter().map(Into::into).collect::<Vec<_>>();
            py.detach(|| $rust_write(&records)).map_err(value_error)
        }
    };
}

strict_fasta_format!(
    parse_a2m,
    write_a2m,
    pdbiox::seq::parse_a2m,
    pdbiox::seq::write_a2m
);

#[pyclass(name = "MatrixProfile", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMatrixProfile(pdbiox::seq::MatrixProfile);

#[pymethods]
impl PyMatrixProfile {
    #[staticmethod]
    fn blosum(level: u8) -> Self {
        Self(pdbiox::seq::MatrixProfile::Blosum(level))
    }

    #[staticmethod]
    fn pam(distance: u16) -> Self {
        Self(pdbiox::seq::MatrixProfile::Pam(distance))
    }

    #[staticmethod]
    fn identity() -> Self {
        Self(pdbiox::seq::MatrixProfile::Identity)
    }

    #[staticmethod]
    fn nuc44() -> Self {
        Self(pdbiox::seq::MatrixProfile::Nuc44)
    }
}

#[pyfunction]
pub(crate) fn load_matrix(
    py: Python<'_>,
    profile: PyMatrixProfile,
) -> PyResult<PySubstitutionMatrix> {
    py.detach(|| pdbiox::seq::load_matrix(profile.0))
        .map(PySubstitutionMatrix)
        .map_err(value_error)
}

#[pyclass(name = "SimilarKmerOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PySimilarKmerOptions(pdbiox::seq::SimilarKmerOptions);

#[pymethods]
impl PySimilarKmerOptions {
    #[new]
    fn new(minimum_score: i32, max_results: usize, max_candidates: usize) -> Self {
        Self(pdbiox::seq::SimilarKmerOptions {
            minimum_score,
            max_results,
            max_candidates,
        })
    }
}

#[pyclass(name = "SimilarKmer", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySimilarKmer {
    word: Vec<u8>,
    #[pyo3(get)]
    score: i32,
}

#[pymethods]
impl PySimilarKmer {
    #[getter]
    fn word<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.word)
    }
}

#[pyfunction]
pub(crate) fn similar_kmers(
    py: Python<'_>,
    query: Vec<u8>,
    alphabet: Vec<u8>,
    matrix: &PySubstitutionMatrix,
    options: PySimilarKmerOptions,
) -> PyResult<Vec<PySimilarKmer>> {
    let matrix = matrix.0.clone();
    py.detach(move || pdbiox::seq::similar_kmers(&query, &alphabet, &matrix, options.0))
        .map(|values| {
            values
                .into_iter()
                .map(|value| PySimilarKmer {
                    word: value.word,
                    score: value.score,
                })
                .collect()
        })
        .map_err(value_error)
}

#[pyclass(name = "AlignmentMode", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyAlignmentMode {
    Global,
    Local,
    SemiGlobal,
}

#[pyclass(name = "RegionOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyRegionOptions(pdbiox::seq::RegionOptions);

#[pymethods]
impl PyRegionOptions {
    #[new]
    #[pyo3(signature = (left, right, mode, *, band=None))]
    fn new(
        left: (usize, usize),
        right: (usize, usize),
        mode: PyAlignmentMode,
        band: Option<usize>,
    ) -> Self {
        Self(pdbiox::seq::RegionOptions {
            left: left.0..left.1,
            right: right.0..right.1,
            mode: mode.into(),
            band,
        })
    }
}

#[pyfunction]
pub(crate) fn align_region(
    py: Python<'_>,
    left: Vec<u8>,
    right: Vec<u8>,
    scoring: PyScoring,
    options: &PyRegionOptions,
) -> PyResult<PyAlignment> {
    let options = options.0.clone();
    py.detach(move || pdbiox::seq::align_region(&left, &right, scoring.inner(), &options))
        .map(Into::into)
        .map_err(value_error)
}

impl From<PyAlignmentMode> for pdbiox::seq::AlignmentMode {
    fn from(value: PyAlignmentMode) -> Self {
        match value {
            PyAlignmentMode::Global => Self::Global,
            PyAlignmentMode::Local => Self::Local,
            PyAlignmentMode::SemiGlobal => Self::SemiGlobal,
        }
    }
}

impl From<pdbiox::seq::FastqRecord> for PyFastqRecord {
    fn from(value: pdbiox::seq::FastqRecord) -> Self {
        Self {
            id: value.id,
            description: value.description,
            sequence: value.sequence,
            quality: value.quality,
        }
    }
}

impl From<PyFastqRecord> for pdbiox::seq::FastqRecord {
    fn from(value: PyFastqRecord) -> Self {
        Self {
            id: value.id,
            description: value.description,
            sequence: value.sequence,
            quality: value.quality,
        }
    }
}

fn value_error(error: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(error.to_string())
}
