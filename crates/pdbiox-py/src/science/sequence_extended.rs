//! Deterministic sequence collections, formats and phylogenetic kernels.

use super::sequence::{PyAlignment, PyScoring};
use numpy::PyReadonlyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[pyclass(name = "MsaOptions", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyMsaOptions(pdbiox::seq::MsaOptions);

#[pymethods]
impl PyMsaOptions {
    #[new]
    #[pyo3(signature = (scoring, *, refinement_passes=0))]
    fn new(scoring: PyScoring, refinement_passes: usize) -> Self {
        Self(
            pdbiox::seq::MsaOptions::progressive(scoring.inner())
                .with_refinement_passes(refinement_passes),
        )
    }
}

#[pyclass(name = "FastaRecord", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFastaRecord {
    #[pyo3(get)]
    id: String,
    #[pyo3(get)]
    description: String,
    sequence: Vec<u8>,
}

#[pymethods]
impl PyFastaRecord {
    #[new]
    fn new(id: String, description: String, sequence: Vec<u8>) -> Self {
        Self {
            id,
            description,
            sequence,
        }
    }

    #[getter]
    fn sequence<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.sequence)
    }
}

macro_rules! parser {
    ($name:ident, $kernel:path) => {
        #[pyfunction]
        pub(crate) fn $name(py: Python<'_>, text: String) -> Vec<PyFastaRecord> {
            py.detach(|| $kernel(&text))
                .into_iter()
                .map(Into::into)
                .collect()
        }
    };
}

macro_rules! writer {
    ($name:ident, $kernel:path) => {
        #[pyfunction]
        pub(crate) fn $name(py: Python<'_>, records: Vec<PyFastaRecord>) -> String {
            let records: Vec<_> = records.into_iter().map(Into::into).collect();
            py.detach(|| $kernel(&records))
        }
    };
}

parser!(parse_fasta, pdbiox::seq::parse_fasta);
parser!(parse_clustal, pdbiox::seq::parse_clustal);
parser!(parse_phylip, pdbiox::seq::parse_phylip);
parser!(parse_stockholm, pdbiox::seq::parse_stockholm);
writer!(write_fasta, pdbiox::seq::write_fasta);
writer!(write_clustal, pdbiox::seq::write_clustal);
writer!(write_phylip, pdbiox::seq::write_phylip);
writer!(write_stockholm, pdbiox::seq::write_stockholm);

#[pyfunction]
pub(crate) fn parse_a3m(py: Python<'_>, text: String) -> PyResult<Vec<PyFastaRecord>> {
    py.detach(move || pdbiox::seq::parse_a3m(&text))
        .map(|records| records.into_iter().map(Into::into).collect())
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_a3m(py: Python<'_>, records: Vec<PyFastaRecord>) -> PyResult<String> {
    let records = records.into_iter().map(Into::into).collect::<Vec<_>>();
    py.detach(|| pdbiox::seq::write_a3m(&records))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn a3m_match_columns(py: Python<'_>, sequence: Vec<u8>) -> Vec<u8> {
    let sequence = sequence.into_boxed_slice();
    py.detach(|| pdbiox::seq::a3m_match_columns(&sequence))
}

#[pyfunction]
pub(crate) fn a2m_match_columns(py: Python<'_>, sequence: Vec<u8>) -> Vec<u8> {
    let sequence = sequence.into_boxed_slice();
    py.detach(|| pdbiox::seq::a2m_match_columns(&sequence))
}

#[pyfunction]
pub(crate) fn kmer_counts(py: Python<'_>, sequence: Vec<u8>, k: usize) -> Vec<(Vec<u8>, u32)> {
    let sequence = sequence.into_boxed_slice();
    py.detach(|| pdbiox::seq::kmer_counts(&sequence, k))
}

#[pyfunction]
pub(crate) fn minimizers(
    py: Python<'_>,
    sequence: Vec<u8>,
    k: usize,
    window: usize,
) -> Vec<(usize, Vec<u8>)> {
    let sequence = sequence.into_boxed_slice();
    py.detach(|| pdbiox::seq::minimizers(&sequence, k, window))
}

#[pyfunction]
pub(crate) fn multiple_sequence_alignment(
    py: Python<'_>,
    sequences: Vec<Vec<u8>>,
    scoring: PyScoring,
) -> PyResult<Vec<Vec<u8>>> {
    let sequences = sequences.into_boxed_slice();
    py.detach(|| {
        let borrowed: Vec<&[u8]> = sequences.iter().map(Vec::as_slice).collect();
        pdbiox::seq::msa(&borrowed, scoring.inner())
    })
    .map_err(msa_error)
}

#[pyfunction]
pub(crate) fn progressive_msa(
    py: Python<'_>,
    sequences: Vec<Vec<u8>>,
    options: PyMsaOptions,
) -> PyResult<Vec<Vec<u8>>> {
    let sequences = sequences.into_boxed_slice();
    py.detach(|| {
        let borrowed: Vec<&[u8]> = sequences.iter().map(Vec::as_slice).collect();
        pdbiox::seq::progressive_msa(&borrowed, options.0)
    })
    .map_err(msa_error)
}

fn msa_error(error: pdbiox::seq::MsaError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

#[pyfunction]
pub(crate) fn align_global_banded(
    py: Python<'_>,
    left: Vec<u8>,
    right: Vec<u8>,
    scoring: PyScoring,
    band: usize,
) -> PyResult<PyAlignment> {
    let left = left.into_boxed_slice();
    let right = right.into_boxed_slice();
    py.detach(move || pdbiox::seq::global_banded(&left, &right, scoring.inner(), band))
        .map(Into::into)
        .map_err(super::sequence::align_error)
}

#[pyfunction]
pub(crate) fn seed_and_extend(
    py: Python<'_>,
    query: Vec<u8>,
    reference: Vec<u8>,
    k: usize,
    scoring: PyScoring,
) -> PyResult<Option<PyAlignment>> {
    let query = query.into_boxed_slice();
    let reference = reference.into_boxed_slice();
    py.detach(move || pdbiox::seq::seed_and_extend(&query, &reference, k, scoring.inner()))
        .map(|alignment| alignment.map(Into::into))
        .map_err(super::sequence::align_error)
}

#[pyclass(name = "SubstitutionMatrix", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySubstitutionMatrix(pub(crate) pdbiox::seq::SubstitutionMatrix);

#[pymethods]
impl PySubstitutionMatrix {
    fn score(&self, left: u8, right: u8) -> i32 {
        self.0.get(left, right)
    }

    #[getter]
    fn name(&self) -> &str {
        self.0.identity().name()
    }

    #[getter]
    fn version(&self) -> &str {
        self.0.identity().version()
    }

    #[getter]
    fn source(&self) -> &str {
        self.0.identity().source()
    }

    #[getter]
    fn retrieved(&self) -> &str {
        self.0.identity().retrieved()
    }
}

#[pyfunction]
pub(crate) fn blosum62() -> PySubstitutionMatrix {
    PySubstitutionMatrix(pdbiox::seq::blosum62())
}

macro_rules! matrix_aligner {
    ($name:ident, $kernel:path) => {
        #[pyfunction]
        pub(crate) fn $name(
            py: Python<'_>,
            left: Vec<u8>,
            right: Vec<u8>,
            matrix: &PySubstitutionMatrix,
            gap_open: i32,
            gap_extend: i32,
        ) -> PyResult<PyAlignment> {
            let matrix = matrix.0.clone();
            py.detach(move || $kernel(&left, &right, &matrix, gap_open, gap_extend))
                .map(Into::into)
                .map_err(super::sequence::align_error)
        }
    };
}

matrix_aligner!(align_global_matrix, pdbiox::seq::global_matrix);
matrix_aligner!(align_local_matrix, pdbiox::seq::local_matrix);
matrix_aligner!(align_semi_global_matrix, pdbiox::seq::semi_global_matrix);

#[pyclass(name = "Tree", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTree(pub(crate) pdbiox::seq::Tree);

#[pyclass(name = "LadderDirection", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyLadderDirection {
    Ascending,
    Descending,
}

#[pymethods]
impl PyTree {
    #[staticmethod]
    fn from_newick(text: &str) -> PyResult<Self> {
        pdbiox::seq::Tree::from_newick(text)
            .map(Self)
            .map_err(|error| PyValueError::new_err(format!("invalid Newick tree: {error:?}")))
    }

    fn to_newick(&self) -> String {
        self.0.to_newick()
    }

    #[getter]
    fn leaf_count(&self) -> usize {
        self.0.leaf_count()
    }

    fn ladderized(&self, direction: PyLadderDirection) -> Self {
        let mut tree = self.0.clone();
        tree.ladderize(direction.into());
        Self(tree)
    }

    fn rerooted_at_leaf(&self, leaf: &str, fraction_from_leaf: f64) -> PyResult<Self> {
        self.0
            .reroot_at_leaf(leaf, fraction_from_leaf)
            .map(Self)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}

impl From<PyLadderDirection> for pdbiox::seq::LadderDirection {
    fn from(value: PyLadderDirection) -> Self {
        match value {
            PyLadderDirection::Ascending => Self::Ascending,
            PyLadderDirection::Descending => Self::Descending,
        }
    }
}

#[pyfunction]
pub(crate) fn neighbor_joining(
    py: Python<'_>,
    labels: Vec<String>,
    distances: PyReadonlyArray2<'_, f64>,
) -> PyResult<Option<PyTree>> {
    tree_from_distances(py, labels, distances, pdbiox::seq::neighbor_joining)
}

#[pyfunction]
pub(crate) fn upgma(
    py: Python<'_>,
    labels: Vec<String>,
    distances: PyReadonlyArray2<'_, f64>,
) -> PyResult<Option<PyTree>> {
    tree_from_distances(py, labels, distances, pdbiox::seq::upgma)
}

fn tree_from_distances(
    py: Python<'_>,
    labels: Vec<String>,
    distances: PyReadonlyArray2<'_, f64>,
    kernel: fn(&[&str], &[Vec<f64>]) -> Option<pdbiox::seq::Tree>,
) -> PyResult<Option<PyTree>> {
    let labels = labels.into_boxed_slice();
    let view = distances.as_array();
    if view.nrows() != view.ncols() || view.nrows() != labels.len() {
        return Err(PyValueError::new_err(
            "distance matrix must be square and match labels",
        ));
    }
    let matrix: Vec<Vec<f64>> = view.rows().into_iter().map(|row| row.to_vec()).collect();
    drop(distances);
    Ok(py
        .detach(|| {
            let names: Vec<&str> = labels.iter().map(String::as_str).collect();
            kernel(&names, &matrix)
        })
        .map(PyTree))
}

impl From<pdbiox::seq::FastaRecord> for PyFastaRecord {
    fn from(value: pdbiox::seq::FastaRecord) -> Self {
        Self {
            id: value.id,
            description: value.description,
            sequence: value.sequence,
        }
    }
}

impl From<PyFastaRecord> for pdbiox::seq::FastaRecord {
    fn from(value: PyFastaRecord) -> Self {
        Self {
            id: value.id,
            description: value.description,
            sequence: value.sequence,
        }
    }
}
