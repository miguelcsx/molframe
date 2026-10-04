//! Substitution matrices, multiple alignment and trees.

use crate::sequence::PyScoring;
use molframe::sequence::{
    self as seq, LadderDirection, MatrixProfile, MsaOptions, SubstitutionMatrix, Tree,
};
use numpy::{PyArray2, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::prelude::*;

/// A residue substitution matrix with the identity of its published source.
#[derive(Clone, Debug)]
#[pyclass(
    name = "SubstitutionMatrix",
    frozen,
    skip_from_py_object,
    module = "molframe.sequence"
)]
pub(crate) struct PySubstitutionMatrix {
    pub(crate) inner: SubstitutionMatrix,
}

impl PySubstitutionMatrix {
    fn load(profile: MatrixProfile) -> PyResult<Self> {
        seq::load_matrix(profile)
            .map(|inner| Self { inner })
            .map_err(crate::error::kernel)
    }
}

#[pymethods]
impl PySubstitutionMatrix {
    /// A BLOSUM matrix: `62` (the usual) or a level from 30 to 90 in steps of 5.
    #[staticmethod]
    fn blosum(level: u8) -> PyResult<Self> {
        if level == 62 {
            return Ok(Self {
                inner: seq::blosum62(),
            });
        }
        Self::load(MatrixProfile::Blosum(level))
    }

    /// A PAM matrix; `distance` is 10 to 500 in steps of 10.
    #[staticmethod]
    fn pam(distance: u16) -> PyResult<Self> {
        Self::load(MatrixProfile::Pam(distance))
    }

    /// The identity matrix: a match scores 1 and a mismatch 0.
    #[staticmethod]
    fn identity() -> PyResult<Self> {
        Self::load(MatrixProfile::Identity)
    }

    /// The NCBI `NUC.4.4` nucleotide matrix.
    #[staticmethod]
    fn nuc44() -> PyResult<Self> {
        Self::load(MatrixProfile::Nuc44)
    }

    /// The score of aligning residue `first` with residue `second`.
    fn score(&self, first: &str, second: &str) -> PyResult<i32> {
        match (first.as_bytes(), second.as_bytes()) {
            ([a], [b]) => Ok(self.inner.get(*a, *b)),
            _ => Err(crate::error::value(
                "score takes two single-letter residues",
            )),
        }
    }

    /// The matrix's name, for example `BLOSUM62`.
    #[getter]
    fn name(&self) -> String {
        self.inner.identity().name().to_owned()
    }

    /// The version of the published table.
    #[getter]
    fn version(&self) -> String {
        self.inner.identity().version().to_owned()
    }

    /// Where the table was retrieved from.
    #[getter]
    fn source(&self) -> String {
        self.inner.identity().source().to_owned()
    }

    fn __repr__(&self) -> String {
        format!("SubstitutionMatrix({})", self.name())
    }
}

/// A progressive multiple alignment of `sequences`, as equal-length rows with
/// `-` for gaps.
///
/// `refinement_passes` is the number of deterministic leave-one-out refinement
/// passes, and `memory_limit` bounds, in bytes, the workspace the alignment may
/// use.
#[pyfunction]
#[pyo3(signature = (sequences, *, scoring, refinement_passes=0, memory_limit=None))]
pub(crate) fn msa(
    py: Python<'_>,
    sequences: Vec<String>,
    scoring: &PyScoring,
    refinement_passes: usize,
    memory_limit: Option<usize>,
) -> PyResult<Vec<String>> {
    let mut options = MsaOptions::progressive(scoring.0).with_refinement_passes(refinement_passes);
    if let Some(bytes) = memory_limit {
        options = options.with_memory_limit(bytes);
    }
    let bytes: Vec<Vec<u8>> = sequences.into_iter().map(String::into_bytes).collect();
    let rows = py
        .detach(|| {
            let views: Vec<&[u8]> = bytes.iter().map(Vec::as_slice).collect();
            seq::progressive_msa(&views, options)
        })
        .map_err(crate::error::kernel)?;
    Ok(rows
        .into_iter()
        .map(|row| String::from_utf8_lossy(&row).into_owned())
        .collect())
}

/// A rooted binary tree with named leaves and branch lengths.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Tree",
    frozen,
    skip_from_py_object,
    module = "molframe.sequence"
)]
pub(crate) struct PyTree {
    inner: Tree,
}

#[pymethods]
impl PyTree {
    /// Reads a Newick string.
    #[staticmethod]
    fn from_newick(text: &str) -> PyResult<Self> {
        Tree::from_newick(text)
            .map(|inner| Self { inner })
            .map_err(crate::error::kernel)
    }

    /// The tree as a Newick string.
    fn newick(&self) -> String {
        self.inner.to_newick()
    }

    /// The number of leaves.
    #[getter]
    fn leaf_count(&self) -> usize {
        self.inner.leaf_count()
    }

    /// The leaf names, in the order a depth-first walk reaches them.
    fn leaves(&self) -> Vec<String> {
        self.inner
            .traverse(seq::TraversalOrder::Preorder)
            .into_iter()
            .filter_map(|node| match node {
                Tree::Leaf { name } => Some(name.clone()),
                Tree::Clade { .. } => None,
            })
            .collect()
    }

    /// A copy with every clade's children put in a canonical order: smaller
    /// clades first (`"ascending"`) or larger first (`"descending"`).
    fn ladderize(&self, direction: &str) -> PyResult<Self> {
        let direction = match direction {
            "ascending" => LadderDirection::Ascending,
            "descending" => LadderDirection::Descending,
            _ => {
                return Err(crate::error::value(
                    "direction must be 'ascending' or 'descending'",
                ));
            }
        };
        let mut copy = self.inner.clone();
        copy.ladderize(direction);
        Ok(Self { inner: copy })
    }

    /// A copy rerooted on the branch above `leaf`, `fraction` of the way up
    /// from the leaf (`0` at the leaf, `1` at its parent).
    fn reroot_at_leaf(&self, leaf: &str, fraction: f64) -> PyResult<Self> {
        self.inner
            .reroot_at_leaf(leaf, fraction)
            .map(|inner| Self { inner })
            .map_err(crate::error::kernel)
    }

    fn __repr__(&self) -> String {
        format!("Tree(leaves={})", self.inner.leaf_count())
    }
}

fn rows_of(distances: &Bound<'_, PyArray2<f64>>) -> PyResult<Vec<Vec<f64>>> {
    let view = distances.readonly();
    let shape = view.shape();
    let flat = view.as_slice().map_err(|_| {
        crate::error::value("distances must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    Ok(flat.chunks(shape[1].max(1)).map(<[f64]>::to_vec).collect())
}

fn build(
    py: Python<'_>,
    labels: &[String],
    distances: &Bound<'_, PyArray2<f64>>,
    builder: fn(&[&str], &[Vec<f64>]) -> Option<Tree>,
) -> PyResult<PyTree> {
    let rows = rows_of(distances)?;
    let names: Vec<&str> = labels.iter().map(String::as_str).collect();
    py.detach(|| builder(&names, &rows))
        .map(|inner| PyTree { inner })
        .ok_or_else(|| {
            crate::error::value("distances must be a square matrix with one row per label")
        })
}

/// Builds a neighbour-joining tree from taxon labels and a symmetric distance matrix.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn neighbor_joining(
    py: Python<'_>,
    labels: Vec<String>,
    distances: &Bound<'_, PyArray2<f64>>,
) -> PyResult<PyTree> {
    build(py, &labels, distances, seq::neighbor_joining)
}

/// Builds a UPGMA tree from taxon labels and a symmetric distance matrix.
#[pyfunction]
#[allow(clippy::needless_pass_by_value)]
pub(crate) fn upgma(
    py: Python<'_>,
    labels: Vec<String>,
    distances: &Bound<'_, PyArray2<f64>>,
) -> PyResult<PyTree> {
    build(py, &labels, distances, seq::upgma)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySubstitutionMatrix>()?;
    module.add_class::<PyTree>()?;
    module.add_function(wrap_pyfunction!(msa, module)?)?;
    module.add_function(wrap_pyfunction!(neighbor_joining, module)?)?;
    module.add_function(wrap_pyfunction!(upgma, module)?)
}
