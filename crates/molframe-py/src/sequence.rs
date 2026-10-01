//! Mechanical adapters for pairwise alignment, FASTA and k-mer counting.

use molframe::sequence::{self as seq, Alignment, FastaRecord, Scoring};
use pyo3::{exceptions::PyValueError, prelude::*};

/// One optimal pairwise alignment.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Alignment",
    frozen,
    skip_from_py_object,
    module = "molframe.sequence"
)]
struct PyAlignment {
    inner: Alignment,
    left: Vec<u8>,
    right: Vec<u8>,
}

#[pymethods]
impl PyAlignment {
    #[getter]
    fn score(&self) -> i32 {
        self.inner.score
    }

    /// Aligned index pairs; `None` marks a gap in that sequence.
    #[getter]
    fn columns(&self) -> Vec<(Option<usize>, Option<usize>)> {
        self.inner
            .columns
            .iter()
            .map(|column| (column.left, column.right))
            .collect()
    }

    /// Both sequences laid out over the alignment columns, gaps as `-`.
    #[getter]
    fn aligned(&self) -> (String, String) {
        let residue = |sequence: &[u8], index: Option<usize>| {
            index
                .and_then(|index| sequence.get(index))
                .map_or('-', |byte| char::from(*byte))
        };
        let mut left = String::with_capacity(self.inner.columns.len());
        let mut right = String::with_capacity(self.inner.columns.len());
        for column in &self.inner.columns {
            left.push(residue(&self.left, column.left));
            right.push(residue(&self.right, column.right));
        }
        (left, right)
    }

    fn __len__(&self) -> usize {
        self.inner.columns.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Alignment(score={}, columns={})",
            self.inner.score,
            self.inner.columns.len()
        )
    }
}

/// One FASTA record.
#[derive(Clone, Debug)]
#[pyclass(
    name = "FastaRecord",
    frozen,
    skip_from_py_object,
    module = "molframe.sequence"
)]
struct PyFastaRecord {
    inner: FastaRecord,
}

#[pymethods]
impl PyFastaRecord {
    #[new]
    #[pyo3(signature = (id, sequence, description=String::new()))]
    fn new(id: String, sequence: &str, description: String) -> Self {
        Self {
            inner: FastaRecord {
                id,
                description,
                sequence: sequence.as_bytes().to_vec(),
            },
        }
    }

    #[getter]
    fn id(&self) -> &str {
        &self.inner.id
    }

    #[getter]
    fn description(&self) -> &str {
        &self.inner.description
    }

    #[getter]
    fn sequence(&self) -> String {
        String::from_utf8_lossy(&self.inner.sequence).into_owned()
    }

    fn __len__(&self) -> usize {
        self.inner.sequence.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "FastaRecord(id={:?}, length={})",
            self.inner.id,
            self.inner.sequence.len()
        )
    }
}

/// Aligns two sequences globally, locally or semi-globally.
#[pyfunction]
#[pyo3(signature = (left, right, *, mode="global", match_score=1, mismatch_score=-1, gap_open=-2, gap_extend=-1))]
fn align(
    py: Python<'_>,
    left: &str,
    right: &str,
    mode: &str,
    match_score: i32,
    mismatch_score: i32,
    gap_open: i32,
    gap_extend: i32,
) -> PyResult<PyAlignment> {
    let scoring = Scoring {
        match_score,
        mismatch_score,
        gap_open,
        gap_extend,
    };
    if gap_open > 0 || gap_extend > 0 {
        return Err(PyValueError::new_err(
            "gap_open and gap_extend must be zero or negative",
        ));
    }
    let solve = match mode {
        "global" => seq::global,
        "local" => seq::local,
        "semi_global" => seq::semi_global,
        _ => {
            return Err(PyValueError::new_err(
                "mode must be 'global', 'local' or 'semi_global'",
            ));
        }
    };
    let (left_bytes, right_bytes) = (left.as_bytes().to_vec(), right.as_bytes().to_vec());
    let inner = py
        .detach(|| solve(&left_bytes, &right_bytes, scoring))
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    Ok(PyAlignment {
        inner,
        left: left_bytes,
        right: right_bytes,
    })
}

/// Parses FASTA text into records, preserving order.
#[pyfunction]
fn parse_fasta(text: &str) -> Vec<PyFastaRecord> {
    seq::parse_fasta(text)
        .into_iter()
        .map(|inner| PyFastaRecord { inner })
        .collect()
}

/// Writes records as FASTA text.
#[pyfunction]
fn write_fasta(records: Vec<PyRef<'_, PyFastaRecord>>) -> String {
    let records: Vec<FastaRecord> = records.iter().map(|record| record.inner.clone()).collect();
    seq::write_fasta(&records)
}

/// Counts each distinct k-mer of a sequence.
#[pyfunction]
fn kmer_counts(sequence: &str, k: usize) -> PyResult<Vec<(String, u32)>> {
    if k == 0 {
        return Err(PyValueError::new_err("k must be positive"));
    }
    Ok(seq::kmer_counts(sequence.as_bytes(), k)
        .into_iter()
        .map(|(kmer, count)| (String::from_utf8_lossy(&kmer).into_owned(), count))
        .collect())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAlignment>()?;
    module.add_class::<PyFastaRecord>()?;
    module.add_function(wrap_pyfunction!(align, module)?)?;
    module.add_function(wrap_pyfunction!(parse_fasta, module)?)?;
    module.add_function(wrap_pyfunction!(write_fasta, module)?)?;
    module.add_function(wrap_pyfunction!(kmer_counts, module)?)
}
