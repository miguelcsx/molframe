//! Pairwise affine-gap alignment with one native call per alignment.

use super::sequence_errors::AlignError;
use pyo3::prelude::*;

#[pyclass(name = "Scoring", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyScoring(pdbiox::seq::Scoring);

impl PyScoring {
    pub(crate) const fn inner(self) -> pdbiox::seq::Scoring {
        self.0
    }
}

#[pyclass(name = "AlignmentColumn", frozen, from_py_object)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PyColumn {
    #[pyo3(get)]
    left: Option<usize>,
    #[pyo3(get)]
    right: Option<usize>,
}

#[pymethods]
impl PyScoring {
    #[new]
    #[pyo3(signature = (match_score=1, mismatch_score=-1, gap_open=-2, gap_extend=-1))]
    fn new(match_score: i32, mismatch_score: i32, gap_open: i32, gap_extend: i32) -> Self {
        Self(pdbiox::seq::Scoring {
            match_score,
            mismatch_score,
            gap_open,
            gap_extend,
        })
    }

    #[staticmethod]
    fn simple() -> Self {
        Self(pdbiox::seq::Scoring::simple())
    }

    #[getter]
    fn match_score(&self) -> i32 {
        self.0.match_score
    }

    #[getter]
    fn mismatch_score(&self) -> i32 {
        self.0.mismatch_score
    }

    #[getter]
    fn gap_open(&self) -> i32 {
        self.0.gap_open
    }

    #[getter]
    fn gap_extend(&self) -> i32 {
        self.0.gap_extend
    }

    fn substitution(&self, left: u8, right: u8) -> i32 {
        self.0.substitution(left, right)
    }
}

#[pyclass(name = "Alignment", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAlignment {
    score: i32,
    columns: Vec<PyColumn>,
}

#[pymethods]
impl PyAlignment {
    #[getter]
    const fn score(&self) -> i32 {
        self.score
    }
    #[getter]
    fn columns(&self) -> Vec<PyColumn> {
        self.columns.clone()
    }
}

macro_rules! aligner {
    ($name:ident, $python:literal, $kernel:path) => {
        #[pyfunction(name = $python)]
        #[pyo3(signature = (left, right, scoring=None))]
        pub(crate) fn $name(
            py: Python<'_>,
            left: Vec<u8>,
            right: Vec<u8>,
            scoring: Option<PyScoring>,
        ) -> PyResult<PyAlignment> {
            let scoring = scoring.map_or_else(pdbiox::seq::Scoring::simple, |value| value.0);
            py.detach(move || $kernel(&left, &right, scoring))
                .map(Into::into)
                .map_err(align_error)
        }
    };
}

aligner!(global, "align_global", pdbiox::seq::global);
aligner!(local, "align_local", pdbiox::seq::local);
aligner!(semi_global, "align_semi_global", pdbiox::seq::semi_global);

pub(super) fn align_error(error: pdbiox::seq::AlignError) -> PyErr {
    AlignError::new_err(error.to_string())
}

impl From<pdbiox::seq::Alignment> for PyAlignment {
    fn from(value: pdbiox::seq::Alignment) -> Self {
        Self {
            score: value.score,
            columns: value
                .columns
                .into_iter()
                .map(|column| PyColumn {
                    left: column.left,
                    right: column.right,
                })
                .collect(),
        }
    }
}
