//! The governed result envelope: a value, how complete it is, and why.

use pyo3::prelude::*;

/// How much of the intended data an analysis used.
#[derive(Clone, Copy, Debug)]
#[pyclass(name = "Coverage", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyCoverage {
    intended: u32,
    used: u32,
    missing: u32,
    ambiguous: u32,
}

#[pymethods]
impl PyCoverage {
    #[getter]
    const fn intended(&self) -> u32 {
        self.intended
    }

    #[getter]
    const fn used(&self) -> u32 {
        self.used
    }

    #[getter]
    const fn missing(&self) -> u32 {
        self.missing
    }

    #[getter]
    const fn ambiguous(&self) -> u32 {
        self.ambiguous
    }

    /// Used over intended; `1.0` when nothing was intended.
    #[getter]
    fn fraction(&self) -> f64 {
        if self.intended == 0 {
            1.0
        } else {
            f64::from(self.used) / f64::from(self.intended)
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Coverage(used={}, intended={}, missing={}, ambiguous={})",
            self.used, self.intended, self.missing, self.ambiguous
        )
    }
}

/// A result carrying its status, coverage, decisions and provenance.
#[derive(Debug)]
#[pyclass(name = "Analysis", frozen, skip_from_py_object, module = "molframe")]
pub(crate) struct PyAnalysis {
    /// `None` exactly when the analysis is indeterminate: there is no value to hold.
    value: Option<Py<PyAny>>,
    reason: Option<String>,
    reads: Option<Vec<molframe::PolicyField>>,
    estimand: Option<String>,
    origin: Option<Vec<u32>>,
    status: &'static str,
    coverage: PyCoverage,
    warnings: Vec<String>,
    assumptions: Vec<String>,
    provenance: String,
    record: molframe::Provenance,
    profile: Option<&'static str>,
}

impl PyAnalysis {
    /// Wraps a native analysis whose value, if it has one, is already a Python object.
    pub(crate) fn new<T>(analysis: &molframe::Analysis<T>, value: Option<Py<PyAny>>) -> Self {
        let coverage = analysis.coverage;
        Self {
            value,
            reason: analysis.indeterminacy().map(ToString::to_string),
            reads: analysis.provenance.policy_reads(),
            estimand: analysis.provenance.estimand().map(str::to_owned),
            origin: analysis.atom_origin().map(<[u32]>::to_vec),
            status: match analysis.status() {
                molframe::Status::Complete => "complete",
                molframe::Status::Partial => "partial",
                molframe::Status::Ambiguous => "ambiguous",
                molframe::Status::Indeterminate => "indeterminate",
                _ => "unknown",
            },
            coverage: PyCoverage {
                intended: coverage.intended,
                used: coverage.used,
                missing: coverage.missing,
                ambiguous: coverage.ambiguous,
            },
            warnings: analysis.warnings.iter().map(ToString::to_string).collect(),
            assumptions: analysis
                .assumptions
                .iter()
                .map(|assumption| {
                    format!(
                        "{:?} = {} ({:?}, {:?} impact)",
                        assumption.field, assumption.value, assumption.source, assumption.impact
                    )
                })
                .collect(),
            provenance: analysis.provenance.to_json(),
            record: analysis.provenance.clone(),
            profile: analysis.provenance.profile.map(molframe::ProfileId::as_str),
        }
    }
}

impl PyAnalysis {
    /// The policy fields the analysis recorded as applied.
    pub(crate) fn reads(&self) -> Option<Vec<molframe::PolicyField>> {
        self.reads.clone()
    }

    /// The record of what produced the analysis.
    pub(crate) fn record(&self) -> molframe::Provenance {
        self.record.clone()
    }

    /// The value, when there is one, without raising.
    pub(crate) fn answer(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.value.as_ref().map(|value| value.clone_ref(py))
    }

    /// The same envelope around a converted value.
    pub(crate) fn with_value(mut self, value: Option<Py<PyAny>>) -> Self {
        self.value = value;
        self
    }
}

#[pymethods]
impl PyAnalysis {
    /// The answer.
    ///
    /// Raises `IndeterminateError` when `status` is `"indeterminate"`: there is
    /// no defensible answer under the policy, so there is nothing to return.
    #[getter]
    fn value(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        match &self.value {
            Some(value) => Ok(value.clone_ref(py)),
            None => Err(crate::error::indeterminate(match self.reason.as_deref() {
                Some(reason) => reason,
                None => "no defensible answer",
            })),
        }
    }

    /// True when there is an answer; false when the policy leaves none.
    #[getter]
    const fn is_determinate(&self) -> bool {
        self.value.is_some()
    }

    /// For each atom the analysis ran over, the input atom it is.
    ///
    /// The atom indices in a result number the system the policy built: the input with
    /// some atoms dropped (a conformer not chosen, hydrogens excluded) or copied (an
    /// assembly). `None` when the result's atoms are the input's own.
    #[getter]
    fn atom_origin<'py>(&self, py: Python<'py>) -> Option<Bound<'py, numpy::PyArray1<u32>>> {
        use numpy::{PyArrayMethods, ToPyArray};
        let array = self.origin.as_ref()?.to_pyarray(py);
        array.readwrite().make_nonwriteable();
        Some(array)
    }

    /// The policy decisions the analysis applied, or `None` when it did not record them.
    #[getter]
    fn policy_reads(&self) -> Option<Vec<&'static str>> {
        self.reads
            .as_ref()
            .map(|fields| fields.iter().map(|field| field.name()).collect())
    }

    /// The quantity the analysis estimates, in words, or `None` when it did not say.
    #[getter]
    fn estimand(&self) -> Option<&str> {
        self.estimand.as_deref()
    }

    /// Why there is no answer, or `None` when there is one.
    #[getter]
    fn indeterminacy(&self) -> Option<&str> {
        self.reason.as_deref()
    }

    /// `complete`, `partial`, `ambiguous` or `indeterminate`.
    #[getter]
    const fn status(&self) -> &'static str {
        self.status
    }

    #[getter]
    const fn coverage(&self) -> PyCoverage {
        self.coverage
    }

    #[getter]
    fn warnings(&self) -> Vec<String> {
        self.warnings.clone()
    }

    /// Decisions made on the caller's behalf, one line each.
    #[getter]
    fn assumptions(&self) -> Vec<String> {
        self.assumptions.clone()
    }

    /// Deterministic JSON record of the input, the policy and the algorithm.
    #[getter]
    fn provenance(&self) -> &str {
        &self.provenance
    }

    /// The named policy profile when no decision was changed, else `None`.
    #[getter]
    const fn profile(&self) -> Option<&'static str> {
        self.profile
    }

    fn __repr__(&self) -> String {
        format!(
            "Analysis(status={}, coverage={:.3})",
            self.status,
            self.coverage.fraction()
        )
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAnalysis>()?;
    module.add_class::<PyCoverage>()
}

#[cfg(test)]
#[path = "analysis_result_tests.rs"]
mod tests;
