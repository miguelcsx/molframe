//! Complete mapping, measurement and verdict execution.

use super::alignment::PyAlignmentKind;
use super::errors::evaluation_error;
use super::measurement::{PyMeasurementOptions, PyMeasurementSet};
use super::specification::PyMotif;
use super::verdict::{PyVerdict, PyVerdictProfile};
use crate::chemistry::PyComponentDictionary;
use crate::query::PyAnalysisPolicy;
use crate::structure::PyStructure;
use pyo3::prelude::*;

#[pyclass(name = "FxEvaluation", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEvaluation(pub(crate) molframe::fx::Evaluation);

#[pymethods]
impl PyEvaluation {
    #[getter]
    const fn mapping_index(&self) -> usize {
        self.0.mapping_index
    }

    #[getter]
    fn alignment(&self) -> PyAlignmentKind {
        self.0.alignment.into()
    }

    #[getter]
    fn measurements(&self) -> PyMeasurementSet {
        PyMeasurementSet(self.0.measurements.clone())
    }

    #[getter]
    fn verdict(&self) -> PyVerdict {
        PyVerdict(self.0.verdict.clone())
    }
}

#[pyclass(name = "EvaluationReport", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEvaluationReport(pub(crate) molframe::fx::EvaluationReport);

#[pymethods]
impl PyEvaluationReport {
    #[getter]
    const fn mapping_ambiguous(&self) -> bool {
        self.0.mapping_ambiguous
    }

    #[getter]
    fn evaluations(&self) -> Vec<PyEvaluation> {
        self.0
            .evaluations
            .iter()
            .cloned()
            .map(PyEvaluation)
            .collect()
    }
}

#[pyfunction]
#[pyo3(signature = (structure, motif, profile, mapping_limit, measurement, *, dictionary=None, policy=None))]
pub(crate) fn evaluate_motif(
    py: Python<'_>,
    structure: &PyStructure,
    motif: &PyMotif,
    profile: &PyVerdictProfile,
    mapping_limit: usize,
    measurement: &PyMeasurementOptions,
    dictionary: Option<&PyComponentDictionary>,
    policy: Option<&PyAnalysisPolicy>,
) -> PyResult<PyEvaluationReport> {
    let structure = structure.structure().clone();
    let motif = motif.0.clone();
    let profile = profile.0.clone();
    let measurement = measurement.0;
    let dictionary = dictionary.map(|value| value.0.clone());
    let policy = policy.map_or_else(molframe::AnalysisPolicy::default, |value| {
        value.inner.clone()
    });
    py.detach(move || {
        let provider = dictionary
            .as_ref()
            .map(|value| value.as_ref() as &dyn molframe::ComponentProvider);
        molframe::fx::evaluate_motif(
            &structure,
            &motif,
            provider,
            &policy,
            &profile,
            mapping_limit,
            measurement,
        )
        .map(PyEvaluationReport)
        .map_err(evaluation_error)
    })
}
