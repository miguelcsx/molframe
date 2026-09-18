//! Frozen compatibility profiles used by functional benchmarks.

use super::verdict::{PyVerdict, PyVerdictStatus};
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "ProfileAtomSet", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyProfileAtomSet {
    FullScaffoldCAlpha,
    MotifBackboneWithOxygen,
    CatalyticBackbone,
    CatalyticHeavyAtoms,
    LigandAndBackbone,
    Unknown,
}

impl From<molframe::fx::ProfileAtomSet> for PyProfileAtomSet {
    fn from(value: molframe::fx::ProfileAtomSet) -> Self {
        match value {
            molframe::fx::ProfileAtomSet::FullScaffoldCAlpha => Self::FullScaffoldCAlpha,
            molframe::fx::ProfileAtomSet::MotifBackboneWithOxygen => Self::MotifBackboneWithOxygen,
            molframe::fx::ProfileAtomSet::CatalyticBackbone => Self::CatalyticBackbone,
            molframe::fx::ProfileAtomSet::CatalyticHeavyAtoms => Self::CatalyticHeavyAtoms,
            molframe::fx::ProfileAtomSet::LigandAndBackbone => Self::LigandAndBackbone,
            _ => Self::Unknown,
        }
    }
}

#[pyclass(name = "ProfileAlignment", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyProfileAlignment {
    MeasuredAtoms,
    CatalyticBackbone,
    PredictionFrame,
    Unknown,
}

impl From<molframe::fx::ProfileAlignment> for PyProfileAlignment {
    fn from(value: molframe::fx::ProfileAlignment) -> Self {
        match value {
            molframe::fx::ProfileAlignment::MeasuredAtoms => Self::MeasuredAtoms,
            molframe::fx::ProfileAlignment::CatalyticBackbone => Self::CatalyticBackbone,
            molframe::fx::ProfileAlignment::PredictionFrame => Self::PredictionFrame,
            _ => Self::Unknown,
        }
    }
}

#[pyclass(name = "ProfileMetric", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyProfileMetric {
    #[pyo3(get)]
    name: String,
    #[pyo3(get)]
    atoms: PyProfileAtomSet,
    #[pyo3(get)]
    alignment: PyProfileAlignment,
}

impl From<molframe::fx::ProfileMetric> for PyProfileMetric {
    fn from(value: molframe::fx::ProfileMetric) -> Self {
        Self {
            name: value.name.to_owned(),
            atoms: value.atoms.into(),
            alignment: value.alignment.into(),
        }
    }
}

#[pyclass(name = "CandidateAggregation", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyCandidateAggregation {
    Any,
}

impl From<molframe::fx::CandidateAggregation> for PyCandidateAggregation {
    fn from(value: molframe::fx::CandidateAggregation) -> Self {
        match value {
            molframe::fx::CandidateAggregation::Any => Self::Any,
        }
    }
}

#[pyclass(name = "CompatibilityProfile", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCompatibilityProfile(pub(crate) molframe::fx::CompatibilityProfile);

#[pymethods]
impl PyCompatibilityProfile {
    #[getter]
    fn id(&self) -> &str {
        self.0.id()
    }

    #[getter]
    fn metrics(&self) -> Vec<PyProfileMetric> {
        self.0.metrics().iter().copied().map(Into::into).collect()
    }

    #[getter]
    fn aggregation(&self) -> PyCandidateAggregation {
        self.0.aggregation().into()
    }

    fn decide_candidate(&self, metrics: BTreeMap<String, f64>) -> PyVerdict {
        let metrics = native_metrics(metrics);
        PyVerdict(self.0.decide_candidate(&metrics))
    }

    fn decide_candidates(&self, candidates: Vec<BTreeMap<String, f64>>) -> PyCompatibilityVerdict {
        let candidates = candidates
            .into_iter()
            .map(native_metrics)
            .collect::<Vec<_>>();
        PyCompatibilityVerdict(self.0.decide_candidates(candidates.iter()))
    }
}

#[pyclass(name = "CompatibilityVerdict", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyCompatibilityVerdict(pub(crate) molframe::fx::CompatibilityVerdict);

#[pymethods]
impl PyCompatibilityVerdict {
    #[getter]
    fn profile(&self) -> &str {
        &self.0.profile
    }

    #[getter]
    fn status(&self) -> PyVerdictStatus {
        self.0.status.into()
    }

    #[getter]
    fn candidates(&self) -> Vec<PyVerdict> {
        self.0.candidates.iter().cloned().map(PyVerdict).collect()
    }
}

#[pyfunction]
pub(crate) fn motifbench_1_0(py: Python<'_>) -> PyCompatibilityProfile {
    py.detach(move || -> PyCompatibilityProfile {
        PyCompatibilityProfile(molframe::fx::motifbench_1_0())
    })
}

#[pyfunction]
pub(crate) fn ame_heavy_atom_1_0(py: Python<'_>) -> PyCompatibilityProfile {
    py.detach(move || -> PyCompatibilityProfile {
        PyCompatibilityProfile(molframe::fx::ame_heavy_atom_1_0())
    })
}

fn native_metrics(metrics: BTreeMap<String, f64>) -> BTreeMap<Box<str>, f64> {
    metrics
        .into_iter()
        .map(|(name, value)| (name.into_boxed_str(), value))
        .collect()
}
