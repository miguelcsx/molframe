//! What an audit reports: effects, interactions, shares and the universes with no answer.

use crate::policy::PyAnalysisPolicy;
use molframe::AnalysisPolicy;
use molframe::audit::Decomposition;
use pyo3::prelude::*;

/// What one decision explains of the variation among runs.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Effect",
    frozen,
    skip_from_py_object,
    module = "molframe.audit"
)]
pub(crate) struct PyEffect {
    #[pyo3(get)]
    field: &'static str,
    #[pyo3(get)]
    mean_change: f64,
    #[pyo3(get)]
    comparisons: usize,
    #[pyo3(get)]
    share: f64,
}

/// What two decisions explain together beyond their separate effects.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Interaction",
    frozen,
    skip_from_py_object,
    module = "molframe.audit"
)]
pub(crate) struct PyInteraction {
    #[pyo3(get)]
    first: &'static str,
    #[pyo3(get)]
    second: &'static str,
    #[pyo3(get)]
    share: f64,
}

/// A completed audit.
#[pyclass(
    name = "AuditResult",
    frozen,
    skip_from_py_object,
    module = "molframe.audit"
)]
pub(crate) struct PyAuditResult {
    pub(crate) certificate: String,
    pub(crate) runs: Vec<Py<PyAny>>,
    pub(crate) policies: Vec<AnalysisPolicy>,
    pub(crate) metric: &'static str,
    pub(crate) read: Vec<&'static str>,
    pub(crate) indeterminate: Vec<usize>,
    pub(crate) indeterminate_fraction: f64,
    pub(crate) decomposition: Option<Decomposition>,
    pub(crate) agreement: Option<f64>,
}

/// A decision's or a class's Shapley share of the variation among runs.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Attribution",
    frozen,
    skip_from_py_object,
    module = "molframe.audit"
)]
pub(crate) struct PyAttribution {
    #[pyo3(get)]
    name: &'static str,
    #[pyo3(get)]
    share: f64,
}

/// A varied decision: its class of uncertainty and the case for varying it.
#[derive(Clone, Debug)]
#[pyclass(
    name = "Decision",
    frozen,
    skip_from_py_object,
    module = "molframe.audit"
)]
pub(crate) struct PyDecision {
    #[pyo3(get)]
    pub(crate) field: &'static str,
    #[pyo3(get)]
    pub(crate) uncertainty: &'static str,
    #[pyo3(get)]
    pub(crate) rationale: String,
    #[pyo3(get)]
    pub(crate) evidence: String,
}

#[pymethods]
impl PyAuditResult {
    /// The analysis of every run, in the plan's order.
    #[getter]
    fn runs(&self, py: Python<'_>) -> Vec<Py<PyAny>> {
        self.runs.iter().map(|run| run.clone_ref(py)).collect()
    }

    /// The policy of every run.
    #[getter]
    fn policies(&self) -> Vec<PyAnalysisPolicy> {
        self.policies
            .iter()
            .cloned()
            .map(PyAnalysisPolicy)
            .collect()
    }

    /// The distance the numbers are in.
    /// The audit as an RO-Crate 1.1 metadata document (JSON-LD): the inputs and their
    /// SHA-256, the algorithm and what it estimates, every decision with its class,
    /// rationale and evidence, the exact policy and provenance of every universe, and what
    /// was measured. Write it as ``ro-crate-metadata.json``.
    ///
    /// It says how far the answer moved when the listed decisions changed; it does not say
    /// which choice is right.
    #[getter]
    fn certificate(&self) -> &str {
        &self.certificate
    }

    #[getter]
    const fn metric(&self) -> &'static str {
        self.metric
    }

    /// The decisions the analysis applied, as it recorded them.
    #[getter]
    fn read(&self) -> Vec<&'static str> {
        self.read.clone()
    }

    /// The runs that had no defensible answer.
    #[getter]
    fn indeterminate(&self) -> Vec<usize> {
        self.indeterminate.clone()
    }

    /// The fraction of defensible universes in which there is no answer.
    #[getter]
    const fn indeterminate_fraction(&self) -> f64 {
        self.indeterminate_fraction
    }

    /// The fraction of runs whose answer equals the first run's, or `None` when
    /// some run has no answer.
    #[getter]
    const fn agreement_with_first(&self) -> Option<f64> {
        self.agreement
    }

    /// What each decision explains alone; `None` while any run has no answer.
    #[getter]
    fn effects(&self) -> Option<Vec<PyEffect>> {
        let decomposition = self.decomposition.as_ref()?;
        Some(
            decomposition
                .main_effects
                .iter()
                .map(|effect| PyEffect {
                    field: effect.field.name(),
                    mean_change: effect.mean_change,
                    comparisons: effect.comparisons,
                    share: effect.share,
                })
                .collect(),
        )
    }

    /// Whether every combination of the varied decisions was run. Only then do the
    /// effect and interaction shares add up; the Shapley shares do not need it.
    #[getter]
    fn balanced(&self) -> Option<bool> {
        self.decomposition
            .as_ref()
            .map(|decomposition| decomposition.balanced)
    }

    /// Each decision's Shapley share of the variation; they sum to one.
    #[getter]
    fn shapley(&self) -> Option<Vec<PyAttribution>> {
        let decomposition = self.decomposition.as_ref()?;
        Some(
            decomposition
                .shapley
                .iter()
                .map(|attribution| PyAttribution {
                    name: attribution.key.name(),
                    share: attribution.share,
                })
                .collect(),
        )
    }

    /// The same attribution over the classes of uncertainty (structural,
    /// interpretive, algorithmic, numerical) among the varied decisions.
    #[getter]
    fn by_class(&self) -> Option<Vec<PyAttribution>> {
        let decomposition = self.decomposition.as_ref()?;
        Some(
            decomposition
                .by_class
                .iter()
                .map(|attribution| PyAttribution {
                    name: attribution.key.name(),
                    share: attribution.share,
                })
                .collect(),
        )
    }

    /// What each pair of decisions explains beyond their separate effects.
    #[getter]
    fn interactions(&self) -> Option<Vec<PyInteraction>> {
        let decomposition = self.decomposition.as_ref()?;
        Some(
            decomposition
                .interactions
                .iter()
                .map(|interaction| PyInteraction {
                    first: interaction.first.name(),
                    second: interaction.second.name(),
                    share: interaction.share,
                })
                .collect(),
        )
    }

    /// The share left to combinations of three or more decisions.
    #[getter]
    fn higher_order(&self) -> Option<f64> {
        self.decomposition
            .as_ref()
            .map(|decomposition| decomposition.higher_order)
    }

    /// Total variation among runs, in squared metric units.
    #[getter]
    fn total_variation(&self) -> Option<f64> {
        self.decomposition
            .as_ref()
            .map(|decomposition| decomposition.total_variation)
    }

    /// Mean distance over every pair of runs.
    #[getter]
    fn mean_distance(&self) -> Option<f64> {
        self.decomposition
            .as_ref()
            .map(|decomposition| decomposition.mean_distance)
    }

    /// Largest distance between any two runs.
    #[getter]
    fn max_distance(&self) -> Option<f64> {
        self.decomposition
            .as_ref()
            .map(|decomposition| decomposition.max_distance)
    }

    fn __repr__(&self) -> String {
        format!(
            "AuditResult(runs={}, metric={}, indeterminate={})",
            self.runs.len(),
            self.metric,
            self.indeterminate.len()
        )
    }
}
