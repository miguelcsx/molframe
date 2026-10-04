//! Policy audits: vary the decisions an analysis makes, run it under every
//! defensible combination, and measure how far the answer moves.
//!
//! The decisions are named in the vocabulary of a written policy. The audit
//! refuses a combination that contradicts itself and a decision the analysis
//! never applied, because a sweep over either would report an answer it has not
//! earned.

use crate::analysis_result::PyAnalysis;
use crate::audit_metric::{MetricKind, Projected};
use crate::policy::PyAnalysisPolicy;
use molframe::AnalysisPolicy;
use molframe::PolicyField;
use molframe::audit::{
    AuditPlan, AuditedRun, Decomposition, PlanError, PolicyDimension, PolicySpace, audit_analyses,
};
use pyo3::prelude::*;

/// One run as the audit sees it: the Python analysis, and what it reduced to.
struct Run {
    original: Py<PyAny>,
    projected: Option<Projected>,
    reads: Option<Vec<PolicyField>>,
}

impl AuditedRun for Run {
    type Value = Projected;

    fn answer(&self) -> Option<&Projected> {
        self.projected.as_ref()
    }

    fn policy_reads(&self) -> Option<Vec<PolicyField>> {
        self.reads.clone()
    }
}

fn field_of(name: &str) -> PyResult<PolicyField> {
    PolicyField::from_name(name).ok_or_else(|| {
        let known: Vec<&str> = PolicyField::ALL.iter().map(|field| field.name()).collect();
        crate::error::value(format!(
            "{name:?} is not a policy decision; the decisions are {}",
            known.join(", ")
        ))
    })
}

fn plan_error(error: &PlanError) -> PyErr {
    crate::error::kernel(error)
}

/// A baseline policy and the decisions to vary, before they are expanded into runs.
#[derive(Clone, Debug)]
#[pyclass(
    name = "PolicySpace",
    frozen,
    from_py_object,
    module = "molframe.audit"
)]
pub(crate) struct PyPolicySpace {
    baseline: AnalysisPolicy,
    dimensions: Vec<(PolicyField, Vec<String>)>,
    max_runs: usize,
}

impl PyPolicySpace {
    fn native(&self) -> PyResult<PolicySpace> {
        let mut space = PolicySpace::new(self.baseline.clone()).with_max_runs(self.max_runs);
        for (field, words) in &self.dimensions {
            let words: Vec<&str> = words.iter().map(String::as_str).collect();
            space =
                space.vary(PolicyDimension::named(*field, &words).map_err(crate::error::kernel)?);
        }
        Ok(space)
    }
}

#[pymethods]
impl PyPolicySpace {
    #[new]
    #[pyo3(signature = (baseline=None, *, max_runs=4096))]
    fn new(baseline: Option<PyRef<'_, PyAnalysisPolicy>>, max_runs: usize) -> Self {
        Self {
            baseline: crate::policy::policy_of(baseline),
            dimensions: Vec::new(),
            max_runs,
        }
    }

    /// A space that also varies `field` over `values`, written as in a policy.
    fn vary(&self, field: &str, values: Vec<String>) -> PyResult<Self> {
        let field = field_of(field)?;
        let words: Vec<&str> = values.iter().map(String::as_str).collect();
        PolicyDimension::named(field, &words).map_err(crate::error::kernel)?;
        let mut next = self.clone();
        next.dimensions.push((field, values));
        Ok(next)
    }

    /// The exact number of runs the space expands to.
    #[getter]
    fn cost(&self) -> PyResult<usize> {
        self.native()?.cost().map_err(|error| plan_error(&error))
    }

    /// Expands the space into runs, refusing contradictory combinations.
    fn plan(&self) -> PyResult<PyAuditPlan> {
        let plan = self.native()?.plan().map_err(|error| plan_error(&error))?;
        Ok(PyAuditPlan(plan))
    }

    fn __repr__(&self) -> String {
        let fields: Vec<&str> = self
            .dimensions
            .iter()
            .map(|(field, _)| field.name())
            .collect();
        format!("PolicySpace(varies=[{}])", fields.join(", "))
    }
}

/// The expanded runs of a policy space.
#[derive(Clone, Debug)]
#[pyclass(name = "AuditPlan", frozen, from_py_object, module = "molframe.audit")]
pub(crate) struct PyAuditPlan(AuditPlan);

#[pymethods]
impl PyAuditPlan {
    /// The number of runs.
    #[getter]
    fn cost(&self) -> usize {
        self.0.cost()
    }

    /// The varied decisions, in the order the plan lays them out.
    #[getter]
    fn fields(&self) -> Vec<&'static str> {
        self.0.fields().iter().map(|field| field.name()).collect()
    }

    /// The policy of every run.
    #[getter]
    fn policies(&self) -> Vec<PyAnalysisPolicy> {
        self.0
            .policies()
            .iter()
            .cloned()
            .map(PyAnalysisPolicy)
            .collect()
    }

    /// Each run's alternative in each varied decision, by position.
    #[getter]
    fn coordinates(&self) -> Vec<Vec<usize>> {
        self.0.coordinates().to_vec()
    }

    fn __len__(&self) -> usize {
        self.0.cost()
    }
}

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
    runs: Vec<Py<PyAny>>,
    policies: Vec<AnalysisPolicy>,
    metric: &'static str,
    read: Vec<&'static str>,
    indeterminate: Vec<usize>,
    indeterminate_fraction: f64,
    decomposition: Option<Decomposition>,
    agreement: Option<f64>,
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

/// Runs `analyse` under every policy of `plan` and measures how far the answers move.
///
/// `analyse` takes an `AnalysisPolicy` and returns a governed `Analysis`.
/// `project` reduces its value to what `metric` compares: items for `set` and
/// `ranking`, a number for `absolute` and `relative`, numbers for `rms` and
/// `correlation`, a category for `flip`, `(nodes, edges)` for the graph metrics;
/// the default is the value itself. A plan that varies a decision the analysis never
/// applied is refused after the first run.
#[pyfunction]
#[pyo3(signature = (plan, analyse, *, metric, project=None))]
fn run(
    py: Python<'_>,
    plan: &PyAuditPlan,
    analyse: &Bound<'_, PyAny>,
    metric: &str,
    project: Option<&Bound<'_, PyAny>>,
) -> PyResult<PyAuditResult> {
    let kind = MetricKind::parse(metric)?;
    let call = |policy: &AnalysisPolicy| -> PyResult<Run> {
        let original = analyse.call1((PyAnalysisPolicy(policy.clone()),))?;
        let analysis: PyRef<'_, PyAnalysis> = original
            .extract()
            .map_err(|_| crate::error::type_error("analyse must return an Analysis"))?;
        let reads = analysis.reads();
        let projected = match analysis.answer(py) {
            Some(value) => {
                let value = value.into_bound(py);
                let reduced = match project {
                    Some(project) => project.call1((value,))?,
                    None => value,
                };
                Some(kind.project(&reduced)?)
            }
            None => None,
        };
        drop(analysis);
        Ok(Run {
            original: original.unbind(),
            projected,
            reads,
        })
    };
    let audit = audit_analyses(&plan.0, call, &kind).map_err(|error| match error {
        molframe::audit::AuditError::Plan(error) => plan_error(&error),
        molframe::audit::AuditError::Analysis(error) => error,
        _ => crate::error::internal("the audit failed in a way this binding does not know"),
    })?;
    Ok(PyAuditResult {
        runs: audit
            .runs
            .iter()
            .map(|run| run.result.original.clone_ref(py))
            .collect(),
        policies: audit.runs.iter().map(|run| run.policy.clone()).collect(),
        metric: audit.metric,
        read: audit.read.iter().map(|field| field.name()).collect(),
        indeterminate: audit.indeterminate.clone(),
        indeterminate_fraction: audit.indeterminate_fraction(),
        agreement: audit.agreement_with_first(),
        decomposition: audit.decomposition,
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyPolicySpace>()?;
    module.add_class::<PyAuditPlan>()?;
    module.add_class::<PyEffect>()?;
    module.add_class::<PyInteraction>()?;
    module.add_class::<PyAuditResult>()?;
    module.add_function(wrap_pyfunction!(run, module)?)
}
