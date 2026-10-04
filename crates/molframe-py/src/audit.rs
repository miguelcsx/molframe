//! Policy audits: vary the decisions an analysis makes, run it under every
//! defensible combination, and measure how far the answer moves.
//!
//! The decisions are named in the vocabulary of a written policy. The audit
//! refuses a combination that contradicts itself and a decision the analysis
//! never applied, because a sweep over either would report an answer it has not
//! earned.

use crate::analysis_result::PyAnalysis;
use crate::audit_metric::{MetricKind, Projected};
use crate::audit_result::{PyAuditResult, PyDecision};
use crate::policy::PyAnalysisPolicy;
use molframe::AnalysisPolicy;
use molframe::PolicyField;
use molframe::audit::{
    AuditPlan, AuditedRun, PlanError, PolicyDimension, PolicySpace, PolicyValue, audit_analyses,
};
use pyo3::prelude::*;

/// One run as the audit sees it: the Python analysis, and what it reduced to.
struct Run {
    original: Py<PyAny>,
    projected: Option<Projected>,
    reads: Option<Vec<PolicyField>>,
    record: molframe::Provenance,
}

impl AuditedRun for Run {
    type Value = Projected;

    fn answer(&self) -> Option<&Projected> {
        self.projected.as_ref()
    }

    fn policy_reads(&self) -> Option<Vec<PolicyField>> {
        self.reads.clone()
    }

    fn provenance(&self) -> Option<&molframe::Provenance> {
        Some(&self.record)
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
    dimensions: Vec<Dimension>,
    exclusions: Vec<((PolicyField, String), (PolicyField, String))>,
    max_runs: usize,
}

/// One varied decision as written: its alternatives, and why they are defensible.
#[derive(Clone, Debug)]
struct Dimension {
    field: PolicyField,
    words: Vec<String>,
    rationale: String,
    evidence: String,
}

impl PyPolicySpace {
    fn native(&self) -> PyResult<PolicySpace> {
        let mut space = PolicySpace::new(self.baseline.clone()).with_max_runs(self.max_runs);
        for dimension in &self.dimensions {
            let words: Vec<&str> = dimension.words.iter().map(String::as_str).collect();
            space = space.vary(
                PolicyDimension::named(dimension.field, &words)
                    .map_err(crate::error::kernel)?
                    .justified(dimension.rationale.as_str(), dimension.evidence.as_str()),
            );
        }
        for ((when_field, when), (then_field, then)) in &self.exclusions {
            space = space.forbid(
                PolicyValue::named(*when_field, when).map_err(crate::error::kernel)?,
                PolicyValue::named(*then_field, then).map_err(crate::error::kernel)?,
            );
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
            exclusions: Vec::new(),
            max_runs,
        }
    }

    /// A space that also varies `field` over `values`, written as in a policy.
    ///
    /// `rationale` says why those alternatives are the defensible ones and `evidence`
    /// what supports that; both travel with the plan and the result.
    #[pyo3(signature = (field, values, *, rationale="", evidence=""))]
    fn vary(
        &self,
        field: &str,
        values: Vec<String>,
        rationale: &str,
        evidence: &str,
    ) -> PyResult<Self> {
        let field = field_of(field)?;
        let words: Vec<&str> = values.iter().map(String::as_str).collect();
        PolicyDimension::named(field, &words).map_err(crate::error::kernel)?;
        let mut next = self.clone();
        next.dimensions.push(Dimension {
            field,
            words: values,
            rationale: rationale.to_owned(),
            evidence: evidence.to_owned(),
        });
        Ok(next)
    }

    /// A space in which a universe holding `when` cannot also hold `then_not`.
    ///
    /// Each is a `(decision, value)` pair written as in a policy. A strict plan refuses a
    /// space with such a pair; a constrained plan drops the universes that hold it.
    fn forbid(&self, when: (String, String), then_not: (String, String)) -> PyResult<Self> {
        let when_field = field_of(&when.0)?;
        let then_field = field_of(&then_not.0)?;
        PolicyValue::named(when_field, &when.1).map_err(crate::error::kernel)?;
        PolicyValue::named(then_field, &then_not.1).map_err(crate::error::kernel)?;
        let mut next = self.clone();
        next.exclusions
            .push(((when_field, when.1), (then_field, then_not.1)));
        Ok(next)
    }

    /// The exact number of runs the space expands to.
    #[getter]
    fn cost(&self) -> PyResult<usize> {
        self.native()?.cost().map_err(|error| plan_error(&error))
    }

    /// Expands the space into runs.
    ///
    /// The plan is the whole product of the alternatives, so a combination that cannot run
    /// or that the space forbids is refused. With `constrained=True` such combinations are
    /// dropped and counted instead, and the plan is no longer the whole product: the effect
    /// and interaction shares do not apply, and the Shapley shares carry the attribution.
    #[pyo3(signature = (*, constrained=false))]
    fn plan(&self, constrained: bool) -> PyResult<PyAuditPlan> {
        let space = self.native()?;
        let plan = if constrained {
            space.plan_constrained()
        } else {
            space.plan()
        }
        .map_err(|error| plan_error(&error))?;
        Ok(PyAuditPlan(plan))
    }

    fn __repr__(&self) -> String {
        let fields: Vec<&str> = self
            .dimensions
            .iter()
            .map(|dimension| dimension.field.name())
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

    /// Whether the plan is the whole product of the alternatives.
    #[getter]
    fn balanced(&self) -> bool {
        self.0.balanced()
    }

    /// How many combinations a constrained plan dropped.
    #[getter]
    fn skipped(&self) -> usize {
        self.0.skipped()
    }

    /// The varied decisions with their class of uncertainty and justification.
    #[getter]
    fn decisions(&self) -> Vec<PyDecision> {
        self.0
            .decisions()
            .iter()
            .map(|decision| PyDecision {
                field: decision.field.name(),
                uncertainty: decision.class.name(),
                rationale: decision.rationale.to_string(),
                evidence: decision.evidence.to_string(),
            })
            .collect()
    }

    fn __len__(&self) -> usize {
        self.0.cost()
    }
}

/// Runs `analyse` under every policy of `plan` and measures how far the answers move.
///
/// `analyse` takes an `AnalysisPolicy` and returns a governed `Analysis`.
/// `project` reduces the whole `Analysis` (its value, and for instance its `atom_origin`) to
/// what `metric` compares: items for `set` and
/// `ranking`, a number for `absolute` and `relative`, numbers for `rms` and
/// `correlation`, a category for `flip`, `(nodes, edges)` for the graph metrics;
/// the default is its value. A plan that varies a decision the analysis never
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
        let record = analysis.record();
        let projected = match analysis.answer(py) {
            Some(value) => {
                let reduced = match project {
                    Some(project) => project.call1((original.clone(),))?,
                    None => value.into_bound(py),
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
            record,
        })
    };
    let audit = audit_analyses(&plan.0, call, &kind).map_err(|error| match error {
        molframe::audit::AuditError::Plan(error) => plan_error(&error),
        molframe::audit::AuditError::Analysis(error) => error,
        _ => crate::error::internal("the audit failed in a way this binding does not know"),
    })?;
    let certificate = molframe::audit::certificate(&plan.0, &audit);
    Ok(PyAuditResult {
        certificate,
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
    module.add_class::<crate::audit_result::PyEffect>()?;
    module.add_class::<crate::audit_result::PyInteraction>()?;
    module.add_class::<crate::audit_result::PyAttribution>()?;
    module.add_class::<PyDecision>()?;
    module.add_class::<PyAuditResult>()?;
    module.add_function(wrap_pyfunction!(run, module)?)
}
