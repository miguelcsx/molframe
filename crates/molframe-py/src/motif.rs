//! Mechanical adapter for declarative functional-geometry evaluation.

use crate::bindings::PyStructure;
use crate::execution::{PyExecutionContext, run};
use crate::motif_stream::{mapping_error, prepare};
use crate::policy::{PyAnalysisPolicy, policy_of};
use molframe::motif::{
    Evaluation, EvaluationReport, MeasurementOptions, VerdictStatus, evaluate_plan,
};

/// Estimated bytes a materialised evaluation takes before its constraints.
const REPORT_BASE_BYTES: usize = 1024;
/// Estimated bytes each constraint adds to a materialised evaluation.
const REPORT_CONSTRAINT_BYTES: usize = 256;
use pyo3::prelude::*;
use std::path::PathBuf;

/// One constraint's measurement, with the intrinsic result kept apart from any verdict.
#[derive(Clone, Debug)]
#[pyclass(
    name = "ConstraintResult",
    frozen,
    skip_from_py_object,
    module = "molframe.motif"
)]
struct PyConstraintResult {
    name: String,
    value: Option<f64>,
    deviation: Option<f64>,
    satisfied: Option<bool>,
}

#[pymethods]
impl PyConstraintResult {
    #[getter]
    fn name(&self) -> &str {
        &self.name
    }

    /// The observed value; `None` when no alternative was measurable.
    #[getter]
    const fn value(&self) -> Option<f64> {
        self.value
    }

    /// Distance from the acceptable region; zero means satisfied.
    #[getter]
    const fn deviation(&self) -> Option<f64> {
        self.deviation
    }

    /// Whether the constraint holds; `None` when it is indeterminate.
    #[getter]
    const fn satisfied(&self) -> Option<bool> {
        self.satisfied
    }

    fn __repr__(&self) -> String {
        format!(
            "ConstraintResult(name={:?}, value={:?}, satisfied={:?})",
            self.name, self.value, self.satisfied
        )
    }
}

/// One mapping of the motif onto the structure and the verdict it earns.
#[derive(Clone, Debug)]
#[pyclass(
    name = "MotifEvaluation",
    frozen,
    skip_from_py_object,
    module = "molframe.motif"
)]
struct PyMotifEvaluation {
    mapping_index: usize,
    profile: String,
    verdict: &'static str,
    constraints: Vec<PyConstraintResult>,
}

#[pymethods]
impl PyMotifEvaluation {
    /// Position of this mapping in the deterministic enumeration.
    #[getter]
    const fn mapping_index(&self) -> usize {
        self.mapping_index
    }

    /// The verdict profile that decided: its exact identifier.
    #[getter]
    fn profile(&self) -> &str {
        &self.profile
    }

    /// `"pass"`, `"fail"` or `"indeterminate"`.
    #[getter]
    const fn verdict(&self) -> &'static str {
        self.verdict
    }

    /// Per-constraint measurements in the motif's order; there is no collapsed score.
    #[getter]
    fn constraints(&self) -> Vec<PyConstraintResult> {
        self.constraints.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "MotifEvaluation(mapping={}, verdict={})",
            self.mapping_index, self.verdict
        )
    }
}

/// Every defensible mapping of a motif, each evaluated.
#[derive(Clone, Debug)]
#[pyclass(
    name = "MotifReport",
    frozen,
    skip_from_py_object,
    module = "molframe.motif"
)]
struct PyMotifReport {
    mapping_ambiguous: bool,
    evaluations: Vec<PyMotifEvaluation>,
}

#[pymethods]
impl PyMotifReport {
    /// Whether chemistry and structure admitted more than one mapping.
    #[getter]
    const fn mapping_ambiguous(&self) -> bool {
        self.mapping_ambiguous
    }

    #[getter]
    fn evaluations(&self) -> Vec<PyMotifEvaluation> {
        self.evaluations.clone()
    }

    fn __len__(&self) -> usize {
        self.evaluations.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "MotifReport(evaluations={}, mapping_ambiguous={})",
            self.evaluations.len(),
            self.mapping_ambiguous
        )
    }
}

const fn verdict_name(status: VerdictStatus) -> &'static str {
    match status {
        VerdictStatus::Pass => "pass",
        VerdictStatus::Fail => "fail",
        VerdictStatus::Indeterminate => "indeterminate",
    }
}

fn evaluation(evaluation: &Evaluation) -> PyMotifEvaluation {
    PyMotifEvaluation {
        mapping_index: evaluation.mapping_index,
        profile: evaluation.verdict.profile.to_string(),
        verdict: verdict_name(evaluation.verdict.status),
        constraints: evaluation
            .measurements
            .constraints
            .iter()
            .map(|measurement| PyConstraintResult {
                name: measurement.name.to_string(),
                value: measurement
                    .value
                    .map(molframe::motif::MeasurementValue::numeric),
                deviation: measurement.deviation,
                satisfied: measurement.satisfied,
            })
            .collect(),
    }
}

fn report(report: &EvaluationReport) -> PyMotifReport {
    PyMotifReport {
        mapping_ambiguous: report.mapping_ambiguous,
        evaluations: report.evaluations.iter().map(evaluation).collect(),
    }
}

/// Evaluates the motif and verdict profile a specification file declares.
///
/// The specification is TOML or JSON. When a Chemical Component Dictionary is
/// given, atom mapping is chemistry-aware, so interchangeable atoms are not
/// penalised for the way a depositor happened to label them.
///
/// `limits` is `(mappings, alternatives)`. `mappings=None` (the default) sets
/// no mapping cap; the call is instead bounded by the context's memory budget,
/// refused before anything is allocated. `alternatives` bounds the
/// equivalent-atom alternatives measured. `constrained=True` keeps only
/// mappings whose distance constraints are satisfied, found with a spatial
/// index, so cost follows the output. To avoid materialising reports at all,
/// use `count`, `matches` or `iter_matches`.
#[pyfunction]
#[pyo3(signature = (
    structure,
    specification,
    *,
    components=None,
    version="unversioned",
    limits=(None, 256),
    constrained=false,
    policy=None,
    context=None
))]
#[allow(clippy::too_many_arguments)]
fn evaluate(
    py: Python<'_>,
    structure: &PyStructure,
    specification: PathBuf,
    components: Option<PathBuf>,
    version: &str,
    limits: (Option<usize>, usize),
    constrained: bool,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
    context: Option<&PyExecutionContext>,
) -> PyResult<PyMotifReport> {
    let (mapping_limit, maximum_alternatives) = limits;
    if mapping_limit == Some(0) || maximum_alternatives == 0 {
        return Err(crate::error::value(
            "limits must be positive: (mappings, alternatives); use None for no mapping cap",
        ));
    }
    let source = structure.inner.clone();
    let policy = policy_of(policy);
    let measurement = MeasurementOptions {
        maximum_alternatives,
        plane_fit: molframe::geometry::EigenOptions::standard(),
    };
    let version = version.to_owned();
    run(py, context, move |ctx| {
        let (specification, plan) = prepare(
            &source,
            &specification,
            components,
            &version,
            &policy,
            constrained,
        )?;
        let total = plan
            .count(ctx.cancellation())
            .map_err(|error| mapping_error(&error))?;
        let total = match usize::try_from(total) {
            Ok(total) => total,
            Err(_) => usize::MAX,
        };
        if let Some(cap) = mapping_limit.filter(|&cap| total > cap) {
            return Err(crate::error::value(format!(
                "motif mapping exceeded the limit of {cap} mappings ({total} found); raise \
                 limits=(mappings, alternatives), pass limits=(None, ...) for no cap, \
                 use constrained=True to keep only distance-satisfying mappings, \
                 or stream with motif.iter_matches / motif.count"
            )));
        }
        let bytes = total.saturating_mul(
            REPORT_BASE_BYTES + REPORT_CONSTRAINT_BYTES * specification.motif.constraints().len(),
        );
        let _reservation = ctx.try_reserve(bytes).map_err(|_| {
            crate::error::named(
                "MemoryBudgetError",
                format!(
                    "evaluating {total} mappings needs about {bytes} bytes, over the execution \
                     memory budget of {} bytes; raise it with ExecutionContext(memory_budget=...) \
                     passed as context=, use constrained=True to keep only distance-satisfying \
                     mappings, or stream with motif.iter_matches / motif.count",
                    ctx.memory_budget().bytes()
                ),
            )
        })?;
        evaluate_plan(
            source.engine(),
            &specification.motif,
            &plan,
            &specification.profile,
            usize::MAX,
            measurement,
            ctx.cancellation(),
        )
        .map(|result| report(&result))
        .map_err(|error| match error {
            molframe::motif::EvaluationError::Mapping(inner) => mapping_error(&inner),
            other => crate::error::failure(other),
        })
    })?
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyConstraintResult>()?;
    module.add_class::<PyMotifEvaluation>()?;
    module.add_class::<PyMotifReport>()?;
    module.add_function(wrap_pyfunction!(evaluate, module)?)?;
    crate::motif_stream::register(module)
}
