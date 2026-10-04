//! Auditing results of any kind, not only sets of items.
//!
//! The question an audit asks is how far the answer moved when a decision
//! changed. What "far" means depends on the answer, so the caller supplies an
//! [`OutcomeMetric`] and the audit supplies the design: which runs, which pairs
//! differ in one decision only, and how the variation splits among decisions and
//! their interactions.

use crate::design::Decomposition;
use crate::metric::OutcomeMetric;
use crate::{AuditPlan, AuditRun, PlanError};
use molframe_core::contract::{
    Analysis, AnalysisPolicy, Assumption, Impact, MeasuredImpact, PolicyField, Provenance,
};

/// A completed audit of results compared by one metric.
#[derive(Clone, Debug)]
pub struct OutcomeAudit<R> {
    /// Completed runs in the plan's deterministic order.
    pub runs: Vec<AuditRun<R>>,
    /// The metric the distances are in.
    pub metric: &'static str,
    /// How the variation among runs splits among decisions and their pairs.
    pub decomposition: Decomposition,
    from_first: Vec<f64>,
}

impl<R> OutcomeAudit<R> {
    /// The fraction of runs whose result equals the first run's.
    ///
    /// The plan's first run takes every dimension's first alternative, so a
    /// caller who lists the decision they would have made first reads this as the
    /// fraction of defensible universes that agree with their own choice.
    #[must_use]
    pub fn agreement_with_first(&self) -> f64 {
        agreement(&self.from_first)
    }

    /// What varying `field` alone did to the result, as a measured impact.
    #[must_use]
    pub fn impact_of(&self, field: PolicyField) -> Option<MeasuredImpact> {
        impact(&self.decomposition, field, self.metric)
    }

    /// Fills in the impact of every assumption whose field was varied.
    pub fn measure(&self, assumptions: &mut [Assumption]) {
        measure(&self.decomposition, self.metric, assumptions);
    }
}

fn agreement(from_first: &[f64]) -> f64 {
    if from_first.is_empty() {
        return 1.0;
    }
    let same = from_first
        .iter()
        .filter(|&&distance| distance == 0.0)
        .count();
    crate::numeric::usize_to_f64(same) / crate::numeric::usize_to_f64(from_first.len())
}

fn impact(
    decomposition: &Decomposition,
    field: PolicyField,
    metric: &'static str,
) -> Option<MeasuredImpact> {
    let effect = decomposition
        .main_effects
        .iter()
        .find(|effect| effect.field == field)?;
    Some(MeasuredImpact {
        metric,
        change: effect.mean_change,
        comparisons: u32::try_from(effect.comparisons).ok()?,
    })
}

fn measure(decomposition: &Decomposition, metric: &'static str, assumptions: &mut [Assumption]) {
    for assumption in assumptions {
        if let Some(measured) = impact(decomposition, assumption.field, metric) {
            assumption.impact = Impact::Measured(measured);
        }
    }
}

/// Runs a plan and compares the results with `metric`.
///
/// # Errors
///
/// Returns the first error `analyse` produces; later policy points are not run.
pub fn audit_outcomes<R, E, A, M>(
    plan: &AuditPlan,
    mut analyse: A,
    metric: &M,
) -> Result<OutcomeAudit<R>, E>
where
    A: FnMut(&AnalysisPolicy) -> Result<R, E>,
    M: OutcomeMetric<R>,
{
    let mut runs = Vec::with_capacity(plan.cost());
    for policy in plan.policies() {
        runs.push(AuditRun {
            policy: policy.clone(),
            result: analyse(policy)?,
        });
    }
    let distance =
        |first: usize, second: usize| metric.distance(&runs[first].result, &runs[second].result);
    let decomposition = plan.decompose(distance);
    let from_first = (0..runs.len())
        .map(|run| metric.distance(&runs[0].result, &runs[run].result))
        .collect();
    Ok(OutcomeAudit {
        runs,
        metric: metric.name(),
        decomposition,
        from_first,
    })
}

/// A run of an analysis the audit can interrogate: whether it has an answer, what the
/// answer is, and which policy fields the analysis applied.
///
/// Implemented for [`Analysis`]; a binding that holds its results in another shape
/// implements it to get the same audit.
pub trait AuditedRun {
    /// The answer, when there is one.
    type Value;

    /// The answer, or `None` when the analysis declined.
    fn answer(&self) -> Option<&Self::Value>;

    /// The policy fields the analysis recorded as applied, when it recorded them.
    fn policy_reads(&self) -> Option<Vec<PolicyField>>;

    /// The record of what produced this run, when it kept one.
    ///
    /// A certificate lists the inputs, algorithm and exact policy of every universe
    /// from these records; a run without one is listed without them.
    fn provenance(&self) -> Option<&Provenance> {
        None
    }
}

impl<T> AuditedRun for Analysis<T> {
    type Value = T;

    fn answer(&self) -> Option<&T> {
        self.value()
    }

    fn policy_reads(&self) -> Option<Vec<PolicyField>> {
        self.provenance.policy_reads()
    }

    fn provenance(&self) -> Option<&Provenance> {
        Some(&self.provenance)
    }
}

/// Why an audit of governed analyses could not be completed.
#[derive(Debug)]
#[non_exhaustive]
pub enum AuditError<E> {
    /// The plan varies a field the analysis does not apply, or is otherwise invalid.
    Plan(PlanError),
    /// The analysis itself failed.
    Analysis(E),
}

impl<E: std::fmt::Display> std::fmt::Display for AuditError<E> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Plan(error) => error.fmt(formatter),
            Self::Analysis(error) => write!(formatter, "analysis failed: {error}"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for AuditError<E> {}

/// A completed audit of governed analyses.
///
/// Some universes may have no defensible answer: an analysis that declines under
/// one assembly choice and answers under another has changed its conclusion in
/// the starkest way there is. That is reported as a fraction, not folded into a
/// distance, and the decomposition is withheld while any run is indeterminate,
/// because a variance split needs an answer in every cell.
#[derive(Debug)]
pub struct AnalysisAudit<R> {
    /// Completed runs in the plan's deterministic order.
    pub runs: Vec<AuditRun<R>>,
    /// The metric the distances are in.
    pub metric: &'static str,
    /// The policy fields the analysis applied across the runs.
    pub read: Vec<PolicyField>,
    /// Zero-based indices of runs with no defensible answer.
    pub indeterminate: Vec<usize>,
    /// The split of variation, when every run has an answer.
    pub decomposition: Option<Decomposition>,
    from_first: Option<Vec<f64>>,
}

impl<R> AnalysisAudit<R> {
    /// The fraction of defensible universes in which the analysis has no answer.
    #[must_use]
    pub fn indeterminate_fraction(&self) -> f64 {
        if self.runs.is_empty() {
            0.0
        } else {
            crate::numeric::usize_to_f64(self.indeterminate.len())
                / crate::numeric::usize_to_f64(self.runs.len())
        }
    }

    /// The fraction of runs whose answer equals the first run's, when all have one.
    #[must_use]
    pub fn agreement_with_first(&self) -> Option<f64> {
        self.from_first.as_deref().map(agreement)
    }

    /// What varying `field` alone did to the answer, when all runs have one.
    #[must_use]
    pub fn impact_of(&self, field: PolicyField) -> Option<MeasuredImpact> {
        impact(self.decomposition.as_ref()?, field, self.metric)
    }

    /// Fills in the measured impact of every assumption whose field was varied.
    pub fn measure(&self, assumptions: &mut [Assumption]) {
        if let Some(decomposition) = &self.decomposition {
            measure(decomposition, self.metric, assumptions);
        }
    }
}

/// Runs a plan over a governed analysis and compares the answers with `metric`.
///
/// The first run decides whether the audit means anything: the analysis records
/// which policy fields it applied, and a plan that varies another is refused
/// before the rest are paid for. An analysis that records nothing applied nothing.
///
/// # Errors
///
/// Returns [`AuditError::Plan`] for a varied field the analysis does not apply,
/// and [`AuditError::Analysis`] for the first failure of `analyse`.
pub fn audit_analyses<R, E, A, M>(
    plan: &AuditPlan,
    mut analyse: A,
    metric: &M,
) -> Result<AnalysisAudit<R>, AuditError<E>>
where
    R: AuditedRun,
    A: FnMut(&AnalysisPolicy) -> Result<R, E>,
    M: OutcomeMetric<R::Value>,
{
    let mut runs: Vec<AuditRun<R>> = Vec::with_capacity(plan.cost());
    let mut read: Vec<PolicyField> = Vec::new();
    for (index, policy) in plan.policies().iter().enumerate() {
        let result = analyse(policy).map_err(AuditError::Analysis)?;
        if let Some(applied) = result.policy_reads() {
            for field in applied {
                if !read.contains(&field) {
                    read.push(field);
                }
            }
        }
        if index == 0 {
            plan.require_read(&read).map_err(AuditError::Plan)?;
        }
        runs.push(AuditRun {
            policy: policy.clone(),
            result,
        });
    }
    let indeterminate: Vec<usize> = runs
        .iter()
        .enumerate()
        .filter_map(|(index, run)| run.result.answer().is_none().then_some(index))
        .collect();
    let (decomposition, from_first) = if indeterminate.is_empty() {
        let values: Vec<&R::Value> = runs.iter().filter_map(|run| run.result.answer()).collect();
        let distance = |first: usize, second: usize| metric.distance(values[first], values[second]);
        let split = plan.decompose(distance);
        let from_first: Vec<f64> = (0..values.len())
            .map(|run| metric.distance(values[0], values[run]))
            .collect();
        (Some(split), Some(from_first))
    } else {
        (None, None)
    };
    Ok(AnalysisAudit {
        runs,
        metric: metric.name(),
        read,
        indeterminate,
        decomposition,
        from_first,
    })
}

#[cfg(test)]
#[path = "outcomes_tests.rs"]
mod tests;
