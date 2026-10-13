use crate::{
    AlignmentKind, MappingError, MappingPlan, MeasurementError, MeasurementOptions, MeasurementSet,
    Motif, Verdict, VerdictProfile, align_intrinsic, measure_constraints,
};
use molframe_chem::ComponentProvider;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::execution::CancellationToken;
use molframe_core::structure::Structure;
use std::sync::Arc;

/// One complete mapping → alignment → measurement → decision path.
#[derive(Clone, Debug, PartialEq)]
pub struct Evaluation {
    /// Mapping ordinal in deterministic enumeration order.
    pub mapping_index: usize,
    /// Residue indices of the mapping, one per motif component in sorted-name
    /// order. Unlike `mapping_index`, this identifies the mapping whatever the
    /// plan's [`crate::Admission`] pruned away.
    pub residues: Box<[u32]>,
    /// Alignment decision used by the measurement stage.
    pub alignment: AlignmentKind,
    /// Per-constraint measurements without a collapsed score.
    pub measurements: MeasurementSet,
    /// Profile decision and its individual rule outcomes.
    pub verdict: Verdict,
}

/// Evaluation of every mapping a plan admits.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct EvaluationReport {
    /// Whether the plan admitted multiple mappings.
    pub mapping_ambiguous: bool,
    /// One evaluation per admitted mapping, in stable order.
    pub evaluations: Vec<Evaluation>,
}

/// Failure in a named evaluation stage.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum EvaluationError {
    /// Chemistry-aware mapping failed.
    #[error("mapping stage failed: {0}")]
    Mapping(#[from] MappingError),
    /// Bounded measurement failed.
    #[error("measurement stage failed: {0}")]
    Measurement(#[from] MeasurementError),
}

/// Runs all four functional-evaluation stages without selecting one ambiguous mapping.
///
/// Internal motif constraints are invariant under rigid transforms, so the
/// alignment stage records [`AlignmentKind::NotRequired`]. Callers evaluating
/// frame-dependent extensions can instead use `align_with_transform` followed
/// by `measure_constraints` and `VerdictProfile::decide` explicitly.
///
/// # Errors
///
/// Returns the stage-specific mapping or measurement error. No partial report
/// is returned after an explicit enumeration bound is exceeded.
pub fn evaluate_motif(
    structure: &Structure,
    motif: &Motif,
    provider: Option<&dyn ComponentProvider>,
    policy: &AnalysisPolicy,
    profile: &VerdictProfile,
    mapping_limit: usize,
    measurement: MeasurementOptions,
) -> Result<EvaluationReport, EvaluationError> {
    let plan = Arc::new(MappingPlan::build(
        structure, motif, provider, policy, false,
    )?);
    evaluate_plan(
        structure,
        motif,
        &plan,
        profile,
        mapping_limit,
        measurement,
        &CancellationToken::new(),
    )
}

/// Evaluates every mapping a prepared [`MappingPlan`] admits.
///
/// The mappings are counted first, so a plan over the limit is refused before
/// any evaluation is allocated.
///
/// # Errors
///
/// Returns the stage-specific error, [`MappingError::LimitExceeded`] when the
/// plan admits more than `mapping_limit` mappings, or
/// [`MappingError::Cancelled`] when `cancel` fires.
pub fn evaluate_plan(
    structure: &Structure,
    motif: &Motif,
    plan: &Arc<MappingPlan>,
    profile: &VerdictProfile,
    mapping_limit: usize,
    measurement: MeasurementOptions,
    cancel: &CancellationToken,
) -> Result<EvaluationReport, EvaluationError> {
    let total = plan.count(cancel)?;
    if total > mapping_limit as u64 {
        return Err(MappingError::LimitExceeded {
            limit: mapping_limit,
        }
        .into());
    }
    let mut stream = plan.stream();
    let mut tuple = vec![0u32; plan.width()];
    let mut evaluations = Vec::with_capacity(match usize::try_from(total) {
        Ok(total) => total,
        Err(_) => 0,
    });
    while stream.next_into(&mut tuple, cancel)? {
        let aligned = align_intrinsic(plan.mapping(&tuple));
        let measurements = measure_constraints(structure, motif, &aligned, measurement)?;
        let verdict = profile.decide(&measurements.metrics());
        evaluations.push(Evaluation {
            mapping_index: evaluations.len(),
            residues: tuple.as_slice().into(),
            alignment: aligned.kind,
            measurements,
            verdict,
        });
    }
    Ok(EvaluationReport {
        mapping_ambiguous: total > 1,
        evaluations,
    })
}

#[cfg(test)]
#[path = "evaluation_tests.rs"]
mod tests;
