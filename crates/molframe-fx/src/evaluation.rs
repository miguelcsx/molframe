use crate::{
    AlignmentKind, MappingError, MeasurementError, MeasurementOptions, MeasurementSet, Motif,
    Verdict, VerdictProfile, align_intrinsic, map_motif, measure_constraints,
};
use molframe_chem::ComponentProvider;
use molframe_core::contract::AnalysisPolicy;
use molframe_core::structure::Structure;

/// One complete mapping → alignment → measurement → decision path.
#[derive(Clone, Debug, PartialEq)]
pub struct Evaluation {
    /// Mapping ordinal in deterministic enumeration order.
    pub mapping_index: usize,
    /// Alignment decision used by the measurement stage.
    pub alignment: AlignmentKind,
    /// Per-constraint measurements without a collapsed score.
    pub measurements: MeasurementSet,
    /// Profile decision and its individual rule outcomes.
    pub verdict: Verdict,
}

/// Evaluation of every defensible motif mapping.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct EvaluationReport {
    /// Whether chemistry and structure admitted multiple mappings.
    pub mapping_ambiguous: bool,
    /// One evaluation per complete mapping, in stable order.
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
    let mappings = map_motif(structure, motif, provider, policy, mapping_limit)?;
    let mut evaluations = Vec::with_capacity(mappings.mappings.len());
    for (mapping_index, mapping) in mappings.mappings.into_iter().enumerate() {
        let aligned = align_intrinsic(mapping);
        let measurements = measure_constraints(structure, motif, &aligned, measurement)?;
        let verdict = profile.decide(&measurements.metrics());
        evaluations.push(Evaluation {
            mapping_index,
            alignment: aligned.kind,
            measurements,
            verdict,
        });
    }
    Ok(EvaluationReport {
        mapping_ambiguous: mappings.ambiguous,
        evaluations,
    })
}

#[cfg(test)]
#[path = "evaluation_tests.rs"]
mod tests;
