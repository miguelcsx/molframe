//! Governed chain, sequence and explicit point-mapping workflows.

use super::common::{GovernedCompareError, complete, covered, float};
use crate::{
    ChainAssignment, ChainMapping, ComparisonAlignment, ComparisonVerdict, DistanceMeasurement,
    InterfaceRmsd, PocketRmsd, PointMapping, ResidueMatch, align_mapping, assign_chains,
    decide_rmsd, interface_rmsd, map_chains, map_sequence_to_structure, measure_mapping,
    pocket_rmsd,
};
use molframe_chem::ComponentProvider;
use molframe_core::Structure;
use molframe_core::contract::{AlignmentPolicy, Analysis, AnalysisPolicy, ParameterValue};
use molframe_seq::Scoring;

/// Maps chains using the policy's identifier namespace.
///
/// # Errors
///
/// Returns CCD/identifier diagnostics or coverage overflow.
pub fn governed_map_chains(
    reference: &Structure,
    target: &Structure,
    provider: &dyn ComponentProvider,
    scoring: Scoring,
    minimum_identity: f64,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<ChainMapping>>, GovernedCompareError> {
    let value = map_chains(
        reference,
        target,
        provider,
        policy.identifiers,
        scoring,
        minimum_identity,
    )
    .map_err(GovernedCompareError::Diagnostic)?;
    let result = complete(
        value,
        reference.data().chains().count(),
        policy,
        "chain-mapping",
    )?;
    Ok(with_chain_mapping_parameters(
        result,
        scoring,
        minimum_identity,
    ))
}

/// Returns governed optimal chain assignment and viable alternatives.
///
/// # Errors
///
/// Returns CCD/identifier diagnostics or coverage overflow.
pub fn governed_assign_chains(
    reference: &Structure,
    target: &Structure,
    provider: &dyn ComponentProvider,
    scoring: Scoring,
    minimum_identity: f64,
    policy: &AnalysisPolicy,
) -> Result<Analysis<ChainAssignment>, GovernedCompareError> {
    let value = assign_chains(
        reference,
        target,
        provider,
        policy.identifiers,
        scoring,
        minimum_identity,
    )
    .map_err(GovernedCompareError::Diagnostic)?;
    let result = complete(
        value,
        reference.data().chains().count(),
        policy,
        "chain-assignment",
    )?;
    Ok(with_chain_mapping_parameters(
        result,
        scoring,
        minimum_identity,
    ))
}

fn with_chain_mapping_parameters<T>(
    mut result: Analysis<T>,
    scoring: Scoring,
    minimum_identity: f64,
) -> Analysis<T> {
    result.provenance = result
        .provenance
        .with_parameter("minimum_identity", float(minimum_identity))
        .with_parameter(
            "match_score",
            ParameterValue::Integer(i64::from(scoring.match_score)),
        )
        .with_parameter(
            "mismatch_score",
            ParameterValue::Integer(i64::from(scoring.mismatch_score)),
        )
        .with_parameter(
            "gap_open",
            ParameterValue::Integer(i64::from(scoring.gap_open)),
        )
        .with_parameter(
            "gap_extend",
            ParameterValue::Integer(i64::from(scoring.gap_extend)),
        );
    result
}

/// Maps a query sequence onto structure residues with explicit scoring.
///
/// # Errors
///
/// Returns CCD diagnostics or coverage overflow.
pub fn governed_map_sequence_to_structure(
    query: &[u8],
    structure: &Structure,
    provider: &dyn ComponentProvider,
    scoring: Scoring,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<ResidueMatch>>, GovernedCompareError> {
    let value = map_sequence_to_structure(query, structure, provider, policy.identifiers, scoring)
        .map_err(GovernedCompareError::Diagnostic)?;
    let used = value.len();
    covered(
        value,
        query.len(),
        used,
        policy,
        "sequence-structure-mapping",
    )
}

fn workflow_alignment(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    policy: &AnalysisPolicy,
) -> Result<ComparisonAlignment, GovernedCompareError> {
    match &policy.alignment {
        AlignmentPolicy::None => Ok(ComparisonAlignment::NotRequired),
        AlignmentPolicy::Explicit(_) | AlignmentPolicy::Global | AlignmentPolicy::Local => {
            Ok(align_mapping(reference, model, mapping)?)
        }
        _ => Err(GovernedCompareError::UnsupportedPolicy("alignment")),
    }
}

/// Executes governed mapping alignment, distance measurement and RMSD verdict.
///
/// `AlignmentPolicy::None` measures the existing coordinate frame;
/// `AlignmentPolicy::Explicit` fits the supplied mapping. Sequence-derived
/// global/local correspondence must first be constructed by the mapping APIs.
///
/// # Errors
///
/// Returns mapping-superposition, unsupported-policy or coverage errors.
pub fn governed_comparison_workflow(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    maximum_rmsd: f64,
    policy: &AnalysisPolicy,
) -> Result<Analysis<(DistanceMeasurement, ComparisonVerdict)>, GovernedCompareError> {
    let alignment = workflow_alignment(reference, model, mapping, policy)?;
    let measurement = measure_mapping(reference, model, mapping, alignment);
    let verdict = decide_rmsd(&measurement, maximum_rmsd);
    let mut result = complete(
        (measurement, verdict),
        mapping.matches().len(),
        policy,
        "comparison-workflow",
    )?;
    result.provenance = result
        .provenance
        .with_parameter("maximum_rmsd", float(maximum_rmsd));
    Ok(result)
}

/// Computes governed RMSD over an explicitly mapped interface.
///
/// # Errors
///
/// Returns alignment or coverage errors.
pub fn governed_interface_rmsd(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    policy: &AnalysisPolicy,
) -> Result<Analysis<InterfaceRmsd>, GovernedCompareError> {
    let alignment = workflow_alignment(reference, model, mapping, policy)?;
    complete(
        interface_rmsd(reference, model, mapping, alignment),
        mapping.matches().len(),
        policy,
        "interface-rmsd",
    )
}

/// Computes governed RMSD over an explicitly mapped binding pocket.
///
/// # Errors
///
/// Returns alignment or coverage errors.
pub fn governed_pocket_rmsd(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    mapping: &PointMapping,
    policy: &AnalysisPolicy,
) -> Result<Analysis<PocketRmsd>, GovernedCompareError> {
    let alignment = workflow_alignment(reference, model, mapping, policy)?;
    complete(
        pocket_rmsd(reference, model, mapping, alignment),
        mapping.matches().len(),
        policy,
        "pocket-rmsd",
    )
}
