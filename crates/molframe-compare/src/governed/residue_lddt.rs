//! Governed residue-aware lDDT.

use super::common::{GovernedCompareError, covered, float};
use crate::mapped::MappedCompareError;
use crate::{ResidueCorrespondence, ResidueLddt, ResidueLddtOptions, SymmetryPolicy, residue_lddt};
use molframe_chem::ComponentProvider;
use molframe_core::ExecutionContext;
use molframe_core::contract::{Analysis, AnalysisPolicy, ParameterValue};
use molframe_core::structure::Structure;

fn text(value: &'static str) -> ParameterValue {
    ParameterValue::Text(value.into())
}

/// Computes a governed residue-aware lDDT over an explicit residue pairing.
///
/// Coverage counts reference residues: a residue the model lacks is intended
/// but not used, and its reference distances still count as not preserved.
/// The provenance records every choice that shapes the number, including the
/// conventions the definition fixes rather than exposes.
///
/// # Errors
///
/// Returns a mapping or scoring failure, or a coverage overflow.
pub fn governed_residue_lddt(
    reference: &Structure,
    model: &Structure,
    residues: &[ResidueCorrespondence],
    provider: &dyn ComponentProvider,
    options: &ResidueLddtOptions,
    policy: &AnalysisPolicy,
    context: &ExecutionContext,
) -> Result<Analysis<ResidueLddt>, GovernedCompareError> {
    let value =
        residue_lddt(reference, model, residues, provider, options, context).map_err(|error| {
            match error {
                MappedCompareError::Mapping(diagnostic) => {
                    GovernedCompareError::Diagnostic(diagnostic)
                }
                MappedCompareError::Compare(error) => GovernedCompareError::Compare(error),
            }
        })?;
    let used = value
        .residues
        .iter()
        .filter(|residue| residue.model.is_some())
        .count();
    let intended = value.residues.len();
    let mut result = covered(value, intended, used, policy, "residue-lddt")?;
    let symmetry = match options.symmetry {
        SymmetryPolicy::Disabled => "disabled".to_owned(),
        SymmetryPolicy::Resolve { automorphism_limit } => {
            format!("topology-only;limit={automorphism_limit}")
        }
    };
    result.provenance = result
        .provenance
        .with_parameter("inclusion_radius", float(options.inclusion_radius))
        .with_parameter(
            "tolerances",
            ParameterValue::Text(
                options
                    .tolerances
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
                    .into(),
            ),
        )
        .with_parameter(
            "sequence_separation",
            ParameterValue::Text(options.sequence_separation.to_string().into()),
        )
        .with_parameter("symmetry", ParameterValue::Text(symmetry.into()))
        .with_parameter(
            "empty_policy",
            ParameterValue::Text(format!("{:?}", options.empty_policy).into()),
        )
        .with_parameter("reference_pairs", text("heavy-atom-inter-residue"))
        .with_parameter("missing_model_atoms", text("not-preserved"))
        .with_parameter("pair_counting", text("both-sides"));
    Ok(result)
}
