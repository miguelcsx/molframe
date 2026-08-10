//! Governed interface-comparison scores.

use super::common::{GovernedCompareError, complete, float};
use crate::{DockQ, DockQOptions, QsOptions, dockq_in_namespace, qs_score_in_namespace};
use pdbiox_core::Structure;
use pdbiox_core::contract::{Analysis, AnalysisPolicy, ParameterValue};

/// Computes governed `DockQ` for explicitly named partners.
///
/// # Errors
///
/// Returns a comparison error or a coverage overflow.
pub fn governed_dockq(
    model: &Structure,
    reference: &Structure,
    receptor: &str,
    ligand: &str,
    options: DockQOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<DockQ>, GovernedCompareError> {
    let value = dockq_in_namespace(
        model,
        reference,
        receptor,
        ligand,
        policy.identifiers,
        options,
    )?;
    let mut result = complete(value, reference.atom_count() as usize, policy, "dockq")?;
    result.provenance = result
        .provenance
        .with_parameter("receptor", ParameterValue::Text(receptor.into()))
        .with_parameter("ligand", ParameterValue::Text(ligand.into()))
        .with_parameter("contact_distance", float(options.contact_distance))
        .with_parameter("ligand_scale", float(options.ligand_scale))
        .with_parameter("interface_scale", float(options.interface_scale));
    Ok(result)
}

/// Computes a governed quaternary-structure interface score.
///
/// # Errors
///
/// Returns a comparison error or a coverage overflow.
pub fn governed_qs_score(
    model: &Structure,
    reference: &Structure,
    first_chain: &str,
    second_chain: &str,
    options: QsOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<f64>, GovernedCompareError> {
    let value = qs_score_in_namespace(
        model,
        reference,
        first_chain,
        second_chain,
        policy.identifiers,
        options,
    )?;
    let mut result = complete(value, reference.atom_count() as usize, policy, "qs-score")?;
    result.provenance = result
        .provenance
        .with_parameter("first_chain", ParameterValue::Text(first_chain.into()))
        .with_parameter("second_chain", ParameterValue::Text(second_chain.into()))
        .with_parameter("contact_distance", float(options.contact_distance))
        .with_parameter(
            "empty_policy",
            ParameterValue::Text(format!("{:?}", options.empty_policy).into()),
        );
    Ok(result)
}
