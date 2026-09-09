//! Governed coordinate, contact-area, CE and contact-map scores.

use super::common::{GovernedCompareError, complete, float};
use crate::governed_parameters::add_ce_parameters;
use crate::{
    CadScore, CeAlignment, CeOptions, CompareError, ContactArea, ContactSimilarity, LddtOptions,
    cad_contact_areas, cad_score, ce_align, ce_alignments, contact_map_similarity, gdt_ha, gdt_ts,
    gdt_with_cutoffs, lddt_with_options, tm_score, weighted_rmsd,
};
use pdbiox_core::ExecutionContext;
use pdbiox_core::contract::{Analysis, AnalysisPolicy, ParameterValue};

type CoordinateScore = fn(&[[f32; 3]], &[[f32; 3]]) -> Result<f64, CompareError>;

/// Constructs governed CAD residue contact areas from atom surfaces.
///
/// # Errors
///
/// Returns surface, length or coverage errors.
pub fn governed_cad_contact_areas(
    positions: &[[f32; 3]],
    radii: &[f32],
    residues: &[u32],
    probe: f32,
    density: f32,
    policy: &AnalysisPolicy,
    context: &ExecutionContext,
) -> Result<Analysis<Vec<ContactArea>>, GovernedCompareError> {
    let value = cad_contact_areas(positions, radii, residues, probe, density, context)?;
    let mut result = complete(value, positions.len(), policy, "cad-contact-areas")?;
    result.provenance = result
        .provenance
        .with_parameter("probe", float(probe))
        .with_parameter("density", float(density))
        .with_parameter(
            "radii",
            ParameterValue::Text("explicit-atom-aligned".into()),
        )
        .with_parameter(
            "area_convention",
            ParameterValue::Text("mean-directed-buried-patches".into()),
        );
    Ok(result)
}

/// Computes governed lDDT for an explicit correspondence.
///
/// # Errors
///
/// Returns a comparison error or a coverage overflow.
pub fn governed_lddt(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    options: &LddtOptions,
    policy: &AnalysisPolicy,
    context: &ExecutionContext,
) -> Result<Analysis<f64>, GovernedCompareError> {
    let value = lddt_with_options(model, reference, options, context)?;
    let mut result = complete(value, reference.len(), policy, "lddt")?;
    result.provenance = result
        .provenance
        .with_parameter("inclusion_radius", float(options.inclusion_radius))
        .with_parameter(
            "minimum_reference_distance",
            float(options.minimum_reference_distance),
        )
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
            "empty_policy",
            ParameterValue::Text(format!("{:?}", options.empty_policy).into()),
        );
    Ok(result)
}

fn governed_superposed(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    policy: &AnalysisPolicy,
    name: &'static str,
    measure: CoordinateScore,
) -> Result<Analysis<f64>, GovernedCompareError> {
    complete(measure(model, reference)?, reference.len(), policy, name)
}

/// Computes a governed TM-score.
///
/// # Errors
///
/// Returns a comparison error or a coverage overflow.
pub fn governed_tm_score(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    policy: &AnalysisPolicy,
) -> Result<Analysis<f64>, GovernedCompareError> {
    governed_superposed(model, reference, policy, "tm-score", tm_score)
}

/// Computes governed weighted RMSD in the caller-selected coordinate frame.
///
/// # Errors
///
/// Returns invalid weights, length mismatch or coverage overflow.
pub fn governed_weighted_rmsd(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    weights: &[f64],
    policy: &AnalysisPolicy,
) -> Result<Analysis<f64>, GovernedCompareError> {
    let value = weighted_rmsd(model, reference, weights)?;
    let mut result = complete(value, reference.len(), policy, "weighted-rmsd")?;
    result.provenance = result.provenance.with_parameter(
        "weights",
        ParameterValue::Text("explicit-atom-aligned".into()),
    );
    Ok(result)
}

/// Computes governed GDT-TS.
///
/// # Errors
///
/// Returns a comparison error or a coverage overflow.
pub fn governed_gdt_ts(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    policy: &AnalysisPolicy,
) -> Result<Analysis<f64>, GovernedCompareError> {
    governed_superposed(model, reference, policy, "gdt-ts", gdt_ts)
}

/// Computes governed GDT-HA.
///
/// # Errors
///
/// Returns a comparison error or a coverage overflow.
pub fn governed_gdt_ha(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    policy: &AnalysisPolicy,
) -> Result<Analysis<f64>, GovernedCompareError> {
    governed_superposed(model, reference, policy, "gdt-ha", gdt_ha)
}

/// Computes governed GDT over explicit, strictly increasing cutoffs.
///
/// # Errors
///
/// Returns a comparison error or coverage overflow.
pub fn governed_gdt_with_cutoffs(
    model: &[[f32; 3]],
    reference: &[[f32; 3]],
    cutoffs: &[f64],
    policy: &AnalysisPolicy,
) -> Result<Analysis<f64>, GovernedCompareError> {
    let value = gdt_with_cutoffs(model, reference, cutoffs)?;
    let mut result = complete(value, reference.len(), policy, "gdt")?;
    result.provenance = result.provenance.with_parameter(
        "cutoffs",
        ParameterValue::Text(
            cutoffs
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
                .into(),
        ),
    );
    Ok(result)
}

/// Computes a governed contact-area-difference score.
///
/// # Errors
///
/// Returns invalid-area or coverage errors.
pub fn governed_cad_score(
    reference: &[ContactArea],
    model: &[ContactArea],
    policy: &AnalysisPolicy,
) -> Result<Analysis<CadScore>, GovernedCompareError> {
    complete(cad_score(reference, model)?, reference.len(), policy, "cad")
}

/// Finds the best governed CE structural alignment.
///
/// # Errors
///
/// Returns CE input/search errors or coverage overflow.
pub fn governed_ce_align(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<CeAlignment>, GovernedCompareError> {
    let value = ce_align(reference, mobile, options)?;
    let mut result = complete(value, reference.len(), policy, "ce-align")?;
    add_ce_parameters(&mut result, options);
    Ok(result)
}

/// Returns every retained governed CE alignment candidate.
///
/// # Errors
///
/// Returns CE input/search errors or coverage overflow.
pub fn governed_ce_alignments(
    reference: &[[f32; 3]],
    mobile: &[[f32; 3]],
    options: CeOptions,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<CeAlignment>>, GovernedCompareError> {
    let value = ce_alignments(reference, mobile, options)?;
    let mut result = complete(value, reference.len(), policy, "ce-alignments")?;
    add_ce_parameters(&mut result, options);
    Ok(result)
}

/// Compares two explicit contact maps with governed provenance.
///
/// # Errors
///
/// Returns coverage overflow.
pub fn governed_contact_map_similarity(
    first: &[(u32, u32)],
    second: &[(u32, u32)],
    policy: &AnalysisPolicy,
) -> Result<Analysis<ContactSimilarity>, GovernedCompareError> {
    complete(
        contact_map_similarity(first, second),
        first.len(),
        policy,
        "contact-map-similarity",
    )
}
