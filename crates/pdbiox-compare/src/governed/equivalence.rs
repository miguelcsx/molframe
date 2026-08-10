//! Governed chemistry-aware atom equivalence and ligand comparison.

use super::common::{GovernedCompareError, complete};
use crate::governed_parameters::integer;
use crate::{EquivalentAtomMapping, LigandRmsd, equivalent_atom_mappings, ligand_symmetry_rmsd};
use pdbiox_chem::Component;
use pdbiox_core::contract::{Analysis, AnalysisPolicy, EquivalencePolicy};

fn require_ccd(policy: &AnalysisPolicy) -> Result<(), GovernedCompareError> {
    match policy.atom_equivalence {
        EquivalencePolicy::Ccd => Ok(()),
        EquivalencePolicy::None | EquivalencePolicy::Explicit => {
            Err(GovernedCompareError::UnsupportedPolicy("atom_equivalence"))
        }
        _ => Err(GovernedCompareError::UnsupportedPolicy("atom_equivalence")),
    }
}

/// Enumerates governed exact component automorphisms.
///
/// # Errors
///
/// Returns an explicit mapping bound or coverage overflow.
pub fn governed_equivalent_atom_mappings(
    component: &Component,
    limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<Vec<EquivalentAtomMapping>>, GovernedCompareError> {
    require_ccd(policy)?;
    let value = equivalent_atom_mappings(component, limit)?;
    let mut result = complete(
        value,
        component.atoms.len(),
        policy,
        "equivalent-atom-mappings",
    )?;
    result.provenance = result
        .provenance
        .with_parameter("mapping_limit", integer(limit));
    Ok(result)
}

/// Computes governed exact symmetry-aware ligand RMSD.
///
/// # Errors
///
/// Returns correspondence, input or coverage errors.
pub fn governed_ligand_symmetry_rmsd(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    component: &Component,
    limit: usize,
    policy: &AnalysisPolicy,
) -> Result<Analysis<LigandRmsd>, GovernedCompareError> {
    require_ccd(policy)?;
    let value = ligand_symmetry_rmsd(reference, model, component, limit)?;
    let mut result = complete(value, reference.len(), policy, "ligand-symmetry-rmsd")?;
    result.provenance = result
        .provenance
        .with_parameter("mapping_limit", integer(limit));
    Ok(result)
}
