//! Exact chemistry-aware atom correspondence and ligand comparison.

use crate::CompareError;
use crate::numeric::usize_to_f64;
use pdbiox_chem::{Component, automorphisms};

/// One exact component-graph atom correspondence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EquivalentAtomMapping {
    /// Target atom position for every reference component atom position.
    pub reference_to_model: Box<[u32]>,
}

/// Enumerates exact chemically valid atom mappings for a component.
///
/// # Errors
///
/// Returns a bound error rather than a truncated mapping set.
pub fn equivalent_atom_mappings(
    component: &Component,
    limit: usize,
) -> Result<Vec<EquivalentAtomMapping>, CompareError> {
    automorphisms(component, limit)
        .map(|mappings| {
            mappings
                .into_iter()
                .map(|reference_to_model| EquivalentAtomMapping { reference_to_model })
                .collect()
        })
        .map_err(|error| CompareError::MappingLimit { limit: error.limit })
}

/// Best symmetry-aware coordinate RMSD and the exact mapping that produced it.
#[derive(Clone, Debug, PartialEq)]
pub struct LigandRmsd {
    /// Root-mean-square distance in ångström.
    pub rmsd: f64,
    /// Chemically valid component automorphism applied to the model.
    pub mapping: EquivalentAtomMapping,
}

/// Compares already aligned ligand coordinates over exact graph automorphisms.
///
/// Coordinates must follow component atom order. Ties retain the
/// lexicographically first mapping, making the result reproducible.
///
/// # Errors
///
/// Returns a length mismatch or an explicit automorphism-bound error.
pub fn ligand_symmetry_rmsd(
    reference: &[[f32; 3]],
    model: &[[f32; 3]],
    component: &Component,
    limit: usize,
) -> Result<LigandRmsd, CompareError> {
    let atoms = component.atoms.len();
    if reference.len() != atoms || model.len() != atoms {
        return Err(CompareError::LengthMismatch {
            model: model.len(),
            reference: reference.len(),
        });
    }
    let mappings = equivalent_atom_mappings(component, limit)?;
    let mut best: Option<LigandRmsd> = None;
    for mapping in mappings {
        let squared: f64 = reference
            .iter()
            .enumerate()
            .map(|(index, position)| {
                let target = mapping.reference_to_model[index] as usize;
                position
                    .iter()
                    .zip(model[target])
                    .map(|(left, right)| {
                        let delta = f64::from(*left - right);
                        delta * delta
                    })
                    .sum::<f64>()
            })
            .sum();
        let candidate = LigandRmsd {
            rmsd: (squared / usize_to_f64(atoms)).sqrt(),
            mapping,
        };
        if best
            .as_ref()
            .is_none_or(|current| candidate.rmsd < current.rmsd)
        {
            best = Some(candidate);
        }
    }
    best.ok_or(CompareError::EmptyComponent)
}

#[cfg(test)]
#[path = "equivalent_tests.rs"]
mod tests;
