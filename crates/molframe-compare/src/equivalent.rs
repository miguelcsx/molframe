//! Exact chemistry-aware atom correspondence and ligand comparison.

use crate::CompareError;
use crate::numeric::usize_to_f64;
use molframe_chem::{Component, automorphisms};

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
    if atoms == 0 {
        return Err(CompareError::EmptyComponent);
    }
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

/// The part of a component whose atoms satisfy `present`.
///
/// Equivalences are those of the fragment that is observed: a structure
/// without hydrogens has no use for the permutations of a methyl group's
/// hydrogens, and counting them would exhaust the automorphism bound on
/// ordinary side chains. Bonds survive only between retained atoms, so an
/// atom whose distinguishing neighbour is missing is not wrongly merged with
/// another branch. Ideal coordinates are dropped because they index the full
/// atom list.
#[must_use]
pub fn component_fragment(component: &Component, present: impl Fn(&str) -> bool) -> Component {
    let atoms: Vec<_> = component
        .atoms
        .iter()
        .filter(|atom| present(&atom.name))
        .cloned()
        .collect();
    let bonds: Vec<_> = component
        .bonds
        .iter()
        .filter(|bond| present(&bond.atom_a) && present(&bond.atom_b))
        .cloned()
        .collect();
    Component {
        atoms: atoms.into(),
        bonds: bonds.into(),
        ideal_coordinates: None,
        model_coordinates: None,
        ..component.clone()
    }
}

/// Symmetry-aware ligand RMSD between atoms identified by name.
///
/// Only the atoms named in both coordinate sets and defined by `component`
/// take part, so a model that lacks hydrogens, or a ligand with missing atoms,
/// is compared over what both have and automorphisms are those of that
/// fragment. Coordinates must already share a frame. The returned mapping
/// indexes the fragment's atoms in component order.
///
/// # Errors
///
/// Returns [`CompareError::EmptyComponent`] when no atom is shared and an
/// explicit bound error when the fragment has more automorphisms than `limit`.
pub fn named_ligand_rmsd(
    reference: &[(&str, [f32; 3])],
    model: &[(&str, [f32; 3])],
    component: &Component,
    limit: usize,
) -> Result<LigandRmsd, CompareError> {
    let first = |atoms: &[(&str, [f32; 3])], name: &str| {
        atoms
            .iter()
            .find(|(candidate, _)| *candidate == name)
            .map(|(_, position)| *position)
    };
    let shared = |name: &str| first(reference, name).is_some() && first(model, name).is_some();
    let fragment = component_fragment(component, shared);
    let gather = |atoms: &[(&str, [f32; 3])]| -> Vec<[f32; 3]> {
        fragment
            .atoms
            .iter()
            .filter_map(|atom| first(atoms, &atom.name))
            .collect()
    };
    ligand_symmetry_rmsd(&gather(reference), &gather(model), &fragment, limit)
}

#[cfg(test)]
#[path = "equivalent_tests.rs"]
mod tests;
