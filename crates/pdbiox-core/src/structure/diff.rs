//! Semantic comparison of immutable structure snapshots.

use super::{Structure, UnitCell};
use crate::index::{AtomIndex, ChainIndex, ModelIndex, ResidueIndex};

/// A value before and after a semantic change.
#[derive(Clone, Debug, PartialEq)]
pub struct ValueDifference<T> {
    /// Value in the left-hand structure.
    pub left: T,
    /// Value in the right-hand structure.
    pub right: T,
}

/// A row count before and after a structural change.
pub type CountDifference = ValueDifference<usize>;

/// Entry-level metadata changes.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MetadataDifference {
    /// Deposited entry identifier.
    pub id: Option<ValueDifference<Option<Box<str>>>>,
    /// Entry title.
    pub title: Option<ValueDifference<Option<Box<str>>>>,
    /// Experimental method.
    pub method: Option<ValueDifference<Option<Box<str>>>>,
    /// Reported resolution in ångström.
    pub resolution: Option<ValueDifference<Option<f32>>>,
    /// Crystallographic unit cell.
    pub cell: Option<ValueDifference<Option<UnitCell>>>,
}

/// Explicit controls for semantic structure comparison.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StructureDifferenceOptions {
    /// Maximum Euclidean coordinate displacement considered unchanged, in ångström.
    pub coordinate_tolerance: f32,
}

/// Why a structure comparison could not be performed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DifferenceError {
    /// Coordinate tolerance must be finite and non-negative.
    InvalidCoordinateTolerance,
}

impl std::fmt::Display for DifferenceError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidCoordinateTolerance => {
                formatter.write_str("coordinate tolerance must be finite and non-negative")
            }
        }
    }
}

impl std::error::Error for DifferenceError {}

/// Complete semantic difference summary for two normalised structures.
#[derive(Clone, Debug, PartialEq)]
pub struct StructureDifference {
    /// Entry-level metadata that changed.
    pub metadata: MetadataDifference,
    /// Model counts when different.
    pub models: Option<CountDifference>,
    /// Entity counts when different.
    pub entities: Option<CountDifference>,
    /// Chain counts when different.
    pub chains: Option<CountDifference>,
    /// Residue counts when different.
    pub residues: Option<CountDifference>,
    /// Atom counts when different.
    pub atoms: Option<CountDifference>,
    /// Bond counts when different.
    pub bonds: Option<CountDifference>,
    /// Corresponding entity rows with different semantic values.
    pub changed_entities: usize,
    /// Corresponding chain rows with different semantic values.
    pub changed_chains: usize,
    /// Corresponding residue rows with different semantic values.
    pub changed_residues: usize,
    /// Corresponding atom rows with different non-coordinate values.
    pub changed_atoms: usize,
    /// Corresponding bond rows that differ.
    pub changed_bonds: usize,
    /// Corresponding model/atom positions beyond the requested tolerance.
    pub changed_positions: usize,
    /// Largest finite corresponding coordinate displacement, in ångström.
    pub maximum_displacement: Option<f32>,
}

impl StructureDifference {
    /// Whether the two normalised structures are semantically equivalent.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.metadata == MetadataDifference::default()
            && self.models.is_none()
            && self.entities.is_none()
            && self.chains.is_none()
            && self.residues.is_none()
            && self.atoms.is_none()
            && self.bonds.is_none()
            && self.changed_entities == 0
            && self.changed_chains == 0
            && self.changed_residues == 0
            && self.changed_atoms == 0
            && self.changed_bonds == 0
            && self.changed_positions == 0
    }
}

/// Compares all normalised structure data with explicit coordinate tolerance.
///
/// # Errors
///
/// Returns [`DifferenceError::InvalidCoordinateTolerance`] when the tolerance
/// is negative or not finite.
pub fn structure_difference(
    left: &Structure,
    right: &Structure,
    options: StructureDifferenceOptions,
) -> Result<StructureDifference, DifferenceError> {
    if !options.coordinate_tolerance.is_finite() || options.coordinate_tolerance < 0.0 {
        return Err(DifferenceError::InvalidCoordinateTolerance);
    }

    let (changed_positions, maximum_displacement) =
        coordinate_changes(left, right, options.coordinate_tolerance);
    Ok(StructureDifference {
        metadata: metadata_difference(left, right),
        models: different_count(left.model_count(), right.model_count()),
        entities: different_count(left.entity_count(), right.entity_count()),
        chains: different_count(left.chain_count(), right.chain_count()),
        residues: different_count(left.residue_count(), right.residue_count()),
        atoms: different_count(left.atom_count() as usize, right.atom_count() as usize),
        bonds: different_count(left.data().bonds.len(), right.data().bonds.len()),
        changed_entities: changed_entity_count(left, right),
        changed_chains: changed_chain_count(left, right),
        changed_residues: changed_residue_count(left, right),
        changed_atoms: changed_atom_count(left, right),
        changed_bonds: left
            .data()
            .bonds
            .iter()
            .zip(right.data().bonds.iter())
            .filter(|(left, right)| left != right)
            .count(),
        changed_positions,
        maximum_displacement,
    })
}

fn changed_entity_count(left: &Structure, right: &Structure) -> usize {
    let left_table = &left.data().topology.entities;
    let right_table = &right.data().topology.entities;
    left_table
        .iter()
        .zip(right_table.iter())
        .filter(|(left_index, right_index)| {
            left_table.kind(*left_index) != right_table.kind(*right_index)
                || resolved(left, left_table.id(*left_index))
                    != resolved(right, right_table.id(*right_index))
                || resolved(left, left_table.description(*left_index))
                    != resolved(right, right_table.description(*right_index))
                || resolved_sequence(left, left_table.canonical_sequence(*left_index))
                    != resolved_sequence(right, right_table.canonical_sequence(*right_index))
        })
        .count()
}

fn changed_chain_count(left: &Structure, right: &Structure) -> usize {
    (0..left.chain_count().min(right.chain_count()))
        .filter(|position| {
            let Ok(position) = u32::try_from(*position) else {
                return true;
            };
            let index = ChainIndex::new(position);
            let left_chain = left.chain(index);
            let right_chain = right.chain(index);
            match (left_chain, right_chain) {
                (Some(left_chain), Some(right_chain)) => {
                    left_chain.label() != right_chain.label()
                        || left_chain.auth_label() != right_chain.auth_label()
                        || left_chain.entity() != right_chain.entity()
                        || left_chain.polymer_kind() != right_chain.polymer_kind()
                }
                _ => true,
            }
        })
        .count()
}

fn changed_residue_count(left: &Structure, right: &Structure) -> usize {
    (0..left.residue_count().min(right.residue_count()))
        .filter(|position| {
            let Ok(position) = u32::try_from(*position) else {
                return true;
            };
            let index = ResidueIndex::new(position);
            let left_residue = left.residue(index);
            let right_residue = right.residue(index);
            match (left_residue, right_residue) {
                (Some(left_residue), Some(right_residue)) => {
                    left_residue.name() != right_residue.name()
                        || left_residue.auth_name() != right_residue.auth_name()
                        || left_residue.label_seq_id() != right_residue.label_seq_id()
                        || left_residue.auth_seq_id() != right_residue.auth_seq_id()
                        || left_residue.ins_code() != right_residue.ins_code()
                        || left_residue.is_het() != right_residue.is_het()
                }
                _ => true,
            }
        })
        .count()
}

fn changed_atom_count(left: &Structure, right: &Structure) -> usize {
    (0..left.atom_count().min(right.atom_count()))
        .filter(|position| {
            let left_atom = left.atom(AtomIndex::new(*position));
            let right_atom = right.atom(AtomIndex::new(*position));
            match (left_atom, right_atom) {
                (Some(left_atom), Some(right_atom)) => {
                    left_atom.name() != right_atom.name()
                        || left_atom.auth_name() != right_atom.auth_name()
                        || left_atom.component_name() != right_atom.component_name()
                        || left_atom.element() != right_atom.element()
                        || left_atom.alt_label() != right_atom.alt_label()
                        || left_atom.b_factor() != right_atom.b_factor()
                        || left_atom.occupancy() != right_atom.occupancy()
                        || left_atom.formal_charge() != right_atom.formal_charge()
                        || left_atom.atom_site_id() != right_atom.atom_site_id()
                }
                _ => true,
            }
        })
        .count()
}

fn coordinate_changes(left: &Structure, right: &Structure, tolerance: f32) -> (usize, Option<f32>) {
    let mut changed = 0usize;
    let mut maximum: Option<f32> = None;
    for model in 0..left.model_count().min(right.model_count()) {
        let Ok(model) = u32::try_from(model) else {
            break;
        };
        let index = ModelIndex::new(model);
        let left_snapshot = left.model_snapshot(index);
        let right_snapshot = right.model_snapshot(index);
        let left_positions: &[[f32; 3]] = match &left_snapshot {
            Some((structure, local)) => structure
                .model_positions(*local)
                .map_or(&[], |positions| positions),
            None => &[],
        };
        let right_positions: &[[f32; 3]] = match &right_snapshot {
            Some((structure, local)) => structure
                .model_positions(*local)
                .map_or(&[], |positions| positions),
            None => &[],
        };
        for (left_position, right_position) in left_positions.iter().zip(right_positions) {
            let displacement = squared_distance(*left_position, *right_position).sqrt();
            if !displacement.is_finite() || displacement > tolerance {
                changed += 1;
            }
            if displacement.is_finite() {
                maximum = Some(maximum.map_or(displacement, |value| value.max(displacement)));
            }
        }
        changed += left_positions.len().abs_diff(right_positions.len());
    }
    (changed, maximum)
}

fn metadata_difference(left: &Structure, right: &Structure) -> MetadataDifference {
    let left_entry = &left.data().entry;
    let right_entry = &right.data().entry;
    MetadataDifference {
        id: different_value(left_entry.id.clone(), right_entry.id.clone()),
        title: different_value(left_entry.title.clone(), right_entry.title.clone()),
        method: different_value(left_entry.method.clone(), right_entry.method.clone()),
        resolution: different_value(left_entry.resolution, right_entry.resolution),
        cell: different_value(left.data().cell, right.data().cell),
    }
}

fn different_count(left: usize, right: usize) -> Option<CountDifference> {
    different_value(left, right)
}

fn different_value<T: PartialEq>(left: T, right: T) -> Option<ValueDifference<T>> {
    (left != right).then_some(ValueDifference { left, right })
}

fn resolved(structure: &Structure, symbol: Option<crate::SymbolId>) -> Option<&str> {
    symbol.and_then(|value| structure.resolve(value))
}

fn resolved_sequence<'a>(
    structure: &'a Structure,
    sequence: &'a [crate::SymbolId],
) -> Vec<Option<&'a str>> {
    sequence
        .iter()
        .map(|symbol| structure.resolve(*symbol))
        .collect()
}

fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    left.into_iter()
        .zip(right)
        .map(|(left, right)| {
            let delta = left - right;
            delta * delta
        })
        .sum()
}

#[cfg(test)]
#[path = "diff_tests.rs"]
mod tests;
