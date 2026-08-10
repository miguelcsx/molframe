//! One-pass projection from structure rows to graph nodes.

use super::{GraphError, NodeLevel};
use pdbiox_core::structure::AtomRef;
use pdbiox_core::{AtomIndex, ResidueIndex, Structure};
use std::ops::Range;

use crate::numeric::f64_to_f32;

const RESIDUE_ATOMS_FEATURE: &str = "residue_atoms";

pub(super) struct Nodes {
    pub(super) positions: Vec<[f32; 3]>,
    pub(super) position_valid: Vec<bool>,
    pub(super) atom_ranges: Vec<Range<u32>>,
    pub(super) atom_to_node: Vec<u32>,
}

pub(super) fn project(structure: &Structure, level: NodeLevel) -> Result<Nodes, GraphError> {
    match level {
        NodeLevel::Atoms => atoms(structure),
        NodeLevel::Residues => residues(structure),
    }
}

fn atoms(structure: &Structure) -> Result<Nodes, GraphError> {
    let count = structure.atom_count();
    if structure.positions().len() != count as usize {
        return Err(GraphError::RaggedCoordinates);
    }
    let mut ranges = Vec::with_capacity(count as usize);
    let mut atom_to_node = Vec::with_capacity(count as usize);
    for atom in 0..count {
        let end = atom.checked_add(1).ok_or(GraphError::IndexOverflow)?;
        ranges.push(atom..end);
        atom_to_node.push(atom);
    }
    let positions = structure.positions().to_vec();
    let position_valid = positions.iter().map(|position| finite(*position)).collect();
    Ok(Nodes {
        positions,
        position_valid,
        atom_ranges: ranges,
        atom_to_node,
    })
}

fn residues(structure: &Structure) -> Result<Nodes, GraphError> {
    let mut positions = Vec::with_capacity(structure.residue_count());
    let mut position_valid = Vec::with_capacity(structure.residue_count());
    let mut atom_ranges = Vec::with_capacity(structure.residue_count());
    let mut atom_to_node = vec![u32::MAX; structure.atom_count() as usize];
    for residue in 0..structure.residue_count() {
        let residue_index = u32::try_from(residue).map_err(|_| GraphError::IndexOverflow)?;
        let Some(range) = structure
            .data()
            .topology
            .residues
            .atoms(ResidueIndex::new(residue_index))
        else {
            return Err(GraphError::MissingFeature {
                feature: RESIDUE_ATOMS_FEATURE,
                row: residue,
            });
        };
        let (position, valid) = centroid(structure, range.clone());
        positions.push(position);
        position_valid.push(valid);
        for atom in range.clone() {
            let Some(target) = atom_to_node.get_mut(atom as usize) else {
                return Err(GraphError::IndexOverflow);
            };
            *target = residue_index;
        }
        atom_ranges.push(range);
    }
    Ok(Nodes {
        positions,
        position_valid,
        atom_ranges,
        atom_to_node,
    })
}

fn centroid(structure: &Structure, atoms: Range<u32>) -> ([f32; 3], bool) {
    let mut sum = [0.0f64; 3];
    let mut count = 0u32;
    for atom in atoms {
        let Some(position) = structure
            .data()
            .atom(AtomIndex::new(atom))
            .and_then(AtomRef::position)
        else {
            continue;
        };
        if !finite(position) {
            continue;
        }
        for axis in 0..3 {
            sum[axis] += f64::from(position[axis]);
        }
        count += 1;
    }
    if count == 0 {
        return ([0.0; 3], false);
    }
    let divisor = f64::from(count);
    (sum.map(|value| f64_to_f32(value / divisor)), true)
}

fn finite(position: [f32; 3]) -> bool {
    position.iter().all(|value| value.is_finite())
}
