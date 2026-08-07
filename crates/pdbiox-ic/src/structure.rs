//! Conversion of a structure model into an internal-coordinate forest.

use crate::{Dihedron, InternalAtom, InternalCoordinates};
use pdbiox_core::{AtomIndex, Code, Diagnostic, ModelIndex, Structure};
use std::collections::VecDeque;

/// Converts one coordinate model using a deterministic bond-graph forest.
///
/// Each connected component keeps enough Cartesian seeds to establish its
/// frame. Deeper atoms store only a bond length, angle and torsion. Missing or
/// degenerate positions remain absent/seeds rather than fabricated geometry.
///
/// # Errors
///
/// Returns `E6003` when the selected model does not exist.
pub fn internal_coordinates(
    structure: &Structure,
    model: ModelIndex,
) -> Result<InternalCoordinates, Diagnostic> {
    let positions = structure
        .model_positions(model)
        .ok_or_else(|| Diagnostic::new(Code::E6003).with_context("model", model.to_string()))?;
    let atom_count = structure.atom_count() as usize;
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    let (parents, order) = spanning_forest(atom_count, adjacency);
    let mut seeds = Vec::new();
    let mut atoms = Vec::new();
    for atom in order {
        let index = AtomIndex::new(atom as u32);
        let Some(position) = positions.get(atom).copied().filter(finite) else {
            continue;
        };
        let references = ancestor_triplet(index, &parents);
        let coordinate = references.and_then(|references| {
            let [i, j, k] = references.map(AtomIndex::as_usize);
            Dihedron::from_points(
                *positions.get(i)?,
                *positions.get(j)?,
                *positions.get(k)?,
                position,
            )
            .map(|coordinate| (references, coordinate))
        });
        match coordinate {
            Some((references, coordinate)) => atoms.push(InternalAtom {
                atom: index,
                references,
                coordinate,
            }),
            None => seeds.push((index, position)),
        }
    }
    Ok(InternalCoordinates::new(atom_count, seeds, atoms))
}

fn spanning_forest(
    atom_count: usize,
    adjacency: &pdbiox_core::BondAdjacency,
) -> (Vec<Option<AtomIndex>>, Vec<usize>) {
    let mut parents = vec![None; atom_count];
    let mut visited = vec![false; atom_count];
    let mut order = Vec::with_capacity(atom_count);
    for root in 0..atom_count {
        if visited[root] {
            continue;
        }
        visited[root] = true;
        let mut queue = VecDeque::from([AtomIndex::new(root as u32)]);
        while let Some(atom) = queue.pop_front() {
            order.push(atom.as_usize());
            for neighbor in adjacency.neighbours(atom) {
                let position = neighbor.as_usize();
                if position >= atom_count || visited[position] {
                    continue;
                }
                visited[position] = true;
                parents[position] = Some(atom);
                queue.push_back(*neighbor);
            }
        }
    }
    (parents, order)
}

fn ancestor_triplet(atom: AtomIndex, parents: &[Option<AtomIndex>]) -> Option<[AtomIndex; 3]> {
    let k = *parents.get(atom.as_usize())?;
    let j = *parents.get(k?.as_usize())?;
    let i = *parents.get(j?.as_usize())?;
    Some([i?, j?, k?])
}

fn finite(position: &[f32; 3]) -> bool {
    position.iter().all(|value| value.is_finite())
}

#[cfg(test)]
#[path = "structure_tests.rs"]
mod tests;
