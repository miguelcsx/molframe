//! Lightweight automatic secondary-structure assignment.
//!
//! A sorted cell grid limits DSSP candidates to CA pairs within 9 Å. The
//! hydrogen-bond pass is therefore O(R*k) for R residues and bounded local
//! density, with one dense backbone column and one candidate table.

use crate::grid::{CellGrid, cell_for};
use molframe_core::SecondaryStructure;
use molframe_core::structure::Structure;
use std::collections::HashSet;

const CA_CUTOFF: f32 = 9.0;
const HBOND_CUTOFF: f64 = -0.5;

#[derive(Clone, Copy, Debug, Default)]
struct Backbone {
    chain: u32,
    ca: Option<[f32; 3]>,
    carbon: Option<[f32; 3]>,
    oxygen: Option<[f32; 3]>,
    nitrogen: Option<[f32; 3]>,
}

/// Assigns unknown residues using a grid-bounded DSSP-like pass.
///
/// File records are intentionally not consulted here; the caller merges this
/// fallback only into rows that remain unknown after parsing.
#[must_use]
pub fn assign_secondary_structure(structure: &Structure) -> Vec<SecondaryStructure> {
    let mut backbones = vec![Backbone::default(); structure.residue_count()];
    for chain in structure.data().chains() {
        for residue in chain.residues() {
            let Some(slot) = backbones.get_mut(residue.index().as_usize()) else {
                continue;
            };
            slot.chain = chain.index().get();
            for atom in residue.atoms() {
                let position = atom.position();
                match atom.name() {
                    Some("CA") => slot.ca = position,
                    Some("N") => slot.nitrogen = position,
                    Some("C") => slot.carbon = position,
                    Some("O") => slot.oxygen = position,
                    _ => {}
                }
            }
        }
    }

    let mut states = backbones
        .iter()
        .map(|backbone| {
            backbone
                .ca
                .map_or(SecondaryStructure::Unknown, |_| SecondaryStructure::Coil)
        })
        .collect::<Vec<_>>();
    let mut entries = Vec::new();
    for (residue, backbone) in backbones.iter().enumerate() {
        let (Some(position), Ok(residue)) = (backbone.ca, u32::try_from(residue)) else {
            continue;
        };
        let Some(cell) = cell_for(position, CA_CUTOFF) else {
            continue;
        };
        entries.push((cell, residue));
    }
    let grid = CellGrid::build(entries);
    let mut hydrogen_bonds = Vec::new();
    grid.for_each_cell(|own, neighbourhood| {
        for &left in &grid.items()[own] {
            let Some(left_backbone) = backbones.get(left as usize) else {
                continue;
            };
            for range in neighbourhood {
                for &right in &grid.items()[range.clone()] {
                    if right <= left {
                        continue;
                    }
                    let Some(right_backbone) = backbones.get(right as usize) else {
                        continue;
                    };
                    if left_backbone.chain != right_backbone.chain
                        || squared_distance(left_backbone.ca, right_backbone.ca)
                            > CA_CUTOFF * CA_CUTOFF
                    {
                        continue;
                    }
                    let (left, right) = (left as usize, right as usize);
                    if hbond_energy(left_backbone, right_backbone) < HBOND_CUTOFF {
                        hydrogen_bonds.push((left, right));
                    }
                    if hbond_energy(right_backbone, left_backbone) < HBOND_CUTOFF {
                        hydrogen_bonds.push((right, left));
                    }
                }
            }
        }
    });
    classify(&mut states, &hydrogen_bonds);
    if hydrogen_bonds.is_empty() {
        assign_zhang_skolnick(&backbones, &mut states);
    }
    states
}

fn classify(states: &mut [SecondaryStructure], bonds: &[(usize, usize)]) {
    let bond_set: HashSet<(usize, usize)> = bonds.iter().copied().collect();
    for &(left, right) in bonds {
        if left.abs_diff(right) == 4 {
            let (start, end) = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            if let Some(slice) = states.get_mut(start..=end) {
                for state in slice {
                    *state = SecondaryStructure::Helix;
                }
            }
        }
    }
    for &(left, right) in bonds {
        if left.abs_diff(right) >= 3
            && bond_set.contains(&(right, left))
            && let Some(state) = states.get_mut(left)
            && *state != SecondaryStructure::Helix
        {
            *state = SecondaryStructure::Strand;
        }
        if left.abs_diff(right) == 3 {
            let (start, end) = if left < right {
                (left, right)
            } else {
                (right, left)
            };
            if let Some(slice) = states.get_mut(start..=end) {
                for state in slice {
                    if *state == SecondaryStructure::Coil {
                        *state = SecondaryStructure::Turn;
                    }
                }
            }
        }
    }
}

fn assign_zhang_skolnick(backbones: &[Backbone], states: &mut [SecondaryStructure]) {
    for index in 0..backbones.len() {
        let Some(current) = backbones.get(index).and_then(|backbone| backbone.ca) else {
            continue;
        };
        let Some(next_three) = backbones.get(index + 3).and_then(|backbone| backbone.ca) else {
            continue;
        };
        if backbones
            .get(index + 3)
            .is_some_and(|backbone| backbone.chain == backbones[index].chain)
            && squared_distance(Some(current), Some(next_three)) < 6.5 * 6.5
            && let Some(state) = states.get_mut(index..=index + 3)
        {
            for value in state {
                *value = SecondaryStructure::Helix;
            }
        }
    }
}

fn hbond_energy(carbonyl: &Backbone, amide: &Backbone) -> f64 {
    let (Some(carbon), Some(oxygen), Some(nitrogen)) =
        (carbonyl.carbon, carbonyl.oxygen, amide.nitrogen)
    else {
        return f64::INFINITY;
    };
    let hydrogen = [
        nitrogen[0] + oxygen[0] - carbon[0],
        nitrogen[1] + oxygen[1] - carbon[1],
        nitrogen[2] + oxygen[2] - carbon[2],
    ];
    let on = distance(oxygen, nitrogen);
    let ch = distance(carbon, hydrogen);
    let oh = distance(oxygen, hydrogen);
    let cn = distance(carbon, nitrogen);
    if on <= f32::EPSILON || ch <= f32::EPSILON || oh <= f32::EPSILON || cn <= f32::EPSILON {
        return f64::INFINITY;
    }
    -27.888_f64
        * (1.0 / f64::from(on) + 1.0 / f64::from(ch) - 1.0 / f64::from(oh) - 1.0 / f64::from(cn))
}

fn squared_distance(left: Option<[f32; 3]>, right: Option<[f32; 3]>) -> f32 {
    let (Some(left), Some(right)) = (left, right) else {
        return f32::INFINITY;
    };
    left.iter().zip(right).map(|(a, b)| (a - b) * (a - b)).sum()
}

fn distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    squared_distance(Some(left), Some(right)).sqrt()
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
