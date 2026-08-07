//! Direct fixed-radius comparison.

use crate::{NeighborPair, PeriodicBox};
use wide::f32x4;

pub(crate) fn pairs(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
) -> Vec<NeighborPair> {
    if periodic.is_none() {
        return pairs_simd(positions, left, right, cutoff_squared);
    }
    pairs_periodic(positions, left, right, cutoff_squared, periodic)
}

fn pairs_periodic(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
) -> Vec<NeighborPair> {
    let mut found = Vec::new();
    for left_atom in left {
        let Some(left_position) = positions.get(*left_atom as usize).copied() else {
            continue;
        };
        if !finite(left_position) {
            continue;
        }
        for right_atom in right {
            if left_atom == right_atom {
                continue;
            }
            let Some(right_position) = positions.get(*right_atom as usize).copied() else {
                continue;
            };
            if !finite(right_position) {
                continue;
            }
            let distance_squared = distance_squared(left_position, right_position, periodic);
            if distance_squared <= cutoff_squared {
                found.push(NeighborPair::new(*left_atom, *right_atom, distance_squared));
            }
        }
    }
    canonicalise(&mut found);
    found
}

fn pairs_simd(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
) -> Vec<NeighborPair> {
    let mut found = Vec::new();
    for left_atom in left {
        let Some(left_position) = positions.get(*left_atom as usize).copied() else {
            continue;
        };
        if !finite(left_position) {
            continue;
        }
        let mut lanes = right.chunks_exact(4);
        for atoms in &mut lanes {
            compare_four(
                positions,
                *left_atom,
                left_position,
                atoms,
                cutoff_squared,
                &mut found,
            );
        }
        for right_atom in lanes.remainder() {
            compare_one(
                positions,
                *left_atom,
                left_position,
                *right_atom,
                cutoff_squared,
                &mut found,
            );
        }
    }
    canonicalise(&mut found);
    found
}

fn compare_four(
    positions: &[[f32; 3]],
    left_atom: u32,
    left: [f32; 3],
    atoms: &[u32],
    cutoff_squared: f32,
    found: &mut Vec<NeighborPair>,
) {
    let Some(atom_ids) = <&[u32; 4]>::try_from(atoms).ok() else {
        return;
    };
    let right = atom_ids.map(|atom| match positions.get(atom as usize).copied() {
        Some(position) => position,
        None => [f32::NAN; 3],
    });
    let dx = f32x4::from(right.map(|position| position[0])) - f32x4::splat(left[0]);
    let dy = f32x4::from(right.map(|position| position[1])) - f32x4::splat(left[1]);
    let dz = f32x4::from(right.map(|position| position[2])) - f32x4::splat(left[2]);
    let squared = (dx * dx + dy * dy + dz * dz).to_array();
    for lane in 0..4 {
        let right_atom = atom_ids[lane];
        let distance_squared = squared[lane];
        if left_atom != right_atom && distance_squared <= cutoff_squared {
            found.push(NeighborPair::new(left_atom, right_atom, distance_squared));
        }
    }
}

fn compare_one(
    positions: &[[f32; 3]],
    left_atom: u32,
    left: [f32; 3],
    right_atom: u32,
    cutoff_squared: f32,
    found: &mut Vec<NeighborPair>,
) {
    if left_atom == right_atom {
        return;
    }
    let Some(right) = positions.get(right_atom as usize).copied() else {
        return;
    };
    if !finite(right) {
        return;
    }
    let squared = distance_squared(left, right, None);
    if squared <= cutoff_squared {
        found.push(NeighborPair::new(left_atom, right_atom, squared));
    }
}

pub(crate) fn distance_squared(
    left: [f32; 3],
    right: [f32; 3],
    periodic: Option<&PeriodicBox>,
) -> f32 {
    match periodic {
        Some(periodic) => periodic.distance_squared(left, right),
        None => left
            .iter()
            .zip(right)
            .map(|(left, right)| {
                let delta = left - right;
                delta * delta
            })
            .sum(),
    }
}

pub(crate) fn finite(position: [f32; 3]) -> bool {
    position.iter().all(|component| component.is_finite())
}

pub(crate) fn canonicalise(found: &mut Vec<NeighborPair>) {
    found.sort_unstable_by_key(|pair| (pair.first, pair.second));
    found.dedup_by(|left, right| left.first == right.first && left.second == right.second);
}

#[cfg(test)]
#[path = "brute_tests.rs"]
mod tests;
