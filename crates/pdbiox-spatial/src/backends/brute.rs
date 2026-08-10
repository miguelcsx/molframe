//! Direct fixed-radius comparison.

use crate::{NeighborPair, PeriodicBox};
use wide::f32x4;

/// Finds all unique pairs between `left` and `right` within `cutoff_squared`.
///
/// Non-periodic searches use four-wide SIMD. Periodic searches delegate
/// distance evaluation to [`PeriodicBox`].
///
/// Time complexity is `O(L * R)`, where `L` and `R` are the selection sizes.
/// Additional space is `O(M)` for the returned matches.
pub(crate) fn pairs(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
) -> Vec<NeighborPair> {
    match periodic {
        Some(periodic) => pairs_periodic(positions, left, right, cutoff_squared, periodic),
        None => pairs_simd(positions, left, right, cutoff_squared),
    }
}

/// Performs direct periodic pair comparison.
///
/// Runtime is `O(L * R)` and only the result vector is allocated.
fn pairs_periodic(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: &PeriodicBox,
) -> Vec<NeighborPair> {
    let mut found = Vec::new();

    for &left_atom in left {
        let Some(left_position) = valid_position(positions, left_atom) else {
            continue;
        };

        append_periodic_pairs(
            positions,
            left_atom,
            left_position,
            right,
            cutoff_squared,
            periodic,
            &mut found,
        );
    }

    canonicalise(&mut found);
    found
}

/// Compares one valid left atom against all periodic right-side atoms.
///
/// Runtime is `O(R)` and no allocation occurs except growth of `found`.
fn append_periodic_pairs(
    positions: &[[f32; 3]],
    left_atom: u32,
    left_position: [f32; 3],
    right: &[u32],
    cutoff_squared: f32,
    periodic: &PeriodicBox,
    found: &mut Vec<NeighborPair>,
) {
    for &right_atom in right {
        if left_atom == right_atom {
            continue;
        }

        let Some(right_position) = valid_position(positions, right_atom) else {
            continue;
        };

        let squared = periodic.distance_squared(left_position, right_position);

        if squared <= cutoff_squared {
            found.push(NeighborPair::new(left_atom, right_atom, squared));
        }
    }
}

/// Performs non-periodic direct comparison using four-wide SIMD.
///
/// Runtime is `O(L * R / 4)` SIMD batches plus scalar remainder work.
/// Additional space is `O(M)` for the returned matches.
fn pairs_simd(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
) -> Vec<NeighborPair> {
    let mut found = Vec::new();

    for &left_atom in left {
        let Some(left_position) = valid_position(positions, left_atom) else {
            continue;
        };

        append_simd_pairs(
            positions,
            left_atom,
            left_position,
            right,
            cutoff_squared,
            &mut found,
        );
    }

    canonicalise(&mut found);
    found
}

/// Compares one valid left atom against all non-periodic right atoms.
///
/// Four targets are processed per SIMD batch and the final remainder is
/// evaluated scalarly. No temporary heap allocation is performed.
fn append_simd_pairs(
    positions: &[[f32; 3]],
    left_atom: u32,
    left_position: [f32; 3],
    right: &[u32],
    cutoff_squared: f32,
    found: &mut Vec<NeighborPair>,
) {
    let mut batches = right.chunks_exact(4);

    for atoms in &mut batches {
        compare_four(
            positions,
            left_atom,
            left_position,
            atoms,
            cutoff_squared,
            found,
        );
    }

    for &right_atom in batches.remainder() {
        compare_one(
            positions,
            left_atom,
            left_position,
            right_atom,
            cutoff_squared,
            found,
        );
    }
}

/// Compares four target atoms with one left position using SIMD.
///
/// Missing/non-finite target positions are represented by NaN lanes, which
/// cannot satisfy the cutoff comparison.
fn compare_four(
    positions: &[[f32; 3]],
    left_atom: u32,
    left: [f32; 3],
    atoms: &[u32],
    cutoff_squared: f32,
    found: &mut Vec<NeighborPair>,
) {
    // The sentinel remains inside SIMD lanes: every comparison with it is false,
    // and only pairs whose comparison mask is true are emitted below.
    const MASKED_POSITION: [f32; 3] = [f32::NAN; 3];
    let Ok(atom_ids) = <&[u32; 4]>::try_from(atoms) else {
        return;
    };

    let right = atom_ids.map(|atom| {
        match usize::try_from(atom)
            .ok()
            .and_then(|index| positions.get(index))
            .copied()
            .filter(|position| finite(*position))
        {
            Some(position) => position,
            None => MASKED_POSITION,
        }
    });

    let dx = f32x4::from(right.map(|position| position[0])) - f32x4::splat(left[0]);
    let dy = f32x4::from(right.map(|position| position[1])) - f32x4::splat(left[1]);
    let dz = f32x4::from(right.map(|position| position[2])) - f32x4::splat(left[2]);

    let squared = (dx * dx + dy * dy + dz * dz).to_array();

    for lane in 0..4 {
        append_if_within(
            left_atom,
            atom_ids[lane],
            squared[lane],
            cutoff_squared,
            found,
        );
    }
}

/// Compares one scalar target atom against a valid left position.
///
/// Runtime and auxiliary space are `O(1)`.
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

    let Some(right) = valid_position(positions, right_atom) else {
        return;
    };

    let squared = euclidean_distance_squared(left, right);

    append_if_within(left_atom, right_atom, squared, cutoff_squared, found);
}

/// Adds one non-self pair when its squared distance satisfies the cutoff.
///
/// Runtime and additional space are `O(1)` apart from result-vector growth.
#[inline]
fn append_if_within(
    left_atom: u32,
    right_atom: u32,
    squared: f32,
    cutoff_squared: f32,
    found: &mut Vec<NeighborPair>,
) {
    if left_atom != right_atom && squared <= cutoff_squared {
        found.push(NeighborPair::new(left_atom, right_atom, squared));
    }
}

/// Returns a finite atom position when `atom` is in bounds.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn valid_position(positions: &[[f32; 3]], atom: u32) -> Option<[f32; 3]> {
    let index = usize::try_from(atom).ok()?;
    positions
        .get(index)
        .copied()
        .filter(|position| finite(*position))
}

/// Computes squared separation under optional periodic boundaries.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
pub(crate) fn distance_squared(
    left: [f32; 3],
    right: [f32; 3],
    periodic: Option<&PeriodicBox>,
) -> f32 {
    match periodic {
        Some(periodic) => periodic.distance_squared(left, right),
        None => euclidean_distance_squared(left, right),
    }
}

/// Computes ordinary three-dimensional squared Euclidean distance.
///
/// The scalar expression avoids iterator machinery in this hot path.
#[inline]
fn euclidean_distance_squared(left: [f32; 3], right: [f32; 3]) -> f32 {
    let dx = left[0] - right[0];
    let dy = left[1] - right[1];
    let dz = left[2] - right[2];

    dx * dx + dy * dy + dz * dz
}

/// Returns whether every coordinate is finite.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
pub(crate) fn finite(position: [f32; 3]) -> bool {
    position[0].is_finite() && position[1].is_finite() && position[2].is_finite()
}

/// Sorts pairs deterministically and removes duplicate atom pairs.
///
/// Runtime is `O(M log M)` for `M` matches and requires no additional heap
/// allocation beyond the existing vector.
pub(crate) fn canonicalise(found: &mut Vec<NeighborPair>) {
    if found.len() < 2 {
        return;
    }

    found.sort_unstable_by_key(|pair| (pair.first, pair.second));
    found.dedup_by(|left, right| left.first == right.first && left.second == right.second);
}

#[cfg(test)]
#[path = "brute_tests.rs"]
mod tests;
