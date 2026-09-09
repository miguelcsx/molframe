//! Direct fixed-radius comparison.

use crate::{NeighborPair, PeriodicBox};

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
    pairs_with_mode::<false>(positions, left, right, cutoff_squared, periodic, true)
}

/// Finds all unique pairs from one selection without evaluating both orders.
pub(crate) fn pairs_same_selection(
    positions: &[[f32; 3]],
    selection: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
) -> Vec<NeighborPair> {
    pairs_with_mode::<true>(
        positions,
        selection,
        selection,
        cutoff_squared,
        periodic,
        true,
    )
}

/// Finds all unique pairs from one selection without sorting the result.
pub(crate) fn pairs_same_selection_unordered(
    positions: &[[f32; 3]],
    selection: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
) -> Vec<NeighborPair> {
    pairs_with_mode::<true>(
        positions,
        selection,
        selection,
        cutoff_squared,
        periodic,
        false,
    )
}

/// Visits oriented cross-selection candidates without allocating.
pub(crate) fn for_each_candidate(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
    mut emit: impl FnMut(u32, u32, f32),
) {
    for_each_candidate_with_mode::<false>(
        positions,
        left,
        right,
        cutoff_squared,
        periodic,
        &mut emit,
    );
}

/// Visits each same-selection unordered pair once without allocation.
pub(crate) fn for_each_pair_same_selection(
    positions: &[[f32; 3]],
    selection: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
    mut emit: impl FnMut(NeighborPair),
) {
    for_each_candidate_with_mode::<true>(
        positions,
        selection,
        selection,
        cutoff_squared,
        periodic,
        &mut |left_atom, right_atom, squared| {
            emit(NeighborPair::new(left_atom, right_atom, squared));
        },
    );
}

/// Dispatches direct comparison with compile-time pair-emission semantics.
fn pairs_with_mode<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
    sort_result: bool,
) -> Vec<NeighborPair> {
    let mut found = Vec::with_capacity(left.len());
    for_each_candidate_with_mode::<UNIQUE>(
        positions,
        left,
        right,
        cutoff_squared,
        periodic,
        &mut |left_atom, right_atom, squared| {
            found.push(NeighborPair::new(left_atom, right_atom, squared));
        },
    );
    if sort_result {
        canonicalise(&mut found);
    }
    found
}

/// Dispatches direct comparison without retaining emitted pairs.
fn for_each_candidate_with_mode<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    match periodic {
        Some(periodic) => {
            for_each_pair_periodic::<UNIQUE>(
                positions,
                left,
                right,
                cutoff_squared,
                periodic,
                emit,
            );
        }
        None => for_each_pair_simd::<UNIQUE>(positions, left, right, cutoff_squared, emit),
    }
}

/// Performs direct periodic pair comparison.
///
/// Runtime is `O(L * R)` with `O(1)` auxiliary space.
fn for_each_pair_periodic<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    periodic: &PeriodicBox,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    for &left_atom in left {
        let Some(left_position) = valid_position(positions, left_atom) else {
            continue;
        };
        let candidates = right_for_left::<UNIQUE>(right, left_atom);

        append_periodic_pairs::<UNIQUE>(
            positions,
            left_atom,
            left_position,
            candidates,
            cutoff_squared,
            periodic,
            emit,
        );
    }
}

/// Compares one valid left atom against all periodic right-side atoms.
///
/// Runtime is `O(R)` and no allocation occurs.
fn append_periodic_pairs<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left_atom: u32,
    left_position: [f32; 3],
    right: &[u32],
    cutoff_squared: f32,
    periodic: &PeriodicBox,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    for &right_atom in right {
        if !allowed_pair::<UNIQUE>(left_atom, right_atom) {
            continue;
        }

        let Some(right_position) = valid_position(positions, right_atom) else {
            continue;
        };

        let squared = periodic.distance_squared(left_position, right_position);

        if squared <= cutoff_squared {
            emit(left_atom, right_atom, squared);
        }
    }
}

/// Performs non-periodic direct comparison using four-wide SIMD.
///
/// Runtime is `O(L * R / 4)` SIMD batches plus scalar remainder work.
/// Auxiliary space is `O(1)` outside the caller-provided reducer.
fn for_each_pair_simd<const UNIQUE: bool>(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff_squared: f32,
    emit: &mut impl FnMut(u32, u32, f32),
) {
    for &left_atom in left {
        let Some(left_position) = valid_position(positions, left_atom) else {
            continue;
        };
        let candidates = right_for_left::<UNIQUE>(right, left_atom);

        super::brute_simd::append_pairs::<UNIQUE>(
            positions,
            left_atom,
            left_position,
            candidates,
            cutoff_squared,
            emit,
        );
    }
}

/// Restricts a same-selection query to the ascending suffix with larger atom
/// IDs. Normal cross-selection searches retain the complete target slice.
#[inline]
fn right_for_left<const UNIQUE: bool>(right: &[u32], left_atom: u32) -> &[u32] {
    if UNIQUE {
        let first = right.partition_point(|&right_atom| right_atom <= left_atom);
        &right[first..]
    } else {
        right
    }
}

/// Returns whether a candidate satisfies the selected pair-emission contract.
#[inline]
pub(super) fn allowed_pair<const UNIQUE: bool>(left_atom: u32, right_atom: u32) -> bool {
    if UNIQUE {
        left_atom < right_atom
    } else {
        left_atom != right_atom
    }
}

/// Returns a finite atom position when `atom` is in bounds.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
pub(super) fn valid_position(positions: &[[f32; 3]], atom: u32) -> Option<[f32; 3]> {
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
pub(super) fn euclidean_distance_squared(left: [f32; 3], right: [f32; 3]) -> f32 {
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
