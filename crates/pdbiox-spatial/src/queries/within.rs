//! `within` selection reduction: which query atoms lie near a target set.

use super::visit::for_each_pairs_within_unsorted;
use crate::{
    NeighborPair, PairQuery, PeriodicBox, SpatialBackend, SpatialError, SpatialSearchOptions,
};
use pdbiox_core::{ExecutionContext, selection::AtomSelection};

/// Selects query atoms within `cutoff` of any target atom.
///
/// Target atoms are included when they also belong to `query`, matching the
/// conventional meaning of `within` rather than `around`.
///
/// # Errors
///
/// Returns the same errors as [`crate::pairs_within`].
pub fn within(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    target: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<AtomSelection, SpatialError> {
    within_with_options(
        positions,
        query,
        target,
        cutoff,
        SpatialSearchOptions::with_backend(backend),
        periodic,
        context,
    )
}

/// Selects query atoms under a complete, inspectable planning profile.
///
/// # Errors
///
/// Returns the same errors as [`crate::pairs_within_with_options`].
pub fn within_with_options(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    target: &AtomSelection,
    cutoff: f32,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<AtomSelection, SpatialError> {
    let mut matched = vec![false; positions.len()];
    let mut reduction_error = None;
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions,
            left: query,
            right: target,
            cutoff,
            options,
            periodic,
            context,
        },
        |pair| {
            if reduction_error.is_none()
                && let Err(error) = mark_pair_matches(&mut matched, query, target, pair)
            {
                reduction_error = Some(error);
            }
        },
    )?;
    if let Some(error) = reduction_error {
        return Err(error);
    }

    collect_matches(query, target, &matched)
}

/// Marks query endpoints participating in query-target neighbour pairs.
///
/// The dense bit vector removes the `O(M log M)` sort/dedup stage previously
/// required by `within`.
fn mark_pair_matches(
    matched: &mut [bool],
    query: &AtomSelection,
    target: &AtomSelection,
    pair: NeighborPair,
) -> Result<(), SpatialError> {
    if query.contains(pair.first) && target.contains(pair.second) {
        mark_atom(matched, pair.first)?;
    }
    if query.contains(pair.second) && target.contains(pair.first) {
        mark_atom(matched, pair.second)?;
    }
    Ok(())
}

/// Collects marked atoms in query-selection order.
///
/// Query-target overlap is included without requiring self-pairs.
fn collect_matches(
    query: &AtomSelection,
    target: &AtomSelection,
    matched: &[bool],
) -> Result<AtomSelection, SpatialError> {
    let mut selected = Vec::new();

    for atom in query {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        let spatial_match = matched
            .get(index)
            .copied()
            .ok_or(SpatialError::AtomOutOfBounds(atom))?;

        if spatial_match || target.contains(atom) {
            selected.push(atom);
        }
    }

    Ok(AtomSelection::from_sorted(selected))
}

/// Marks one atom in a dense match bitmap.
///
/// # Errors
///
/// Returns an out-of-bounds error for an invalid pair endpoint.
fn mark_atom(matched: &mut [bool], atom: u32) -> Result<(), SpatialError> {
    let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
    let Some(slot) = matched.get_mut(index) else {
        return Err(SpatialError::AtomOutOfBounds(atom));
    };

    *slot = true;
    Ok(())
}
