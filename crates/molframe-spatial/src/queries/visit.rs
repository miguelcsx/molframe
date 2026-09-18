//! Streaming fixed-radius query execution.

use super::reduce::PairQuery;
use super::search::validate_cutoff;
use crate::brute;
use crate::{CellList, KdTree, NeighborPair, SpatialBackend, SpatialError};

/// Visits fixed-radius pairs without requiring a result vector.
///
/// The callback receives the same unordered pair set as the materialising
/// unsorted query. Every backend streams matches directly. A one-shot neighbour
/// list query uses its cell index without retaining a reusable pair cache.
/// This is intended for reductions that do not consume pair order.
///
/// # Errors
///
/// Returns the same validation and backend errors as the unsorted pair query.
pub fn for_each_pairs_within_unsorted<F>(query: &PairQuery<'_>, emit: F) -> Result<(), SpatialError>
where
    F: FnMut(NeighborPair),
{
    let PairQuery {
        positions,
        left,
        right,
        cutoff,
        ..
    } = query;
    let cutoff = *cutoff;
    validate_cutoff(cutoff)?;
    let indices = super::indices::IndexWorkspace::new(left, right, positions.len(), query.context)?;
    visit_indices(query, indices.left(), indices.right(), emit)
}

// Reuse validation and index buffers when the reduction executor falls back
// to sequential visitation instead of expanding the selections a second time.
pub(super) fn visit_indices(
    query: &PairQuery<'_>,
    left_indices: &[u32],
    right_indices: &[u32],
    mut emit: impl FnMut(NeighborPair),
) -> Result<(), SpatialError> {
    let &PairQuery {
        positions,
        cutoff,
        options,
        periodic,
        context,
        ..
    } = query;
    let plan = options.plan(left_indices.len(), right_indices.len(), cutoff)?;

    let same_selection = left_indices == right_indices;
    if matches!(
        plan.backend,
        SpatialBackend::CellList | SpatialBackend::NeighborList
    ) {
        let index = CellList::build_in(
            positions,
            right_indices,
            cutoff,
            periodic.copied(),
            options.cell_grid,
            context,
        )?;
        return if same_selection {
            crate::backends::cell::for_each_pairs_same_selection_unordered(
                &index,
                left_indices,
                cutoff,
                emit,
            )
        } else {
            index.for_each_candidate(left_indices, cutoff, |left_atom, right_atom, squared| {
                emit_unique_cross_candidate(
                    left_atom,
                    right_atom,
                    squared,
                    left_indices,
                    right_indices,
                    &mut emit,
                );
            })
        };
    }
    if plan.backend == SpatialBackend::BruteForce {
        let cutoff_squared = cutoff * cutoff;
        if same_selection {
            brute::for_each_pair_same_selection(
                positions,
                left_indices,
                cutoff_squared,
                periodic,
                emit,
            );
        } else {
            brute::for_each_candidate(
                positions,
                left_indices,
                right_indices,
                cutoff_squared,
                periodic,
                |left_atom, right_atom, squared| {
                    emit_unique_cross_candidate(
                        left_atom,
                        right_atom,
                        squared,
                        left_indices,
                        right_indices,
                        &mut emit,
                    );
                },
            );
        }
        return Ok(());
    }

    let index = KdTree::build_in(
        positions,
        right_indices,
        periodic.copied(),
        options.kd_periodic,
        context,
    )?;
    index.for_each_candidate(
        left_indices,
        cutoff,
        context,
        |left_atom, right_atom, squared| {
            emit_unique_cross_candidate(
                left_atom,
                right_atom,
                squared,
                left_indices,
                right_indices,
                &mut emit,
            );
        },
    )
}

/// Emits the canonical orientation once when cross-selections overlap.
///
/// Duplicate orientations exist only when both endpoints belong to both
/// selections. Since validated index slices are ascending, two binary searches
/// remove the reverse orientation without pair-proportional storage.
pub(super) fn emit_unique_cross_candidate(
    left_atom: u32,
    right_atom: u32,
    squared: f32,
    left: &[u32],
    right: &[u32],
    emit: &mut impl FnMut(NeighborPair),
) {
    let reverse_is_duplicate = left_atom > right_atom
        && right.binary_search(&left_atom).is_ok()
        && left.binary_search(&right_atom).is_ok();
    if !reverse_is_duplicate {
        emit(NeighborPair::new(left_atom, right_atom, squared));
    }
}
