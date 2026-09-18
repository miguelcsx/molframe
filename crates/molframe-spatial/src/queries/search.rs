//! Backend selection and public fixed-radius operations.

use crate::brute;
use crate::{
    CellList, KdTree, NeighborList, NeighborPair, PeriodicBox, SpatialBackend, SpatialError,
    SpatialSearchOptions,
};
use molframe_core::selection::AtomSelection;
use molframe_core::{CoordinateGeneration, ExecutionContext};
use std::borrow::Cow;

#[derive(Clone, Copy)]
enum PairOrder {
    Sorted,
    Unsorted,
}

/// How a pair enumeration should be produced.
#[derive(Clone, Copy)]
struct PairRequest {
    options: SpatialSearchOptions,
    order: PairOrder,
}

impl PairRequest {
    const fn sorted(options: SpatialSearchOptions) -> Self {
        Self {
            options,
            order: PairOrder::Sorted,
        }
    }
    const fn unsorted(options: SpatialSearchOptions) -> Self {
        Self {
            options,
            order: PairOrder::Unsorted,
        }
    }
}

#[derive(Clone, Copy)]
struct PairDispatch {
    plan: crate::SpatialPlan,
    options: SpatialSearchOptions,
    order: PairOrder,
}

/// Enumerates unique unordered pairs no further apart than `cutoff`.
///
/// Results are sorted by `(first, second)` for every backend. Self-pairs are
/// excluded, including when the two selections overlap.
///
/// # Errors
///
/// Returns an error for an invalid cutoff or an out-of-range atom index.
pub fn pairs_within(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_within_with_options(
        positions,
        left,
        right,
        cutoff,
        SpatialSearchOptions::with_backend(backend),
        periodic,
        context,
    )
}

/// Enumerates unique same-selection pairs without sorting the returned vector.
///
/// The deterministic ordering contract of [`pairs_within`] remains unchanged.
/// This variant is intended for reductions such as histograms and counts that
/// do not consume pair order. For different left and right selections it keeps
/// the canonical implementation so overlapping selections remain duplicate-free.
///
/// # Errors
///
/// Returns the same validation and backend errors as [`pairs_within`].
pub fn pairs_within_unsorted(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_within_unsorted_with_options(
        positions,
        left,
        right,
        cutoff,
        SpatialSearchOptions::with_backend(backend),
        periodic,
        context,
    )
}

/// Enumerates pairs under a complete, inspectable planning profile.
///
/// # Errors
///
/// Returns an error for invalid options, cutoff or atom indices.
pub fn pairs_within_with_options(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_with_order(
        positions,
        left,
        right,
        cutoff,
        periodic,
        PairRequest::sorted(options),
        context,
    )
}

/// Executes a fixed-radius query without sorting identical-selection results.
///
/// # Errors
///
/// Returns the same validation and backend errors as [`pairs_within`].
pub fn pairs_within_unsorted_with_options(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_with_order(
        positions,
        left,
        right,
        cutoff,
        periodic,
        PairRequest::unsorted(options),
        context,
    )
}

fn pairs_with_order(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    periodic: Option<&PeriodicBox>,
    request: PairRequest,
    context: &ExecutionContext,
) -> Result<Vec<NeighborPair>, SpatialError> {
    validate_cutoff(cutoff)?;

    let indices = super::indices::IndexWorkspace::new(left, right, positions.len(), context)?;
    let left_indices = indices.left();
    let right_indices = indices.right();
    let plan = request.options.plan(
        left_indices.len(),
        right_indices.len(),
        periodic.is_some(),
        cutoff,
    )?;

    dispatch_pairs(
        positions,
        left_indices,
        right_indices,
        cutoff,
        PairDispatch {
            plan,
            options: request.options,
            order: request.order,
        },
        periodic,
        context,
    )
}

/// Dispatches already-validated index slices to one spatial backend.
///
/// Every plan contains a concrete backend selected before dispatch.
fn dispatch_pairs(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff: f32,
    dispatch: PairDispatch,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<Vec<NeighborPair>, SpatialError> {
    let same_selection = left == right;

    match dispatch.plan.backend {
        SpatialBackend::CellList => {
            let index = CellList::build_in(
                positions,
                right,
                cutoff,
                periodic,
                dispatch.options.cell_grid,
                context,
            )?;
            if same_selection {
                let sorted = matches!(dispatch.order, PairOrder::Sorted);
                if context.worker_budget() > 1 {
                    crate::backends::cell::pairs_same_selection_parallel(
                        &index, left, cutoff, context, sorted,
                    )
                } else if sorted {
                    crate::backends::cell::pairs_same_selection(&index, left, cutoff)
                } else {
                    crate::backends::cell::pairs_same_selection_unordered(&index, left, cutoff)
                }
            } else {
                index.pairs(left, cutoff)
            }
        }
        SpatialBackend::KdTree => KdTree::build_in(
            positions,
            right,
            periodic,
            dispatch.options.kd_periodic,
            context,
        )?
        .pairs(left, cutoff),
        SpatialBackend::NeighborList => neighbor_list_pairs(
            positions,
            left,
            right,
            cutoff,
            dispatch.plan.neighbor_skin,
            dispatch.options,
            periodic,
        ),
        SpatialBackend::BruteForce => {
            let cutoff_squared = cutoff * cutoff;
            if same_selection {
                if matches!(dispatch.order, PairOrder::Sorted) {
                    Ok(brute::pairs_same_selection(
                        positions,
                        left,
                        cutoff_squared,
                        periodic,
                    ))
                } else {
                    Ok(brute::pairs_same_selection_unordered(
                        positions,
                        left,
                        cutoff_squared,
                        periodic,
                    ))
                }
            } else {
                Ok(brute::pairs(
                    positions,
                    left,
                    right,
                    cutoff_squared,
                    periodic,
                ))
            }
        }
        SpatialBackend::Auto => Err(SpatialError::InvalidOption(
            crate::SpatialOption::PeriodicBackend,
        )),
    }
}

/// Executes the one-shot neighbour-list backend.
///
/// The helper keeps backend dispatch small while preserving the explicit
/// `NeighborList` backend semantics.
fn neighbor_list_pairs(
    positions: &[[f32; 3]],
    left: &[u32],
    right: &[u32],
    cutoff: f32,
    skin: Option<f32>,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<NeighborPair>, SpatialError> {
    let Some(skin) = skin else {
        return Err(SpatialError::InvalidOption(
            crate::SpatialOption::NeighborSkinMinimum,
        ));
    };
    let list = NeighborList::build_with_options(
        positions,
        left,
        right,
        cutoff,
        periodic,
        CoordinateGeneration::INITIAL,
        crate::NeighborListOptions {
            skin,
            cell_grid: options.cell_grid,
        },
    )?;

    list.pairs(positions, cutoff, periodic)
}

/// Validates an atom selection, borrowing existing sparse index storage.
///
/// Runtime is `O(N)`. Sparse and empty selections need no allocation; other
/// representations expand into one contiguous index buffer.
pub(super) fn checked_indices(
    selection: &AtomSelection,
    position_count: usize,
) -> Result<Cow<'_, [u32]>, SpatialError> {
    if let AtomSelection::Sparse(indices) = selection {
        for &atom in indices {
            if usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?
                >= position_count
            {
                return Err(SpatialError::AtomOutOfBounds(atom));
            }
        }
        return Ok(Cow::Borrowed(indices));
    }
    if selection.is_empty() {
        return Ok(Cow::Borrowed(&[]));
    }
    let capacity =
        usize::try_from(selection.len()).map_err(|_| SpatialError::NumericRangeExceeded)?;
    let mut indices = Vec::with_capacity(capacity);

    for atom in selection {
        let position = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        if position >= position_count {
            return Err(SpatialError::AtomOutOfBounds(atom));
        }

        indices.push(atom);
    }

    Ok(Cow::Owned(indices))
}

/// Validates a non-negative finite cutoff.
///
/// Runtime and auxiliary space are `O(1)`.
pub(super) fn validate_cutoff(cutoff: f32) -> Result<(), SpatialError> {
    if cutoff.is_finite() && cutoff >= 0.0 {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

#[cfg(test)]
#[path = "search_tests.rs"]
mod tests;
