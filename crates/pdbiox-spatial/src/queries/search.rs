//! Backend selection and public fixed-radius operations.

use crate::brute;
use crate::{
    CellList, KdTree, NeighborList, NeighborPair, PeriodicBox, SpatialBackend, SpatialError,
    SpatialSearchOptions,
};
use pdbiox_core::CoordinateGeneration;
use pdbiox_core::selection::AtomSelection;

#[derive(Clone, Copy)]
enum PairOrder {
    Sorted,
    Unsorted,
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
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_within_with_options(
        positions,
        left,
        right,
        cutoff,
        SpatialSearchOptions::with_backend(backend),
        periodic,
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
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_within_unsorted_with_options(
        positions,
        left,
        right,
        cutoff,
        SpatialSearchOptions::with_backend(backend),
        periodic,
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
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_with_order(
        positions,
        left,
        right,
        cutoff,
        options,
        periodic,
        PairOrder::Sorted,
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
) -> Result<Vec<NeighborPair>, SpatialError> {
    pairs_with_order(
        positions,
        left,
        right,
        cutoff,
        options,
        periodic,
        PairOrder::Unsorted,
    )
}

/// Visits fixed-radius same-selection pairs without requiring a result vector.
///
/// The callback receives the same unordered pair set as
/// [`pairs_within_unsorted`]. Cell-list execution streams matches directly;
/// other backends retain their existing implementation and forward its result
/// vector. This is intended for reductions that do not consume pair order.
///
/// # Errors
///
/// Returns the same validation and backend errors as [`pairs_within_unsorted`].
pub fn for_each_pairs_within_unsorted<F>(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
    mut emit: F,
) -> Result<(), SpatialError>
where
    F: FnMut(NeighborPair),
{
    validate_cutoff(cutoff)?;
    let left_indices = checked_indices(left, positions.len())?;
    let right_indices = checked_indices(right, positions.len())?;
    let plan = options.plan(
        left_indices.len(),
        right_indices.len(),
        periodic.is_some(),
        cutoff,
    )?;

    if plan.backend == SpatialBackend::CellList && left_indices == right_indices {
        let index = CellList::build_with_options(
            positions,
            &right_indices,
            cutoff,
            periodic,
            options.cell_grid,
        )?;
        return crate::backends::cell::for_each_pairs_same_selection_unordered(
            &index,
            &left_indices,
            cutoff,
            emit,
        );
    }

    let pairs = dispatch_pairs(
        positions,
        &left_indices,
        &right_indices,
        cutoff,
        PairDispatch {
            plan,
            options,
            order: PairOrder::Unsorted,
        },
        periodic,
    )?;
    for pair in pairs {
        emit(pair);
    }
    Ok(())
}

fn pairs_with_order(
    positions: &[[f32; 3]],
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
    order: PairOrder,
) -> Result<Vec<NeighborPair>, SpatialError> {
    validate_cutoff(cutoff)?;

    let left = checked_indices(left, positions.len())?;
    let right = checked_indices(right, positions.len())?;
    let plan = options.plan(left.len(), right.len(), periodic.is_some(), cutoff)?;

    dispatch_pairs(
        positions,
        &left,
        &right,
        cutoff,
        PairDispatch {
            plan,
            options,
            order,
        },
        periodic,
    )
}

/// Selects query atoms within `cutoff` of any target atom.
///
/// Target atoms are included when they also belong to `query`, matching the
/// conventional meaning of `within` rather than `around`.
///
/// # Errors
///
/// Returns the same errors as [`pairs_within`].
pub fn within(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    target: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
) -> Result<AtomSelection, SpatialError> {
    within_with_options(
        positions,
        query,
        target,
        cutoff,
        SpatialSearchOptions::with_backend(backend),
        periodic,
    )
}

/// Selects query atoms under a complete, inspectable planning profile.
///
/// # Errors
///
/// Returns the same errors as [`pairs_within_with_options`].
pub fn within_with_options(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    target: &AtomSelection,
    cutoff: f32,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<AtomSelection, SpatialError> {
    let pairs = pairs_within_with_options(positions, query, target, cutoff, options, periodic)?;

    let mut matched = vec![false; positions.len()];

    mark_matches(&mut matched, query, target, pairs)?;

    collect_matches(query, target, &matched)
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
) -> Result<Vec<NeighborPair>, SpatialError> {
    let same_selection = left == right;

    match dispatch.plan.backend {
        SpatialBackend::CellList => {
            let index = CellList::build_with_options(
                positions,
                right,
                cutoff,
                periodic,
                dispatch.options.cell_grid,
            )?;
            if same_selection {
                if matches!(dispatch.order, PairOrder::Sorted) {
                    crate::backends::cell::pairs_same_selection(&index, left, cutoff)
                } else {
                    crate::backends::cell::pairs_same_selection_unordered(&index, left, cutoff)
                }
            } else {
                index.pairs(left, cutoff)
            }
        }
        SpatialBackend::KdTree => {
            KdTree::build_with_options(positions, right, periodic, dispatch.options.kd_periodic)?
                .pairs(left, cutoff)
        }
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

/// Copies and validates an atom selection.
///
/// Runtime is `O(N)` with exactly one output allocation.
fn checked_indices(
    selection: &AtomSelection,
    position_count: usize,
) -> Result<Vec<u32>, SpatialError> {
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

    Ok(indices)
}

/// Marks query endpoints participating in query-target neighbour pairs.
///
/// The dense bit vector removes the `O(M log M)` sort/dedup stage previously
/// required by `within`.
fn mark_matches(
    matched: &mut [bool],
    query: &AtomSelection,
    target: &AtomSelection,
    pairs: Vec<NeighborPair>,
) -> Result<(), SpatialError> {
    for pair in pairs {
        if query.contains(pair.first) && target.contains(pair.second) {
            mark_atom(matched, pair.first)?;
        }

        if query.contains(pair.second) && target.contains(pair.first) {
            mark_atom(matched, pair.second)?;
        }
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

/// Validates a non-negative finite cutoff.
///
/// Runtime and auxiliary space are `O(1)`.
fn validate_cutoff(cutoff: f32) -> Result<(), SpatialError> {
    if cutoff.is_finite() && cutoff >= 0.0 {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

#[cfg(test)]
#[path = "search_tests.rs"]
mod tests;
