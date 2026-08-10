//! Backend selection and public fixed-radius operations.

use crate::brute;
use crate::{
    CellList, KdTree, NeighborList, NeighborPair, PeriodicBox, SpatialBackend, SpatialError,
    SpatialSearchOptions,
};
use pdbiox_core::CoordinateGeneration;
use pdbiox_core::selection::AtomSelection;

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
    validate_cutoff(cutoff)?;

    let left = checked_indices(left, positions.len())?;
    let right = checked_indices(right, positions.len())?;
    let plan = options.plan(left.len(), right.len(), periodic.is_some(), cutoff)?;

    dispatch_pairs(positions, &left, &right, cutoff, plan, options, periodic)
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
    plan: crate::SpatialPlan,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
) -> Result<Vec<NeighborPair>, SpatialError> {
    match plan.backend {
        SpatialBackend::CellList => {
            CellList::build_with_options(positions, right, cutoff, periodic, options.cell_grid)?
                .pairs(left, cutoff)
        }
        SpatialBackend::KdTree => {
            KdTree::build_with_options(positions, right, periodic, options.kd_periodic)?
                .pairs(left, cutoff)
        }
        SpatialBackend::NeighborList => neighbor_list_pairs(
            positions,
            left,
            right,
            cutoff,
            plan.neighbor_skin,
            options,
            periodic,
        ),
        SpatialBackend::BruteForce => Ok(brute::pairs(
            positions,
            left,
            right,
            cutoff * cutoff,
            periodic,
        )),
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
