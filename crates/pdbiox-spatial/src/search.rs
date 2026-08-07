//! Backend selection and public fixed-radius operations.

use crate::brute;
use crate::{
    CellList, KdTree, NeighborList, NeighborPair, PeriodicBox, SpatialBackend, SpatialError,
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
    if !cutoff.is_finite() || cutoff < 0.0 {
        return Err(SpatialError::InvalidCutoff);
    }
    let left = checked_indices(left, positions.len())?;
    let right = checked_indices(right, positions.len())?;
    let chosen = match backend {
        SpatialBackend::Auto => choose(left.len(), right.len(), periodic.is_some()),
        other => other,
    };
    let cutoff_squared = cutoff * cutoff;
    let pairs = match chosen {
        SpatialBackend::BruteForce => {
            brute::pairs(positions, &left, &right, cutoff_squared, periodic)
        }
        SpatialBackend::CellList => {
            CellList::build(positions, &right, cutoff, periodic)?.pairs(&left, cutoff)?
        }
        SpatialBackend::KdTree => {
            KdTree::build(positions, &right, periodic)?.pairs(&left, cutoff)?
        }
        SpatialBackend::NeighborList => NeighborList::build(
            positions,
            &left,
            &right,
            cutoff,
            default_skin(cutoff),
            periodic,
            CoordinateGeneration::INITIAL,
        )?
        .pairs(positions, cutoff, periodic)?,
        SpatialBackend::Auto => Vec::new(),
    };
    Ok(pairs)
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
    let pairs = pairs_within(positions, query, target, cutoff, backend, periodic)?;
    let mut selected: Vec<u32> = query.intersect(target).iter().collect();
    for pair in pairs {
        if query.contains(pair.first) && target.contains(pair.second) {
            selected.push(pair.first);
        }
        if query.contains(pair.second) && target.contains(pair.first) {
            selected.push(pair.second);
        }
    }
    selected.sort_unstable();
    selected.dedup();
    Ok(AtomSelection::from_sorted(selected))
}

fn checked_indices(
    selection: &AtomSelection,
    position_count: usize,
) -> Result<Vec<u32>, SpatialError> {
    let mut indices = Vec::with_capacity(selection.len() as usize);
    for atom in selection {
        if atom as usize >= position_count {
            return Err(SpatialError::AtomOutOfBounds(atom));
        }
        indices.push(atom);
    }
    Ok(indices)
}

fn choose(left: usize, right: usize, periodic: bool) -> SpatialBackend {
    if periodic || left.saturating_mul(right) <= 250_000 {
        SpatialBackend::BruteForce
    } else if right >= 20_000 && left < right / 8 {
        SpatialBackend::KdTree
    } else {
        SpatialBackend::CellList
    }
}

fn default_skin(cutoff: f32) -> f32 {
    (cutoff * 0.2).max(0.5)
}

#[cfg(test)]
#[path = "search_tests.rs"]
mod tests;
