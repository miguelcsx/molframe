//! Construction of the dense and sparse cell grids from target coordinates.
//!
//! Dense grids count, prefix-sum and fill a flat table in `O(T + C)`; sparse
//! grids sort members by linear cell id in `O(T log T)` so memory stays `O(T)`.

use super::{
    CellId, Grid, Layout, choose_layout, finite, finite_bounds, linear, linear_index, raw_cell,
};
use crate::{CellGridOptions, SpatialError};

/// Builds the contiguous grid representation for finite target coordinates.
///
/// Building is `O(T + C)` for `T` targets and `C` cells and uses contiguous
/// count, offset and member arrays.
pub(super) fn build_grid(
    positions: &[[f32; 3]],
    targets: &[u32],
    cutoff: f32,
    options: CellGridOptions,
) -> Result<Option<Grid>, SpatialError> {
    let Some(bounds) = finite_bounds(positions, targets) else {
        return Ok(None);
    };
    let (edge, dims, count) = match choose_layout(&bounds, cutoff, options)? {
        Layout::Dense { edge, dims, cells } => (edge, dims, cells),
        Layout::Sparse { edge, dims } => {
            return build_sparse_grid(
                positions,
                targets,
                bounds.min,
                edge,
                dims,
                bounds.finite_count,
            )
            .map(Some);
        }
    };

    let mut counts = count_members(positions, targets, bounds.min, edge, dims, count)?;
    let mut non_empty_cells = Vec::with_capacity(bounds.finite_count.min(count));
    let mut largest_cell = 0;
    for (cell, &count) in counts.iter().enumerate() {
        largest_cell = largest_cell.max(count);
        if count != 0 {
            non_empty_cells.push(match CellId::try_from(cell) {
                Ok(id) => id,
                Err(_) => CellId::MAX,
            });
        }
    }
    let offsets = prefix_offsets(&counts)?;
    let cell_count = counts.len();
    counts.copy_from_slice(&offsets[..cell_count]);
    let members = fill_members(
        positions,
        targets,
        bounds.min,
        edge,
        dims,
        &mut counts,
        bounds.finite_count,
    )?;

    Ok(Some(Grid {
        origin: bounds.min,
        edge,
        dims,
        offsets,
        members,
        non_empty_cells,
        largest_cell,
        sparse: false,
    }))
}

/// Builds a grid that stores only its occupied cells.
///
/// Runtime is `O(T log T)` and memory is `O(T)` however large the bounding
/// box: members are sorted by linear cell id, stably so each cell keeps target
/// order exactly as the dense layout does.
fn build_sparse_grid(
    positions: &[[f32; 3]],
    targets: &[u32],
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    finite_count: usize,
) -> Result<Grid, SpatialError> {
    let mut entries: Vec<(CellId, u32)> = Vec::with_capacity(finite_count);
    for &target in targets {
        let Ok(index) = usize::try_from(target) else {
            return Err(SpatialError::NumericRangeExceeded);
        };
        let Some(position) = positions.get(index).copied() else {
            continue;
        };
        if !finite(position) {
            continue;
        }
        entries.push((
            linear(dims, raw_cell(origin, edge, dims, position)?),
            target,
        ));
    }
    entries.sort_by_key(|&(cell, _)| cell);

    let mut non_empty_cells: Vec<CellId> = Vec::with_capacity(entries.len());
    let mut offsets: Vec<usize> = Vec::with_capacity(entries.len() + 1);
    let mut members: Vec<u32> = Vec::with_capacity(entries.len());
    let mut largest_cell = 0;
    for (position, &(cell, target)) in entries.iter().enumerate() {
        if non_empty_cells.last() != Some(&cell) {
            non_empty_cells.push(match CellId::try_from(cell) {
                Ok(id) => id,
                Err(_) => CellId::MAX,
            });
            offsets.push(position);
        }
        members.push(target);
    }
    offsets.push(members.len());
    for pair in offsets.windows(2) {
        largest_cell = largest_cell.max(pair[1] - pair[0]);
    }
    drop(entries);
    non_empty_cells.shrink_to_fit();
    offsets.shrink_to_fit();

    Ok(Grid {
        origin,
        edge,
        dims,
        offsets,
        members,
        non_empty_cells,
        largest_cell,
        sparse: true,
    })
}

/// Computes one count per grid cell.
///
/// Runtime is `O(T)` with `O(C)` count storage.
fn count_members(
    positions: &[[f32; 3]],
    targets: &[u32],
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    count: usize,
) -> Result<Vec<usize>, SpatialError> {
    let mut counts = vec![0usize; count];

    for &target in targets {
        let Ok(index) = usize::try_from(target) else {
            return Err(SpatialError::NumericRangeExceeded);
        };
        let Some(position) = positions.get(index).copied() else {
            continue;
        };

        if !finite(position) {
            continue;
        }

        let cell = linear_index(dims, raw_cell(origin, edge, dims, position)?)?;

        if let Some(slot) = counts.get_mut(cell) {
            *slot += 1;
        }
    }

    Ok(counts)
}

/// Converts per-cell counts into CSR-style prefix offsets.
///
/// Runtime is `O(C)` and output space is `O(C)`.
fn prefix_offsets(counts: &[usize]) -> Result<Vec<usize>, SpatialError> {
    let capacity = counts
        .len()
        .checked_add(1)
        .ok_or(SpatialError::NumericRangeExceeded)?;
    let mut offsets = Vec::with_capacity(capacity);
    let mut total = 0usize;

    offsets.push(total);

    for &count in counts {
        total = total
            .checked_add(count)
            .ok_or(SpatialError::NumericRangeExceeded)?;
        offsets.push(total);
    }

    Ok(offsets)
}

/// Fills the contiguous CSR-style member array.
///
/// Runtime is `O(T + C)` and output space is exactly one `u32` per finite
/// target occurrence.
fn fill_members(
    positions: &[[f32; 3]],
    targets: &[u32],
    origin: [f32; 3],
    edge: f32,
    dims: [usize; 3],
    cursors: &mut [usize],
    finite_count: usize,
) -> Result<Vec<u32>, SpatialError> {
    let mut members = vec![0u32; finite_count];

    for &target in targets {
        let Ok(index) = usize::try_from(target) else {
            return Err(SpatialError::NumericRangeExceeded);
        };
        let Some(position) = positions.get(index).copied() else {
            continue;
        };

        if !finite(position) {
            continue;
        }

        let cell = linear_index(dims, raw_cell(origin, edge, dims, position)?)?;

        let Some(cursor) = cursors.get_mut(cell) else {
            continue;
        };

        let Some(slot) = members.get_mut(*cursor) else {
            continue;
        };

        *slot = target;
        *cursor += 1;
    }

    Ok(members)
}
