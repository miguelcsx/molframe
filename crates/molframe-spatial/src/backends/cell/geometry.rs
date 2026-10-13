//! Cell-grid bounds and memory-bounded geometry.

use crate::brute::finite;
use crate::numeric::floor_usize;
use crate::{CellGridOptions, SpatialError};

#[derive(Clone, Copy)]
pub(super) struct Bounds {
    pub(super) min: [f32; 3],
    pub(super) max: [f32; 3],
    pub(super) finite_count: usize,
}

pub(super) fn finite_bounds(positions: &[[f32; 3]], targets: &[u32]) -> Option<Bounds> {
    let mut min = [f32::INFINITY; 3];
    let mut max = [f32::NEG_INFINITY; 3];
    let mut finite_count = 0usize;

    for &target in targets {
        let index = usize::try_from(target).ok()?;
        let position = positions.get(index).copied()?;
        if !finite(position) {
            continue;
        }
        finite_count += 1;
        update_bounds(&mut min, &mut max, position);
    }

    (finite_count != 0).then_some(Bounds {
        min,
        max,
        finite_count,
    })
}

#[inline]
fn update_bounds(min: &mut [f32; 3], max: &mut [f32; 3], position: [f32; 3]) {
    for axis in 0..3 {
        min[axis] = min[axis].min(position[axis]);
        max[axis] = max[axis].max(position[axis]);
    }
}

#[cfg(test)]
pub(super) fn grid_geometry(
    min: [f32; 3],
    max: [f32; 3],
    cutoff: f32,
    options: CellGridOptions,
) -> Result<(f32, [usize; 3], usize), SpatialError> {
    let mut edge = cutoff;
    let mut dims = dimensions(min, max, edge)?;
    let mut count = cell_count(dims)?;

    while count > options.maximum_cell_count {
        edge *= options.edge_growth_factor;
        dims = dimensions(min, max, edge)?;
        count = cell_count(dims)?;
    }

    Ok((edge, dims, count))
}

fn dimensions(min: [f32; 3], max: [f32; 3], edge: f32) -> Result<[usize; 3], SpatialError> {
    Ok([
        dimension(min[0], max[0], edge)?,
        dimension(min[1], max[1], edge)?,
        dimension(min[2], max[2], edge)?,
    ])
}

fn dimension(min: f32, max: f32, edge: f32) -> Result<usize, SpatialError> {
    let cells = floor_usize((f64::from(max) - f64::from(min)) / f64::from(edge))
        .ok_or(SpatialError::NumericRangeExceeded)?;
    cells
        .checked_add(1)
        .map(|dimension| dimension.max(1))
        .ok_or(SpatialError::NumericRangeExceeded)
}

#[inline]
fn cell_count(dims: [usize; 3]) -> Result<usize, SpatialError> {
    dims[0]
        .checked_mul(dims[1])
        .and_then(|plane| plane.checked_mul(dims[2]))
        .ok_or(SpatialError::NumericRangeExceeded)
}

/// Cells a dense grid may hold per indexed atom, plus a small floor.
///
/// A dense grid stores one offset per cell, so a sparse or huge-extent
/// structure would pay for empty space. Beyond this ratio the index keeps only
/// its occupied cells.
fn dense_cell_allowance(atoms: usize, options: CellGridOptions) -> usize {
    atoms
        .saturating_mul(8)
        .saturating_add(4096)
        .min(options.maximum_cell_count)
}

/// How the cells of a Cartesian grid are stored.
#[derive(Clone, Copy)]
pub(super) enum Layout {
    /// One offset slot per cell of the bounding box.
    Dense {
        edge: f32,
        dims: [usize; 3],
        cells: usize,
    },
    /// Only occupied cells, found by binary search over their linear ids.
    Sparse { edge: f32, dims: [usize; 3] },
}

/// Chooses a memory-bounded layout: dense when the bounding box is small
/// relative to the atom count, sparse otherwise. Both use cells at least
/// `cutoff` wide, so they answer every query identically.
pub(super) fn choose_layout(
    bounds: &Bounds,
    cutoff: f32,
    options: CellGridOptions,
) -> Result<Layout, SpatialError> {
    if cutoff == 0.0 {
        return Ok(Layout::Dense {
            edge: 1.0,
            dims: [1; 3],
            cells: 1,
        });
    }
    if let Ok(dims) = dimensions(bounds.min, bounds.max, cutoff)
        && let Ok(cells) = cell_count(dims)
        && cells <= dense_cell_allowance(bounds.finite_count, options)
    {
        return Ok(Layout::Dense {
            edge: cutoff,
            dims,
            cells,
        });
    }
    // Cells keep their cutoff edge whatever the extent: a wide-id sparse grid
    // never coarsens, so a far outlier cannot swallow the cluster it sits beside.
    let mut dims = [1_usize; 3];
    for (axis, dim) in dims.iter_mut().enumerate() {
        let span = (f64::from(bounds.max[axis]) - f64::from(bounds.min[axis])) / f64::from(cutoff);
        // Two spare cells absorb the rounding difference between this
        // double-precision span and the single-precision cell lookup.
        *dim = floor_usize(span)
            .and_then(|cells| cells.checked_add(2))
            .ok_or(SpatialError::NumericRangeExceeded)?;
    }
    Ok(Layout::Sparse { edge: cutoff, dims })
}
