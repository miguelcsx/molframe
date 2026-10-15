//! Atom depth: how far each atom sits from the molecular surface.
//!
//! A deeply buried atom is far from the nearest surface point; an exposed one is
//! right on it. Given a set of surface points — from the accessible-surface
//! sampler, say — this returns each atom's distance to the closest one. The
//! surface points are bucketed on a grid and searched shell by shell outward, so
//! the nearest is found without measuring every atom against every point.
//!
//! Runs in roughly `O(atoms · shells)` once the grid is built, which is far below
//! the product of the two counts when the points are spread out.

use std::collections::HashMap;
use std::collections::hash_map::Entry;

use crate::numeric::{f64_to_f32, f64_to_i64, i64_to_f64, i128_to_i64};

/// Explicit spatial resolution for atom-depth queries.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomDepthOptions {
    /// Surface-point bucket edge length in angstroms.
    pub cell_size: f64,
}

/// Invalid atom-depth spatial configuration.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[error("atom-depth cell size must be finite and greater than zero")]
pub struct AtomDepthError;

/// One signed grid coordinate.
type CellIndex = (i64, i64, i64);

/// Bounds of occupied surface-grid cells.
#[derive(Clone, Copy)]
struct CellBounds {
    min: CellIndex,
    max: CellIndex,
}

/// Surface-point buckets and their occupied coordinate extent.
struct SurfaceGrid {
    buckets: HashMap<CellIndex, Vec<usize>>,
    bounds: CellBounds,
}

/// The distance from each atom to the nearest surface point.
///
/// Atoms are returned in input order. When there are no surface points every
/// depth is infinite, since there is no surface to measure against.
///
/// # Examples
///
/// ```
/// use molframe_surface::{AtomDepthOptions, atom_depths};
///
/// let depths = atom_depths(
///     &[[0.0, 0.0, 0.0]],
///     &[[0.0, 0.0, 5.0], [10.0, 0.0, 0.0]],
///     AtomDepthOptions { cell_size: 3.0 },
/// )?;
/// assert!((depths[0] - 5.0).abs() < 1e-4);
/// # Ok::<(), molframe_surface::AtomDepthError>(())
/// ```
///
/// # Errors
///
/// Returns [`AtomDepthError`] when `cell_size` is non-finite or not positive.
pub fn atom_depths(
    atoms: &[[f32; 3]],
    surface: &[[f32; 3]],
    options: AtomDepthOptions,
) -> Result<Vec<f32>, AtomDepthError> {
    if !options.cell_size.is_finite() || options.cell_size <= 0.0 {
        return Err(AtomDepthError);
    }
    let Some(grid) = build_grid(surface, options.cell_size) else {
        return Ok(vec![f32::INFINITY; atoms.len()]);
    };

    Ok(atoms
        .iter()
        .copied()
        .map(|atom| nearest(atom, surface, &grid, options.cell_size))
        .collect())
}

/// Builds surface-point buckets while ignoring non-finite points.
///
/// Non-finite points can never yield a finite Euclidean distance, so omitting
/// them preserves nearest-distance semantics while preventing pathological cell
/// coordinates.
fn build_grid(surface: &[[f32; 3]], cell_size: f64) -> Option<SurfaceGrid> {
    let mut buckets = HashMap::<CellIndex, Vec<usize>>::with_capacity(surface.len());
    let mut bounds = None;

    for (index, &point) in surface.iter().enumerate() {
        if !finite(point) {
            continue;
        }

        let cell = cell_of(point, cell_size);
        update_bounds(&mut bounds, cell);

        match buckets.entry(cell) {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => entry.insert(Vec::new()),
        }
        .push(index);
    }

    bounds.map(|bounds| SurfaceGrid { buckets, bounds })
}

/// Extends optional occupied-cell bounds with `cell`.
///
/// Runtime and auxiliary space are `O(1)`.
fn update_bounds(bounds: &mut Option<CellBounds>, cell: CellIndex) {
    match bounds {
        Some(bounds) => {
            bounds.min.0 = bounds.min.0.min(cell.0);
            bounds.min.1 = bounds.min.1.min(cell.1);
            bounds.min.2 = bounds.min.2.min(cell.2);

            bounds.max.0 = bounds.max.0.max(cell.0);
            bounds.max.1 = bounds.max.1.max(cell.1);
            bounds.max.2 = bounds.max.2.max(cell.2);
        }
        None => {
            *bounds = Some(CellBounds {
                min: cell,
                max: cell,
            });
        }
    }
}

/// The grid cell an point falls in.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn cell_of(point: [f32; 3], cell_size: f64) -> CellIndex {
    (
        f64_to_i64((f64::from(point[0]) / cell_size).floor()),
        f64_to_i64((f64::from(point[1]) / cell_size).floor()),
        f64_to_i64((f64::from(point[2]) / cell_size).floor()),
    )
}

/// The distance from an atom to the nearest surface point, via shell search.
///
/// Shell cells are streamed directly rather than materialised in a temporary
/// `Vec`, removing `O(r²)` allocation at every search radius.
fn nearest(atom: [f32; 3], surface: &[[f32; 3]], grid: &SurfaceGrid, cell_size: f64) -> f32 {
    if !finite(atom) {
        return f32::INFINITY;
    }

    let centre = cell_of(atom, cell_size);
    let max_ring = furthest_required_ring(centre, grid.bounds);

    let mut best = f64::INFINITY;
    let mut ring = 0i64;

    loop {
        // No point in a shell this far out can beat the best already found.
        if shell_lower_bound_squared(ring, cell_size) > best {
            break;
        }

        scan_shell(centre, ring, atom, surface, &grid.buckets, &mut best);

        if ring >= max_ring {
            break;
        }

        ring += 1;
    }

    if best.is_finite() {
        f64_to_f32(best.sqrt())
    } else {
        f32::INFINITY
    }
}

/// Scans the grid cells at Chebyshev distance exactly `ring`.
///
/// The six faces are emitted without duplicates and without a heap allocation.
/// A shell contains `O(r²)` cells.
fn scan_shell(
    centre: CellIndex,
    ring: i64,
    atom: [f32; 3],
    surface: &[[f32; 3]],
    grid: &HashMap<CellIndex, Vec<usize>>,
    best: &mut f64,
) {
    if ring == 0 {
        scan_cell(centre, atom, surface, grid, best);
        return;
    }

    let inner = ring - 1;

    for dx in -ring..=ring {
        for dy in -ring..=ring {
            scan_offset_cell(centre, (dx, dy, -ring), atom, surface, grid, best);
            scan_offset_cell(centre, (dx, dy, ring), atom, surface, grid, best);
        }
    }

    for dx in -ring..=ring {
        for dz in -inner..=inner {
            scan_offset_cell(centre, (dx, -ring, dz), atom, surface, grid, best);
            scan_offset_cell(centre, (dx, ring, dz), atom, surface, grid, best);
        }
    }

    for dy in -inner..=inner {
        for dz in -inner..=inner {
            scan_offset_cell(centre, (-ring, dy, dz), atom, surface, grid, best);
            scan_offset_cell(centre, (ring, dy, dz), atom, surface, grid, best);
        }
    }
}

/// Resolves and scans one shell cell offset from `centre`.
///
/// Overflowing integer coordinates are ignored because they cannot represent a
/// valid neighbouring grid cell.
#[inline]
fn scan_offset_cell(
    centre: CellIndex,
    delta: CellIndex,
    atom: [f32; 3],
    surface: &[[f32; 3]],
    grid: &HashMap<CellIndex, Vec<usize>>,
    best: &mut f64,
) {
    let Some(cell) = offset_cell(centre, delta) else {
        return;
    };

    scan_cell(cell, atom, surface, grid, best);
}

/// Tests all surface points stored in one occupied grid cell.
///
/// Runtime is linear in the bucket's local point count.
fn scan_cell(
    cell: CellIndex,
    atom: [f32; 3],
    surface: &[[f32; 3]],
    grid: &HashMap<CellIndex, Vec<usize>>,
    best: &mut f64,
) {
    let Some(indices) = grid.get(&cell) else {
        return;
    };

    for &index in indices {
        let Some(&point) = surface.get(index) else {
            continue;
        };

        let squared = squared_distance(atom, point);

        if squared < *best {
            *best = squared;
        }
    }
}

/// Adds a signed grid-cell offset without integer overflow.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn offset_cell(centre: CellIndex, delta: CellIndex) -> Option<CellIndex> {
    Some((
        centre.0.checked_add(delta.0)?,
        centre.1.checked_add(delta.1)?,
        centre.2.checked_add(delta.2)?,
    ))
}

/// Returns the largest shell radius required to cover all occupied grid cells.
///
/// Unlike the original origin-based global span, this bound is relative to the
/// queried atom and therefore remains correct even when the atom lies far
/// outside the surface-point coordinate extent.
fn furthest_required_ring(centre: CellIndex, bounds: CellBounds) -> i64 {
    let maximum = [
        cell_distance(centre.0, bounds.min.0),
        cell_distance(centre.0, bounds.max.0),
        cell_distance(centre.1, bounds.min.1),
        cell_distance(centre.1, bounds.max.1),
        cell_distance(centre.2, bounds.min.2),
        cell_distance(centre.2, bounds.max.2),
    ]
    .into_iter()
    .max();
    let Some(distance) = maximum else {
        return 0;
    };
    distance
}

/// Computes the saturated absolute difference between two cell coordinates.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn cell_distance(left: i64, right: i64) -> i64 {
    let difference = (i128::from(left) - i128::from(right)).abs();

    i128_to_i64(difference.min(i128::from(i64::MAX)))
}

/// Returns a conservative squared lower bound for any point in `ring`.
///
/// The atom can lie anywhere inside its centre cell, so a shell `r` has a
/// guaranteed Cartesian separation of at least `(r - 1) * CELL` on one axis.
#[inline]
fn shell_lower_bound_squared(ring: i64, cell_size: f64) -> f64 {
    if ring <= 1 {
        return 0.0;
    }

    let distance = (i64_to_f64(ring) - 1.0) * cell_size;
    distance * distance
}

/// Returns squared Euclidean distance between two points.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f64 {
    let dx = f64::from(left[0]) - f64::from(right[0]);
    let dy = f64::from(left[1]) - f64::from(right[1]);
    let dz = f64::from(left[2]) - f64::from(right[2]);

    dx * dx + dy * dy + dz * dz
}

/// Returns whether all coordinates are finite.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn finite(point: [f32; 3]) -> bool {
    point[0].is_finite() && point[1].is_finite() && point[2].is_finite()
}

#[cfg(test)]
#[path = "depth_tests.rs"]
mod tests;
