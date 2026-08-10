//! Solvent-excluded surface from a probe-centre accessibility grid.
//!
//! Probe centres outside the union of probe-grown atoms are flooded from the
//! boundary. The molecular region is the complement of all points a probe
//! centred in that exterior can reach. This is the rolling-probe definition of
//! the solvent-excluded surface, including contact and re-entrant regions.
//! A six-tetrahedra decomposition polygonises the boundary without the axis
//! bias of exposed voxel faces.

use crate::SurfaceGridOptions;
use crate::accessible_area::SasaError;
use crate::cavity::{EXTERIOR, Grid};
use crate::numeric::isize_to_f64;

#[path = "ses_mesh.rs"]
mod mesh;

/// One triangle of a solvent-excluded surface mesh.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceTriangle {
    /// Three triangle vertices.
    pub vertices: [[f32; 3]; 3],
    /// Unit face normal.
    pub normal: [f32; 3],
    /// Triangle area.
    pub area: f64,
}

/// Deterministic triangle mesh and its total area.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct SolventExcludedSurface {
    /// Boundary triangles in grid traversal order.
    pub triangles: Vec<SurfaceTriangle>,
    /// Sum of triangle areas.
    pub area: f64,
}

/// Computes a rolling-probe solvent-excluded surface (Connolly surface).
///
/// Smaller `resolution` values converge on the continuous surface and cost more
/// memory. Atom painting scatters only over each atom's local grid box; the
/// remaining passes are linear in grid cells times the fixed probe stencil.
/// This convenience wrapper uses [`SurfaceGridOptions::standard`]; use
/// [`solvent_excluded_surface_with_options`] for a caller-declared cell budget.
///
/// # Errors
///
/// Returns [`SasaError`] for mismatched or invalid radii/probe/resolution, or
/// when the bounded grid would be too large.
pub fn solvent_excluded_surface(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    resolution: f32,
) -> Result<SolventExcludedSurface, SasaError> {
    solvent_excluded_surface_with_options(
        positions,
        radii,
        probe,
        SurfaceGridOptions::standard(resolution),
    )
}

/// Computes a rolling-probe surface with explicit grid resource controls.
///
/// # Errors
///
/// Returns [`SasaError`] for invalid geometry, radii, probe, grid resolution,
/// or allocation ceiling.
pub fn solvent_excluded_surface_with_options(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    grid_options: SurfaceGridOptions,
) -> Result<SolventExcludedSurface, SasaError> {
    let grid_options = grid_options.validate()?;
    validate(positions, radii, probe)?;

    if positions.is_empty() {
        return Ok(SolventExcludedSurface::default());
    }

    let expanded = expanded_radii(radii, probe);

    let grid = Grid::new(positions, &expanded, grid_options)?;

    let mut probe_centres = grid.paint(positions, &expanded);
    grid.flood_exterior(&mut probe_centres);

    let distance = distance_from_exterior(&grid, &probe_centres);

    Ok(mesh::polygonise(&grid, &distance, f64::from(probe)))
}

/// Validates solvent-excluded-surface input.
///
/// Runtime is `O(R)` for radius validation and requires no allocation.
fn validate(positions: &[[f32; 3]], radii: &[f32], probe: f32) -> Result<(), SasaError> {
    if positions.len() != radii.len() {
        return Err(SasaError::LengthMismatch {
            positions: positions.len(),
            radii: radii.len(),
        });
    }

    if !probe.is_finite() || probe < 0.0 {
        return Err(SasaError::InvalidProbe);
    }

    if radii
        .iter()
        .any(|radius| !radius.is_finite() || *radius < 0.0)
    {
        return Err(SasaError::InvalidRadius);
    }

    Ok(())
}

/// Builds double-precision probe-grown radii.
///
/// Runtime and output space are `O(R)`.
fn expanded_radii(radii: &[f32], probe: f32) -> Vec<f64> {
    let probe = f64::from(probe);

    radii
        .iter()
        .map(|radius| f64::from(*radius) + probe)
        .collect()
}

/// One chamfer-neighbour displacement and its Cartesian cost.
#[derive(Clone, Copy)]
struct ChamferOffset {
    delta: [isize; 3],
    cost: f64,
}

/// The thirteen previously visited neighbours in a 3×3×3 chamfer stencil.
const FORWARD_DELTAS: [[isize; 3]; 13] = [
    [-1, -1, -1],
    [0, -1, -1],
    [1, -1, -1],
    [-1, 0, -1],
    [0, 0, -1],
    [1, 0, -1],
    [-1, 1, -1],
    [0, 1, -1],
    [1, 1, -1],
    [-1, -1, 0],
    [0, -1, 0],
    [1, -1, 0],
    [-1, 0, 0],
];

/// Computes the chamfer distance transform from all exterior cells.
///
/// The two passes each visit every grid cell and a fixed thirteen-neighbour
/// stencil, giving `O(cells)` time and `O(cells)` distance storage.
fn distance_from_exterior(grid: &Grid, state: &[u8]) -> Vec<f64> {
    let mut distance = initial_distances(state);

    let forward = chamfer_offsets(false, grid.step);
    relax_forward(grid, &mut distance, &forward);

    let backward = chamfer_offsets(true, grid.step);
    relax_backward(grid, &mut distance, &backward);

    distance
}

/// Creates the initial distance field with exterior cells as zero.
///
/// Runtime and output space are `O(cells)`.
fn initial_distances(state: &[u8]) -> Vec<f64> {
    state
        .iter()
        .map(|cell| {
            if *cell == EXTERIOR {
                0.0
            } else {
                f64::INFINITY
            }
        })
        .collect()
}

/// Builds the fixed thirteen-element chamfer stencil.
///
/// A stack-allocated array replaces the previous heap-allocated `Vec`.
fn chamfer_offsets(reverse: bool, step: f64) -> [ChamferOffset; 13] {
    std::array::from_fn(|index| {
        let delta = FORWARD_DELTAS[index];
        let delta = if reverse {
            [-delta[0], -delta[1], -delta[2]]
        } else {
            delta
        };

        let squared = isize_to_f64(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);

        ChamferOffset {
            delta,
            cost: squared.sqrt() * step,
        }
    })
}

/// Executes the forward chamfer relaxation pass.
///
/// Runtime is `13 · O(cells)` with `O(1)` auxiliary storage.
fn relax_forward(grid: &Grid, distance: &mut [f64], offsets: &[ChamferOffset; 13]) {
    for z in 0..grid.dims[2] {
        for y in 0..grid.dims[1] {
            for x in 0..grid.dims[0] {
                relax(grid, distance, x, y, z, offsets);
            }
        }
    }
}

/// Executes the backward chamfer relaxation pass.
///
/// Runtime is `13 · O(cells)` with `O(1)` auxiliary storage.
fn relax_backward(grid: &Grid, distance: &mut [f64], offsets: &[ChamferOffset; 13]) {
    for z in (0..grid.dims[2]).rev() {
        for y in (0..grid.dims[1]).rev() {
            for x in (0..grid.dims[0]).rev() {
                relax(grid, distance, x, y, z, offsets);
            }
        }
    }
}

/// Relaxes one distance-transform grid cell against a fixed stencil.
///
/// Runtime is `O(13)` and auxiliary space is `O(1)`.
fn relax(
    grid: &Grid,
    distance: &mut [f64],
    x: usize,
    y: usize,
    z: usize,
    offsets: &[ChamferOffset; 13],
) {
    let index = grid.index(x, y, z);
    let Some(&initial) = distance.get(index) else {
        return;
    };

    let mut best = initial;

    for offset in offsets {
        let Some(neighbour) = neighbour_index(grid, x, y, z, offset.delta) else {
            continue;
        };

        let Some(&neighbour_distance) = distance.get(neighbour) else {
            continue;
        };

        best = best.min(neighbour_distance + offset.cost);
    }

    if let Some(cell) = distance.get_mut(index) {
        *cell = best;
    }
}

/// Resolves one signed neighbour displacement to a linear grid index.
///
/// Runtime and auxiliary space are `O(1)`.
fn neighbour_index(grid: &Grid, x: usize, y: usize, z: usize, delta: [isize; 3]) -> Option<usize> {
    let nx = x.checked_add_signed(delta[0])?;
    let ny = y.checked_add_signed(delta[1])?;
    let nz = z.checked_add_signed(delta[2])?;

    if nx >= grid.dims[0] || ny >= grid.dims[1] || nz >= grid.dims[2] {
        return None;
    }

    Some(grid.index(nx, ny, nz))
}

#[cfg(test)]
#[path = "ses_tests.rs"]
mod tests;
