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
use crate::cavity::{FloodWorkspace, Grid};

#[path = "ses_distance.rs"]
mod distance;
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
/// exact separable Euclidean distance transform and polygonisation passes are
/// linear in the grid size.
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

    crate::workspace::ensure_surface_input(positions.len(), grid_options.max_workspace_bytes)?;
    let expanded = expanded_radii(radii, probe)?;
    let grid = Grid::new(positions, &expanded, grid_options)?;
    let longest_line = grid
        .dims
        .into_iter()
        .max()
        .ok_or(SasaError::GridDimensionsOverflow)?;
    crate::workspace::ensure_ses_grid(
        positions.len(),
        grid.cell_count(),
        longest_line,
        grid_options.max_workspace_bytes,
    )?;

    let mut probe_centres = grid.paint(positions, &expanded)?;
    drop(expanded);
    let mut flood = FloodWorkspace::new(grid.cell_count())?;
    grid.flood_exterior(&mut probe_centres, &mut flood);
    drop(flood);

    let distance = distance::distance_from_exterior(&grid, probe_centres)?;
    let level = probe;
    let triangle_capacity = mesh::triangle_capacity(&grid, &distance, level)?;
    crate::workspace::ensure_ses_output(
        positions.len(),
        grid.cell_count(),
        triangle_capacity,
        grid_options.max_workspace_bytes,
    )?;

    mesh::polygonise(&grid, &distance, level, triangle_capacity)
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
fn expanded_radii(radii: &[f32], probe: f32) -> Result<Vec<f64>, SasaError> {
    let probe = f64::from(probe);
    let mut expanded = crate::workspace::empty_with_capacity(radii.len())?;
    for &radius in radii {
        expanded.push(f64::from(radius) + probe);
    }
    Ok(expanded)
}

#[cfg(test)]
#[path = "ses_tests.rs"]
mod tests;
