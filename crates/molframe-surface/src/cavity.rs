//! Enclosed cavities found by flooding a grid.
//!
//! A cavity is empty space the solvent cannot reach: a pocket sealed off from the
//! outside by the atoms around it. The molecule is laid on a grid, each cell
//! marked solid if it lies inside a probe-grown atom, and the exterior is flooded
//! inward from the grid boundary. Empty cells the flood never reaches are
//! enclosed, and each connected group of them is one cavity.
//!
//! The solid cells are painted by scattering from each atom's own box, so that
//! part costs `O(atoms · cells per atom)` rather than a cell-by-cell scan of all
//! atoms; the flood and the labelling are linear in the number of cells.

use crate::SurfaceGridOptions;
use crate::accessible_area::SasaError;
use crate::numeric::{f64_to_usize, usize_to_f64};

#[path = "cavity_flood.rs"]
mod flood;
pub(crate) use flood::FloodWorkspace;

/// An enclosed cavity: its volume and a point inside it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cavity {
    /// The cavity volume in cubic ångström.
    pub volume: f64,
    /// A representative point inside the cavity.
    pub representative: [f32; 3],
    /// How many grid cells the cavity occupies.
    pub cells: usize,
}

/// Finds the enclosed cavities of a set of atoms at a grid resolution.
///
/// `radii` are the bare atomic radii and the probe is added internally; a cell
/// is solid when its centre lies within a grown radius of some atom. Cavities
/// are returned largest first. This convenience wrapper uses the named
/// [`SurfaceGridOptions::standard`] allocation profile; use
/// [`cavities_with_options`] to declare a different resource ceiling.
///
/// Runs in `O(atoms · cells per atom + cells)` time.
///
/// # Errors
///
/// Returns [`SasaError::LengthMismatch`], [`SasaError::InvalidProbe`],
/// [`SasaError::InvalidRadius`] for bad inputs, [`SasaError::NoPoints`] for a
/// non-positive resolution, and [`SasaError::GridTooLarge`] when the grid would
/// exceed the cell limit.
pub fn cavities(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    resolution: f32,
) -> Result<Vec<Cavity>, SasaError> {
    cavities_with_options(
        positions,
        radii,
        probe,
        SurfaceGridOptions::standard(resolution),
    )
}

/// Finds enclosed cavities with an explicit grid resolution and allocation ceiling.
///
/// # Errors
///
/// Returns [`SasaError`] for invalid geometry, radii, probe, grid resolution,
/// or allocation ceiling.
pub fn cavities_with_options(
    positions: &[[f32; 3]],
    radii: &[f32],
    probe: f32,
    grid_options: SurfaceGridOptions,
) -> Result<Vec<Cavity>, SasaError> {
    let grid_options = grid_options.validate()?;
    validate_inputs(positions, radii, probe)?;

    if positions.is_empty() {
        return Ok(Vec::new());
    }

    crate::workspace::ensure_surface_input(positions.len(), grid_options.max_workspace_bytes)?;
    let expanded = expanded_radii(radii, probe)?;
    let grid = Grid::new(positions, &expanded, grid_options)?;
    crate::workspace::ensure_cavity_grid(
        positions.len(),
        grid.cell_count(),
        grid_options.max_workspace_bytes,
    )?;

    let mut state = grid.paint(positions, &expanded)?;
    drop(expanded);
    let mut flood = FloodWorkspace::new(grid.cell_count())?;
    grid.flood_exterior(&mut state, &mut flood);

    grid.collect_cavities(&mut state, &mut flood)
}

/// Validates cavity-analysis input dimensions and scalar parameters.
///
/// Runtime is `O(R)` for `R` radii and requires no allocation.
fn validate_inputs(positions: &[[f32; 3]], radii: &[f32], probe: f32) -> Result<(), SasaError> {
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

/// Converts bare radii to double-precision probe-grown radii.
///
/// Runtime and output space are `O(R)`.
fn expanded_radii(radii: &[f32], probe: f32) -> Result<Vec<f64>, SasaError> {
    let probe = f64::from(probe);
    let mut expanded = crate::workspace::empty_with_capacity(radii.len())?;

    for &radius in radii {
        if !radius.is_finite() || radius < 0.0 {
            return Err(SasaError::InvalidRadius);
        }

        expanded.push(f64::from(radius) + probe);
    }

    Ok(expanded)
}

/// The regular grid the molecule sits in.
pub(crate) struct Grid {
    pub(crate) origin: [f64; 3],
    pub(crate) step: f64,
    pub(crate) dims: [usize; 3],
}

/// Cell markings during the flood.
pub(crate) const EMPTY: u8 = 0;
pub(super) const SOLID: u8 = 1;
pub(crate) const EXTERIOR: u8 = 2;
pub(super) const CAVITY: u8 = 3;

impl Grid {
    /// Builds a grid enclosing the grown atoms with a one-cell margin.
    ///
    /// Runtime is `O(A)` for `A` atoms and additional space is `O(1)`.
    pub(crate) fn new(
        positions: &[[f32; 3]],
        expanded: &[f64],
        options: SurfaceGridOptions,
    ) -> Result<Grid, SasaError> {
        let step = f64::from(options.resolution);
        let (origin, dims) =
            crate::cavity_geometry::geometry(positions, expanded, step, options.max_cells)?;

        Ok(Grid { origin, step, dims })
    }

    /// Converts three-dimensional grid coordinates into a linear index.
    ///
    /// Runtime and additional space are `O(1)`.
    #[inline]
    pub(crate) fn index(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.dims[1] + y) * self.dims[0] + x
    }

    /// Returns the Cartesian centre of one grid cell.
    ///
    /// Runtime and additional space are `O(1)`.
    #[inline]
    pub(crate) fn centre(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.axis_centre(0, x),
            self.axis_centre(1, y),
            self.axis_centre(2, z),
        ]
    }

    /// Marks every cell inside a grown atom as solid by scattering from atoms.
    ///
    /// Runtime is proportional to the cells in the atoms' local bounding boxes.
    /// The only allocation is the final `O(cells)` state buffer.
    pub(crate) fn paint(
        &self,
        positions: &[[f32; 3]],
        expanded: &[f64],
    ) -> Result<Vec<u8>, SasaError> {
        let mut state = crate::workspace::filled(self.cell_count(), EMPTY)?;

        for (position, &radius) in positions.iter().zip(expanded) {
            if radius <= 0.0 {
                continue;
            }

            self.paint_atom(&mut state, position.map(f64::from), radius);
        }

        Ok(state)
    }

    /// Returns the total number of cells in the validated grid.
    ///
    /// Construction has already guaranteed that the multiplication fits the
    /// configured grid limit.
    #[inline]
    pub(crate) fn cell_count(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Returns one cell-centre coordinate along `axis`.
    ///
    /// Runtime and additional space are `O(1)`.
    #[inline]
    fn axis_centre(&self, axis: usize, coordinate: usize) -> f64 {
        self.origin[axis] + (usize_to_f64(coordinate) + 0.5) * self.step
    }

    /// Marks the cells intersecting one atom's grown sphere.
    ///
    /// Squared Z and Y contributions are rejected before entering the X loop,
    /// reducing arithmetic for cells outside the sphere.
    fn paint_atom(&self, state: &mut [u8], centre: [f64; 3], radius: f64) {
        let Some((lo, hi)) = self.cell_box(centre, radius) else {
            return;
        };

        let radius_squared = radius * radius;
        let first_x = self.axis_centre(0, lo[0]) - centre[0];
        let first_y = self.axis_centre(1, lo[1]) - centre[1];
        let mut dz = self.axis_centre(2, lo[2]) - centre[2];

        for z in lo[2]..=hi[2] {
            let remaining_z = radius_squared - dz * dz;

            if remaining_z < 0.0 {
                dz += self.step;
                continue;
            }

            let mut dy = first_y;
            for y in lo[1]..=hi[1] {
                let remaining_xy = remaining_z - dy * dy;

                if remaining_xy < 0.0 {
                    dy += self.step;
                    continue;
                }

                let mut dx = first_x;
                for (index, _) in (self.index(lo[0], y, z)..).zip(lo[0]..=hi[0]) {
                    if dx * dx <= remaining_xy
                        && let Some(cell) = state.get_mut(index)
                    {
                        *cell = SOLID;
                    }
                    dx += self.step;
                }
                dy += self.step;
            }
            dz += self.step;
        }
    }

    /// The inclusive cell-index box covering sphere-centre candidates.
    ///
    /// Bounds account for the half-cell centre offset, making the candidate box
    /// tighter than the previous whole-cell-coordinate bound.
    fn cell_box(&self, centre: [f64; 3], radius: f64) -> Option<([usize; 3], [usize; 3])> {
        let mut lo = [0usize; 3];
        let mut hi = [0usize; 3];

        for axis in 0..3 {
            let (start, end) = self.cell_axis_range(axis, centre[axis], radius)?;

            lo[axis] = start;
            hi[axis] = end;
        }

        Some((lo, hi))
    }

    /// Returns the inclusive cell-centre index interval touched on one axis.
    ///
    /// Runtime and additional space are `O(1)`.
    fn cell_axis_range(&self, axis: usize, centre: f64, radius: f64) -> Option<(usize, usize)> {
        let last = self.dims[axis].checked_sub(1)?;

        let start = ((centre - radius - self.origin[axis]) / self.step - 0.5).ceil();

        let end = ((centre + radius - self.origin[axis]) / self.step - 0.5).floor();

        if end < 0.0 || start > usize_to_f64(last) {
            return None;
        }

        let start = if start <= 0.0 {
            0
        } else {
            f64_to_usize(start).min(last)
        };

        let end = if end <= 0.0 {
            0
        } else {
            f64_to_usize(end).min(last)
        };

        (start <= end).then_some((start, end))
    }
}

#[cfg(test)]
#[path = "cavity_tests.rs"]
mod tests;
