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
use crate::numeric::{f64_to_f32, f64_to_usize, usize_to_f64};

/// One integer grid coordinate.
type Cell = (usize, usize, usize);

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

    let expanded = expanded_radii(radii, probe)?;
    let grid = Grid::new(positions, &expanded, grid_options)?;

    let mut state = grid.paint(positions, &expanded);
    grid.flood_exterior(&mut state);

    Ok(grid.collect_cavities(&mut state))
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
    let mut expanded = Vec::with_capacity(radii.len());

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
const SOLID: u8 = 1;
pub(crate) const EXTERIOR: u8 = 2;
const CAVITY: u8 = 3;

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
    pub(crate) fn paint(&self, positions: &[[f32; 3]], expanded: &[f64]) -> Vec<u8> {
        let mut state = vec![EMPTY; self.cell_count()];

        for (position, &radius) in positions.iter().zip(expanded) {
            if radius <= 0.0 {
                continue;
            }

            self.paint_atom(&mut state, position.map(f64::from), radius);
        }

        state
    }

    /// Floods the exterior inward from every empty boundary cell.
    ///
    /// Boundary seeding touches only the six faces rather than scanning every
    /// grid cell. The subsequent flood is `O(cells)` and allocates only its DFS
    /// stack.
    pub(crate) fn flood_exterior(&self, state: &mut [u8]) {
        let mut stack = Vec::new();
        self.seed_boundary(state, &mut stack);

        while let Some((x, y, z)) = stack.pop() {
            self.expand_marked_cell(state, &mut stack, x, y, z, EXTERIOR);
        }
    }

    /// Labels each connected group of unreached empty cells as one cavity.
    ///
    /// The input state is reused as the visitation bitmap, eliminating the
    /// previous full-grid `state.to_vec()` allocation. Runtime is `O(cells)`.
    fn collect_cavities(&self, state: &mut [u8]) -> Vec<Cavity> {
        let cell_volume = self.step * self.step * self.step;
        let mut cavities = Vec::new();

        for z in 0..self.dims[2] {
            for y in 0..self.dims[1] {
                for x in 0..self.dims[0] {
                    let index = self.index(x, y, z);

                    if state.get(index).copied() != Some(EMPTY) {
                        continue;
                    }

                    let count = self.label_cavity(state, (x, y, z));
                    let representative = self.centre(x, y, z);

                    cavities.push(Cavity {
                        volume: usize_to_f64(count) * cell_volume,
                        representative: [
                            f64_to_f32(representative[0]),
                            f64_to_f32(representative[1]),
                            f64_to_f32(representative[2]),
                        ],
                        cells: count,
                    });
                }
            }
        }

        cavities.sort_by(|left, right| right.volume.total_cmp(&left.volume));

        cavities
    }

    /// Returns the total number of cells in the validated grid.
    ///
    /// Construction has already guaranteed that the multiplication fits the
    /// configured grid limit.
    #[inline]
    fn cell_count(&self) -> usize {
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

        for z in lo[2]..=hi[2] {
            let dz = self.axis_centre(2, z) - centre[2];
            let remaining_z = radius_squared - dz * dz;

            if remaining_z < 0.0 {
                continue;
            }

            for y in lo[1]..=hi[1] {
                let dy = self.axis_centre(1, y) - centre[1];
                let remaining_xy = remaining_z - dy * dy;

                if remaining_xy < 0.0 {
                    continue;
                }

                for x in lo[0]..=hi[0] {
                    let dx = self.axis_centre(0, x) - centre[0];

                    if dx * dx <= remaining_xy {
                        let index = self.index(x, y, z);

                        if let Some(cell) = state.get_mut(index) {
                            *cell = SOLID;
                        }
                    }
                }
            }
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

    /// Seeds all empty boundary cells exactly once where practical.
    ///
    /// Face traversal is `O(nx·ny + nx·nz + ny·nz)` instead of `O(cells)`.
    fn seed_boundary(&self, state: &mut [u8], stack: &mut Vec<Cell>) {
        let [nx, ny, nz] = self.dims;

        for z in 0..nz {
            for y in 0..ny {
                self.mark_and_push(state, stack, (0, y, z), EXTERIOR);
                self.mark_and_push(state, stack, (nx - 1, y, z), EXTERIOR);
            }
        }

        for z in 0..nz {
            for x in 1..nx - 1 {
                self.mark_and_push(state, stack, (x, 0, z), EXTERIOR);
                self.mark_and_push(state, stack, (x, ny - 1, z), EXTERIOR);
            }
        }

        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                self.mark_and_push(state, stack, (x, y, 0), EXTERIOR);
                self.mark_and_push(state, stack, (x, y, nz - 1), EXTERIOR);
            }
        }
    }

    /// Labels one enclosed connected component and returns its cell count.
    ///
    /// Each member is marked before it enters the stack, so no cell is pushed
    /// more than once.
    fn label_cavity(&self, state: &mut [u8], seed: Cell) -> usize {
        let mut stack = Vec::new();

        self.mark_and_push(state, &mut stack, seed, CAVITY);

        let mut count = 0usize;

        while let Some((x, y, z)) = stack.pop() {
            count += 1;

            self.expand_marked_cell(state, &mut stack, x, y, z, CAVITY);
        }

        count
    }

    /// Marks one empty cell and appends it to `stack`.
    ///
    /// Non-empty or invalid cells are ignored.
    #[inline]
    fn mark_and_push(&self, state: &mut [u8], stack: &mut Vec<Cell>, cell: Cell, mark: u8) {
        let index = self.index(cell.0, cell.1, cell.2);

        let Some(value) = state.get_mut(index) else {
            return;
        };

        if *value != EMPTY {
            return;
        }

        *value = mark;
        stack.push(cell);
    }

    /// Visits the six in-grid neighbours of one marked cell.
    ///
    /// Unlike the original `neighbours()` implementation, this performs no
    /// heap allocation for a temporary six-element vector.
    fn expand_marked_cell(
        &self,
        state: &mut [u8],
        stack: &mut Vec<Cell>,
        x: usize,
        y: usize,
        z: usize,
        mark: u8,
    ) {
        if x > 0 {
            self.mark_and_push(state, stack, (x - 1, y, z), mark);
        }

        if x + 1 < self.dims[0] {
            self.mark_and_push(state, stack, (x + 1, y, z), mark);
        }

        if y > 0 {
            self.mark_and_push(state, stack, (x, y - 1, z), mark);
        }

        if y + 1 < self.dims[1] {
            self.mark_and_push(state, stack, (x, y + 1, z), mark);
        }

        if z > 0 {
            self.mark_and_push(state, stack, (x, y, z - 1), mark);
        }

        if z + 1 < self.dims[2] {
            self.mark_and_push(state, stack, (x, y, z + 1), mark);
        }
    }
}

#[cfg(test)]
#[path = "cavity_tests.rs"]
mod tests;
