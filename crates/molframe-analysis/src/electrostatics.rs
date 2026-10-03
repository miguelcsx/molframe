//! Finite-range screened Coulomb contact potentials, not Poisson–Boltzmann.

use crate::numeric::{f64_to_f32, usize_to_f64};
use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockExecutionError, BlockPlan, try_for_each_block_in};
use molframe_core::structure::Structure;
use molframe_spatial::CellList;

/// A Cartesian affine grid with x-fastest storage.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GridSpec {
    /// Row-major affine acting on column vectors `(x, y, z, 1)` in ångström.
    pub voxel_to_world: [[f64; 4]; 4],
    /// Voxel counts along canonical x, y, z.
    pub dimensions: [usize; 3],
}

/// A dense potential field in thermal-voltage units (kT/e at 298 K).
#[derive(Clone, Debug, PartialEq)]
pub struct ScalarGrid {
    /// Geometry used for evaluation.
    pub spec: GridSpec,
    /// Values in `(z * ny + y) * nx + x` order.
    pub values: Vec<f64>,
}

/// Why a contact potential could not be evaluated.
#[derive(Clone, Debug, thiserror::Error, PartialEq)]
pub enum PotentialError {
    /// Coordinates must be dense and charges must align with them.
    #[error("contact potential requires one dense coordinate and charge per atom")]
    MisalignedInput,
    /// Every coordinate and charge must be finite.
    #[error("contact potential input must be finite")]
    NonFiniteInput,
    /// Grid geometry must be finite, affine, nondegenerate and nonempty.
    #[error("invalid contact potential grid")]
    InvalidGrid,
    /// The grid cannot be indexed on this platform.
    #[error("contact potential grid exceeds the platform index range")]
    GridTooLarge,
    /// Spatial indexing failed.
    #[error("contact potential spatial query failed: {0}")]
    Spatial(molframe_spatial::SpatialError),
    /// A worker could not complete its output.
    #[error("a contact potential worker panicked")]
    WorkerPanicked,
}

/// Voxels per block: small enough to balance uneven local densities.
const VOXELS_PER_BLOCK: usize = 4096;

// Exact elementary charge and Boltzmann constant; CODATA vacuum permittivity.
const THERMAL_COULOMB: f64 = 1.602_176_634e-19 * 1.602_176_634e-19
    / (4.0 * std::f64::consts::PI * 8.854_187_812_8e-12 * 1.0e-10 * 1.380_649e-23 * 298.0);
/// Default finite interaction radius in ångström.
pub const CONTACT_POTENTIAL_CUTOFF: f32 = 12.0;

/// Computes a screened contact potential using the default worker budget.
///
/// Charges are in elementary-charge units. Each atom contributes
/// `q * e² / (4π ε₀ kT Å) / (4 * max(r, 1 Å)²)` within `cutoff`.
/// The distance-dependent dielectric is `4r`; softening makes the value finite
/// at atom centres. This is a local contact field, not a solvent boundary-value
/// calculation or a Poisson–Boltzmann electrostatic potential.
///
/// # Errors
///
/// Returns [`PotentialError`] for misaligned/nonfinite inputs, invalid geometry,
/// spatial-index failure, or worker failure.
///
/// # Examples
///
/// ```
/// use molframe_analysis::{GridSpec, ScalarGrid, PotentialError, contact_potential};
/// use molframe_core::structure::Structure;
/// fn field(structure: &Structure, charges: &[f64]) -> Result<ScalarGrid, PotentialError> {
///     contact_potential(structure, charges, GridSpec {
///         dimensions: [2, 1, 1],
///         voxel_to_world: [[5.,0.,0.,5.],[0.,1.,0.,0.],[0.,0.,1.,0.],[0.,0.,0.,1.]],
///     }, 12.0)
/// }
/// ```
pub fn contact_potential(
    structure: &Structure,
    charges: &[f64],
    spec: GridSpec,
    cutoff: f32,
) -> Result<ScalarGrid, PotentialError> {
    contact_potential_in(
        structure,
        charges,
        spec,
        cutoff,
        &ExecutionContext::default(),
    )
}

/// [`contact_potential`] with an explicit shared worker budget.
///
/// Fixed voxel blocks and deterministic neighbor order make output independent
/// of worker count. Coordinates are borrowed; no expanded atom/voxel cloud is
/// constructed. Index memory is O(atoms + cells); output memory is O(voxels).
///
/// # Errors
///
/// Returns the same errors as [`contact_potential`].
pub fn contact_potential_in(
    structure: &Structure,
    charges: &[f64],
    spec: GridSpec,
    cutoff: f32,
    context: &ExecutionContext,
) -> Result<ScalarGrid, PotentialError> {
    let positions = structure.positions();
    if !structure.data().coords.is_dense()
        || positions.len() != structure.atom_count() as usize
        || positions.len() != charges.len()
    {
        return Err(PotentialError::MisalignedInput);
    }
    if positions.iter().flatten().any(|value| !value.is_finite())
        || charges.iter().any(|value| !value.is_finite())
    {
        return Err(PotentialError::NonFiniteInput);
    }
    let count = spec.validate()?;
    if !cutoff.is_finite() || cutoff < 0.0 {
        return Err(PotentialError::Spatial(
            molframe_spatial::SpatialError::InvalidCutoff,
        ));
    }
    let search_radius = conservative_radius(spec, positions, cutoff)?;
    let targets: Vec<_> = (0..structure.atom_count()).collect();
    let index = CellList::build(positions, &targets, search_radius, None)
        .map_err(PotentialError::Spatial)?;
    let cutoff_squared = f64::from(cutoff).powi(2);
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| PotentialError::GridTooLarge)?;
    // Voxel blocks are appended as they complete, in block order, so the peak
    // is the grid plus a bounded window of blocks rather than the grid twice.
    try_for_each_block_in(
        BlockPlan::new(count, VOXELS_PER_BLOCK),
        context,
        VOXELS_PER_BLOCK * size_of::<f64>(),
        |_, range| {
            range
                .map(|voxel| {
                    let point = spec.point(voxel);
                    let query = point.map(f64_to_f32);
                    let mut potential = 0.0;
                    index
                        .for_each_neighbor(query, search_radius, |atom, _| {
                            let atom = atom as usize;
                            let squared: f64 = point
                                .iter()
                                .zip(positions[atom])
                                .map(|(left, right)| (left - f64::from(right)).powi(2))
                                .sum();
                            if squared <= cutoff_squared {
                                potential +=
                                    charges[atom] * THERMAL_COULOMB / (4.0 * squared.max(1.0));
                            }
                        })
                        .map_err(PotentialError::Spatial)?;
                    Ok(potential)
                })
                .collect::<Result<Vec<f64>, PotentialError>>()
        },
        |block| {
            values.extend_from_slice(&block);
            Ok(())
        },
    )
    .map_err(|error| match error {
        BlockExecutionError::Memory(error) => {
            PotentialError::Spatial(molframe_spatial::SpatialError::Memory(error))
        }
        BlockExecutionError::Cancelled => {
            PotentialError::Spatial(molframe_spatial::SpatialError::Cancelled)
        }
        BlockExecutionError::Worker(_) => PotentialError::WorkerPanicked,
        BlockExecutionError::Operation(error) => error,
    })?;
    Ok(ScalarGrid { spec, values })
}

impl GridSpec {
    fn validate(self) -> Result<usize, PotentialError> {
        let m = self.voxel_to_world;
        if self.dimensions.contains(&0)
            || m.iter().flatten().any(|value| !value.is_finite())
            || m[3]
                .iter()
                .zip([0.0, 0.0, 0.0, 1.0])
                .any(|(value, expected)| (value - expected).abs() > 0.0)
        {
            return Err(PotentialError::InvalidGrid);
        }
        let determinant = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if !determinant.is_finite() || determinant.abs() <= 0.0 {
            return Err(PotentialError::InvalidGrid);
        }
        self.dimensions
            .into_iter()
            .try_fold(1usize, |count, dimension| {
                count
                    .checked_mul(dimension)
                    .filter(|count| *count <= isize::MAX as usize / size_of::<f64>())
                    .ok_or(PotentialError::GridTooLarge)
            })
    }

    fn point(self, voxel: usize) -> [f64; 3] {
        let [nx, ny, _] = self.dimensions;
        let indices = [voxel % nx, voxel / nx % ny, voxel / nx / ny].map(usize_to_f64);
        std::array::from_fn(|row| {
            self.voxel_to_world[row][3]
                + (0..3)
                    .map(|axis| self.voxel_to_world[row][axis] * indices[axis])
                    .sum::<f64>()
        })
    }
}

fn conservative_radius(
    spec: GridSpec,
    positions: &[[f32; 3]],
    cutoff: f32,
) -> Result<f32, PotentialError> {
    let extent = spec.dimensions.map(|dimension| usize_to_f64(dimension - 1));
    let mut scale = positions
        .iter()
        .flatten()
        .map(|value| f64::from(value.abs()))
        .fold(0.0, f64::max);
    for row in spec.voxel_to_world.iter().take(3) {
        // Bounds every affine grid point, including skew and reflected axes.
        let bound = row[3].abs()
            + (0..3)
                .map(|axis| row[axis].abs() * extent[axis])
                .sum::<f64>();
        if !bound.is_finite() || bound > f64::from(f32::MAX) {
            return Err(PotentialError::InvalidGrid);
        }
        scale = scale.max(bound);
    }
    // Cell-list distances are f32; inflate only candidate search, then test
    // the scientific cutoff in f64. This prevents boundary false negatives.
    let radius =
        f64_to_f32(f64::from(cutoff) + (scale + f64::from(cutoff)) * f64::from(f32::EPSILON) * 8.0);
    if !radius.is_finite() {
        return Err(PotentialError::InvalidGrid);
    }
    Ok(radius)
}

#[cfg(test)]
#[path = "electrostatics_tests.rs"]
mod tests;
