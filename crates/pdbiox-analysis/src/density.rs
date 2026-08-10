//! One- and three-dimensional densities with explicit grids and weights.

use crate::numeric::{f32_to_usize, usize_to_f32};

/// Cartesian axis for a linear density profile.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CartesianAxis {
    /// X axis.
    X,
    /// Y axis.
    Y,
    /// Z axis.
    Z,
}

impl CartesianAxis {
    const fn index(self) -> usize {
        match self {
            Self::X => 0,
            Self::Y => 1,
            Self::Z => 2,
        }
    }
}

/// Explicit one-dimensional histogram definition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearDensityOptions {
    /// Profile axis.
    pub axis: CartesianAxis,
    /// Inclusive coordinate lower bound.
    pub minimum: f32,
    /// Exclusive coordinate upper bound.
    pub maximum: f32,
    /// Number of equal-width bins.
    pub bins: usize,
}

/// One linear-density interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearDensityBin {
    /// Inclusive lower coordinate.
    pub lower: f32,
    /// Exclusive upper coordinate.
    pub upper: f32,
    /// Sum of caller-provided weights in the interval.
    pub weight: f64,
    /// Weight per unit length.
    pub density: f64,
}

/// Explicit rectilinear density-grid definition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DensityGridSpec {
    /// Cartesian lower corner of the grid.
    pub origin: [f32; 3],
    /// Positive voxel edge lengths.
    pub spacing: [f32; 3],
    /// Voxel counts along x, y, and z.
    pub shape: [usize; 3],
}

/// A dense x-major, then y, then z density grid.
#[derive(Clone, Debug, PartialEq)]
pub struct DensityGrid {
    /// Grid definition used for accumulation.
    pub spec: DensityGridSpec,
    /// Weight per voxel volume in flat `(x * ny + y) * nz + z` order.
    pub density: Vec<f64>,
    /// Total weight outside the half-open grid extent.
    pub excluded_weight: f64,
}

/// Why density accumulation could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq)]
pub enum DensityError {
    /// Positions and weights are atom-aligned inputs.
    #[error("expected {positions} weights, found {weights}")]
    LengthMismatch {
        /// Number of positions.
        positions: usize,
        /// Number of weights.
        weights: usize,
    },
    /// Histogram bounds or bin count are invalid.
    #[error("density bounds must be finite and increasing and bins must be non-zero")]
    InvalidBins,
    /// Grid origin, spacing, or shape is invalid.
    #[error("grid origin must be finite; spacing and shape must be positive")]
    InvalidGrid,
    /// Grid shape overflows addressable memory.
    #[error("grid shape overflows the platform index range")]
    GridTooLarge,
    /// Coordinates and weights must be finite.
    #[error("coordinates and weights must be finite")]
    NonFiniteInput,
}

/// Computes a weighted linear density profile.
///
/// Callers choose the meaning and units of `weights` (occupancy, atom count,
/// mass, or charge); no chemistry or unit convention is inferred.
///
/// # Errors
///
/// Returns [`DensityError`] for misaligned, non-finite, or invalid inputs.
pub fn linear_density(
    positions: &[[f32; 3]],
    weights: &[f64],
    options: LinearDensityOptions,
) -> Result<Vec<LinearDensityBin>, DensityError> {
    validate_aligned(positions, weights)?;
    if !options.minimum.is_finite()
        || !options.maximum.is_finite()
        || options.maximum <= options.minimum
        || options.bins == 0
    {
        return Err(DensityError::InvalidBins);
    }
    let width = (options.maximum - options.minimum) / usize_to_f32(options.bins);
    let mut sums = vec![0.0; options.bins];
    for (position, weight) in positions.iter().zip(weights) {
        let coordinate = position[options.axis.index()];
        if coordinate >= options.minimum && coordinate < options.maximum {
            let Some(index) = f32_to_usize((coordinate - options.minimum) / width) else {
                continue;
            };
            if let Some(sum) = sums.get_mut(index) {
                *sum += weight;
            }
        }
    }
    Ok(sums
        .into_iter()
        .enumerate()
        .map(|(index, weight)| {
            let lower = options.minimum + usize_to_f32(index) * width;
            LinearDensityBin {
                lower,
                upper: lower + width,
                weight,
                density: weight / f64::from(width),
            }
        })
        .collect())
}

/// Accumulates caller-provided weights into a rectilinear 3D density map.
///
/// # Errors
///
/// Returns [`DensityError`] for invalid, oversized, misaligned, or non-finite input.
pub fn density_map(
    positions: &[[f32; 3]],
    weights: &[f64],
    spec: DensityGridSpec,
) -> Result<DensityGrid, DensityError> {
    validate_aligned(positions, weights)?;
    validate_grid(spec)?;
    let voxel_count = spec
        .shape
        .into_iter()
        .try_fold(1_usize, usize::checked_mul)
        .ok_or(DensityError::GridTooLarge)?;
    let mut density = vec![0.0; voxel_count];
    let mut excluded_weight = 0.0;
    for (position, weight) in positions.iter().zip(weights) {
        if let Some([x, y, z]) = voxel_index(*position, spec) {
            density[(x * spec.shape[1] + y) * spec.shape[2] + z] += weight;
        } else {
            excluded_weight += weight;
        }
    }
    let volume = spec.spacing.into_iter().map(f64::from).product::<f64>();
    for value in &mut density {
        *value /= volume;
    }
    Ok(DensityGrid {
        spec,
        density,
        excluded_weight,
    })
}

fn validate_aligned(positions: &[[f32; 3]], weights: &[f64]) -> Result<(), DensityError> {
    if positions.len() != weights.len() {
        return Err(DensityError::LengthMismatch {
            positions: positions.len(),
            weights: weights.len(),
        });
    }
    if positions.iter().flatten().all(|value| value.is_finite())
        && weights.iter().all(|value| value.is_finite())
    {
        Ok(())
    } else {
        Err(DensityError::NonFiniteInput)
    }
}

fn validate_grid(spec: DensityGridSpec) -> Result<(), DensityError> {
    if spec.origin.iter().all(|value| value.is_finite())
        && spec
            .spacing
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
        && spec.shape.iter().all(|count| *count > 0)
    {
        Ok(())
    } else {
        Err(DensityError::InvalidGrid)
    }
}

fn voxel_index(position: [f32; 3], spec: DensityGridSpec) -> Option<[usize; 3]> {
    let mut index = [0_usize; 3];
    for axis in 0..3 {
        let relative = (position[axis] - spec.origin[axis]) / spec.spacing[axis];
        if relative < 0.0 || relative >= usize_to_f32(spec.shape[axis]) {
            return None;
        }
        index[axis] = f32_to_usize(relative)?;
    }
    Some(index)
}

#[cfg(test)]
#[path = "density_tests.rs"]
mod tests;
