//! Deterministic streamline integration over caller-supplied Cartesian vector grids.
//!
//! Trilinear sampling is constant time and one streamline costs
//! `O(max_steps)`. The integrator follows the normalized field, so its step
//! size is a spatial distance rather than a field-dependent time increment.

use num_traits::ToPrimitive;

#[cfg(test)]
#[path = "vector_field_tests.rs"]
mod tests;

/// Dense x-major, then y, then z Cartesian vector field.
#[derive(Clone, Debug, PartialEq)]
pub struct VectorFieldGrid {
    /// Cartesian lower corner.
    pub origin: [f32; 3],
    /// Positive sample spacing.
    pub spacing: [f32; 3],
    /// Sample counts along x, y, and z.
    pub shape: [usize; 3],
    /// Vectors in flat `(x * ny + y) * nz + z` order.
    pub vectors: Vec<[f32; 3]>,
}

impl VectorFieldGrid {
    /// Validates one dense vector field.
    ///
    /// # Errors
    ///
    /// Returns [`VectorFieldError`] for malformed geometry, length or values.
    pub fn new(
        origin: [f32; 3],
        spacing: [f32; 3],
        shape: [usize; 3],
        vectors: Vec<[f32; 3]>,
    ) -> Result<Self, VectorFieldError> {
        if !origin.iter().all(|value| value.is_finite())
            || !spacing
                .iter()
                .all(|value| value.is_finite() && *value > 0.0)
            || !shape.iter().all(|count| *count >= 2)
        {
            return Err(VectorFieldError::InvalidGrid);
        }
        let expected = shape
            .into_iter()
            .try_fold(1_usize, usize::checked_mul)
            .ok_or(VectorFieldError::GridTooLarge)?;
        if vectors.len() != expected {
            return Err(VectorFieldError::LengthMismatch {
                expected,
                actual: vectors.len(),
            });
        }
        if !vectors.iter().flatten().all(|value| value.is_finite()) {
            return Err(VectorFieldError::NonFiniteInput);
        }
        Ok(Self {
            origin,
            spacing,
            shape,
            vectors,
        })
    }

    /// Trilinearly samples a point inside the closed sample domain.
    #[must_use]
    pub fn sample(&self, point: [f32; 3]) -> Option<[f32; 3]> {
        if !point.iter().all(|value| value.is_finite()) {
            return None;
        }
        let mut lower = [0_usize; 3];
        let mut fraction = [0.0_f32; 3];
        for axis in 0..3 {
            let coordinate = (point[axis] - self.origin[axis]) / self.spacing[axis];
            let maximum = self.shape[axis].checked_sub(1)?.to_f32()?;
            if coordinate < 0.0 || coordinate > maximum {
                return None;
            }
            let floor = coordinate.floor().to_usize()?;
            lower[axis] = floor.min(self.shape[axis] - 2);
            fraction[axis] = coordinate - lower[axis].to_f32()?;
        }
        let mut sampled = [0.0_f32; 3];
        for dx in 0..=1 {
            for dy in 0..=1 {
                for dz in 0..=1 {
                    let weight = axis_weight(fraction[0], dx)
                        * axis_weight(fraction[1], dy)
                        * axis_weight(fraction[2], dz);
                    let value = self.value(lower[0] + dx, lower[1] + dy, lower[2] + dz)?;
                    for axis in 0..3 {
                        sampled[axis] += value[axis] * weight;
                    }
                }
            }
        }
        Some(sampled)
    }

    fn value(&self, x: usize, y: usize, z: usize) -> Option<[f32; 3]> {
        self.vectors
            .get((x * self.shape[1] + y) * self.shape[2] + z)
            .copied()
    }
}

/// Direction in which a seed is integrated.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum StreamlineDirection {
    /// Follow the field.
    #[default]
    Forward,
    /// Follow the negated field.
    Backward,
    /// Integrate both directions and join them at the seed.
    Both,
}

/// Bounded spatial streamline integration controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StreamlineOptions {
    /// Spatial RK4 step length.
    pub step_size: f32,
    /// Hard step ceiling per direction.
    pub max_steps: usize,
    /// Hard accumulated arc-length ceiling per direction.
    pub max_length: f32,
    /// Field magnitude below which integration terminates.
    pub min_speed: f32,
    /// Integration direction.
    pub direction: StreamlineDirection,
}

impl Default for StreamlineOptions {
    fn default() -> Self {
        Self {
            step_size: 0.25,
            max_steps: 4096,
            max_length: 1000.0,
            min_speed: 1.0e-6,
            direction: StreamlineDirection::Forward,
        }
    }
}

/// Why a vector-field operation could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq)]
pub enum VectorFieldError {
    /// Grid geometry is malformed.
    #[error(
        "vector grids require finite origins, positive spacing, and at least two samples per axis"
    )]
    InvalidGrid,
    /// Grid dimensions overflow addressable memory.
    #[error("vector grid shape overflows the platform index range")]
    GridTooLarge,
    /// Vector storage does not match the declared grid.
    #[error("expected {expected} vectors, found {actual}")]
    LengthMismatch {
        /// Declared number of samples.
        expected: usize,
        /// Supplied number of samples.
        actual: usize,
    },
    /// A vector or seed is not finite.
    #[error("vector fields and seeds must be finite")]
    NonFiniteInput,
    /// Integration controls are malformed or unbounded.
    #[error("streamline controls must be finite, positive, and bounded")]
    InvalidOptions,
}

/// Integrates one deterministic polyline for every seed.
///
/// A seed outside the grid yields a one-point line rather than changing the
/// output cardinality. This keeps caller provenance aligned with seed order.
///
/// # Errors
///
/// Returns [`VectorFieldError`] for non-finite seeds or invalid controls.
pub fn integrate_streamlines(
    field: &VectorFieldGrid,
    seeds: &[[f32; 3]],
    options: StreamlineOptions,
) -> Result<Vec<Vec<[f32; 3]>>, VectorFieldError> {
    validate_options(options)?;
    if !seeds.iter().flatten().all(|value| value.is_finite()) {
        return Err(VectorFieldError::NonFiniteInput);
    }
    Ok(seeds
        .iter()
        .map(|seed| integrate_seed(field, *seed, options))
        .collect())
}

fn integrate_seed(
    field: &VectorFieldGrid,
    seed: [f32; 3],
    options: StreamlineOptions,
) -> Vec<[f32; 3]> {
    match options.direction {
        StreamlineDirection::Forward => integrate_direction(field, seed, options, 1.0),
        StreamlineDirection::Backward => integrate_direction(field, seed, options, -1.0),
        StreamlineDirection::Both => {
            let mut backward = integrate_direction(field, seed, options, -1.0);
            let forward = integrate_direction(field, seed, options, 1.0);
            backward.reverse();
            if !backward.is_empty() {
                backward.pop();
            }
            backward.extend(forward);
            backward
        }
    }
}

fn integrate_direction(
    field: &VectorFieldGrid,
    seed: [f32; 3],
    options: StreamlineOptions,
    sign: f32,
) -> Vec<[f32; 3]> {
    let mut points = Vec::with_capacity(options.max_steps.saturating_add(1));
    points.push(seed);
    let mut current = seed;
    let mut length = 0.0_f32;
    for _ in 0..options.max_steps {
        let Some(next) = rk4_step(field, current, options.step_size * sign, options.min_speed)
        else {
            break;
        };
        let distance = norm(sub(next, current));
        if !distance.is_finite()
            || distance <= f32::EPSILON
            || length + distance > options.max_length
        {
            break;
        }
        points.push(next);
        current = next;
        length += distance;
    }
    points
}

fn rk4_step(
    field: &VectorFieldGrid,
    point: [f32; 3],
    step: f32,
    min_speed: f32,
) -> Option<[f32; 3]> {
    let k1 = unit(field.sample(point)?, min_speed)?;
    let k2 = unit(field.sample(add_scaled(point, k1, step * 0.5))?, min_speed)?;
    let k3 = unit(field.sample(add_scaled(point, k2, step * 0.5))?, min_speed)?;
    let k4 = unit(field.sample(add_scaled(point, k3, step))?, min_speed)?;
    let mut next = point;
    for axis in 0..3 {
        next[axis] += step * (k1[axis] + 2.0 * k2[axis] + 2.0 * k3[axis] + k4[axis]) / 6.0;
    }
    field.sample(next).map(|_| next)
}

fn validate_options(options: StreamlineOptions) -> Result<(), VectorFieldError> {
    if options.step_size.is_finite()
        && options.step_size > 0.0
        && options.max_steps > 0
        && options.max_steps <= 1_000_000
        && options.max_length.is_finite()
        && options.max_length > 0.0
        && options.min_speed.is_finite()
        && options.min_speed >= 0.0
    {
        Ok(())
    } else {
        Err(VectorFieldError::InvalidOptions)
    }
}

fn unit(value: [f32; 3], minimum: f32) -> Option<[f32; 3]> {
    let magnitude = norm(value);
    if magnitude <= minimum || !magnitude.is_finite() {
        return None;
    }
    Some(value.map(|component| component / magnitude))
}

fn norm(value: [f32; 3]) -> f32 {
    value
        .into_iter()
        .map(|component| component * component)
        .sum::<f32>()
        .sqrt()
}

fn sub(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|axis| left[axis] - right[axis])
}

fn add_scaled(point: [f32; 3], direction: [f32; 3], scale: f32) -> [f32; 3] {
    std::array::from_fn(|axis| point[axis] + direction[axis] * scale)
}

fn axis_weight(fraction: f32, upper: usize) -> f32 {
    if upper == 0 { 1.0 - fraction } else { fraction }
}
