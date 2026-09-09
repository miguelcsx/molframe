//! Native pore-radius profiles from explicit atoms, radii, and search geometry.
//!
//! The public result remains the exhaustive circular-grid optimum. Large
//! searches build a balanced bounding-volume tree in `O(N log N)` time and
//! `O(N)` memory, then prune atoms whose conservative surface-clearance bound
//! cannot improve a candidate. Query cost is `O(log N)` when bounds separate
//! the atoms and `O(N)` in the degenerate worst case.

use crate::numeric::{f32_to_usize, usize_to_f32};

#[path = "pore/index.rs"]
mod spatial_index;
use spatial_index::{ClearanceIndex, IndexBuildError};

const INDEX_MIN_DISTANCE_EVALUATIONS: usize = 65_536;

/// Search cylinder and sampling resolution for a pore profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PoreProfileOptions {
    /// Point on the pore axis at axial coordinate zero.
    pub axis_origin: [f32; 3],
    /// Axis direction; normalized internally.
    pub axis_direction: [f32; 3],
    /// First axial coordinate to sample.
    pub start: f32,
    /// Last axial coordinate to sample, inclusive when more than one sample is requested.
    pub end: f32,
    /// Number of axial slices.
    pub samples: usize,
    /// Maximum transverse displacement searched around the axis.
    pub search_radius: f32,
    /// Transverse Cartesian grid spacing.
    pub grid_spacing: f32,
    /// Radius of a probe subtracted from the geometric clearance.
    pub probe_radius: f32,
    /// Hard ceiling for output and spatial-search workspace.
    pub memory_limit_bytes: usize,
}

impl PoreProfileOptions {
    /// Default ceiling, chosen to keep the operation-owned working set near 100 MB.
    pub const DEFAULT_MEMORY_LIMIT_BYTES: usize = 100_000_000;

    /// Creates an exhaustive-grid profile request under the default 100 MB ceiling.
    #[must_use]
    pub const fn new(
        axis: ([f32; 3], [f32; 3]),
        start: f32,
        end: f32,
        samples: usize,
        search_radius: f32,
        grid_spacing: f32,
        probe_radius: f32,
    ) -> Self {
        Self {
            axis_origin: axis.0,
            axis_direction: axis.1,
            start,
            end,
            samples,
            search_radius,
            grid_spacing,
            probe_radius,
            memory_limit_bytes: Self::DEFAULT_MEMORY_LIMIT_BYTES,
        }
    }

    /// Sets the operation-owned memory ceiling.
    #[must_use]
    pub const fn with_memory_limit(mut self, bytes: usize) -> Self {
        self.memory_limit_bytes = bytes;
        self
    }
}

/// Best probe centre and clearance at one axial slice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PoreSample {
    /// Coordinate along the requested axis.
    pub axial_coordinate: f32,
    /// Cartesian centre with maximum sampled clearance.
    pub centre: [f32; 3],
    /// Maximum non-negative probe-centre clearance.
    pub radius: f32,
}

/// Why a pore profile could not be computed.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq)]
pub enum PoreError {
    /// Positions and radii must be atom-aligned.
    #[error("expected {positions} radii, found {radii}")]
    LengthMismatch {
        /// Number of positions.
        positions: usize,
        /// Number of radii.
        radii: usize,
    },
    /// At least one atom is needed to bound a pore.
    #[error("at least one atom is required to bound a pore")]
    EmptyAtoms,
    /// Atom positions and radii must be finite and radii non-negative.
    #[error("atom positions and radii must be finite and radii non-negative")]
    InvalidAtoms,
    /// Axis and sampling options are invalid.
    #[error(
        "axis and sampling parameters must be finite, positive where required, and non-degenerate"
    )]
    InvalidOptions,
    /// The requested transverse grid cannot be represented.
    #[error("the transverse pore-search grid is too large")]
    GridTooLarge,
    /// Output and exact spatial-search storage exceed the caller's ceiling.
    #[error("pore-profile storage requires {required} bytes, over the {limit} byte limit")]
    MemoryLimit {
        /// Required operation-owned bytes, or `usize::MAX` after overflow or allocation refusal.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
}

/// Computes a deterministic pore-radius profile without an external executable.
///
/// Each axial slice searches the complete explicit circular transverse grid and
/// reports the point maximizing the minimum atom-surface clearance. Spatial
/// pruning changes only which atoms can be proven irrelevant, not the sampled
/// points or their stable first-maximum tie break. Decreasing `grid_spacing`
/// converges on the continuous plane optimization.
///
/// # Errors
///
/// Returns [`PoreError`] for invalid atoms, axis, bounds, or sampling resolution.
pub fn pore_profile(
    positions: &[[f32; 3]],
    radii: &[f32],
    options: PoreProfileOptions,
) -> Result<Vec<PoreSample>, PoreError> {
    validate(positions, radii, options)?;
    let axis = normalize(options.axis_direction);
    let first_basis = perpendicular(axis);
    let second_basis = cross(axis, first_basis);
    let transverse_steps = (2.0 * options.search_radius / options.grid_spacing).ceil();
    let transverse_steps = f32_to_usize(transverse_steps).ok_or(PoreError::GridTooLarge)?;
    let side = transverse_steps
        .checked_add(1)
        .ok_or(PoreError::GridTooLarge)?;
    let candidate_upper_bound = side
        .checked_mul(side)
        .and_then(|count| count.checked_mul(options.samples))
        .ok_or(PoreError::GridTooLarge)?;
    let output_bytes =
        options
            .samples
            .checked_mul(size_of::<PoreSample>())
            .ok_or(PoreError::MemoryLimit {
                required: usize::MAX,
                limit: options.memory_limit_bytes,
            })?;
    if output_bytes > options.memory_limit_bytes {
        return Err(PoreError::MemoryLimit {
            required: output_bytes,
            limit: options.memory_limit_bytes,
        });
    }
    let distance_evaluations = candidate_upper_bound.saturating_mul(positions.len());
    let should_index = positions.len() > spatial_index::LEAF_ATOMS
        && distance_evaluations >= INDEX_MIN_DISTANCE_EVALUATIONS;
    let clearance_index = if should_index {
        let index_bytes =
            ClearanceIndex::required_bytes(positions.len()).ok_or(PoreError::MemoryLimit {
                required: usize::MAX,
                limit: options.memory_limit_bytes,
            })?;
        let required = output_bytes
            .checked_add(index_bytes)
            .ok_or(PoreError::MemoryLimit {
                required: usize::MAX,
                limit: options.memory_limit_bytes,
            })?;
        if required > options.memory_limit_bytes {
            return Err(PoreError::MemoryLimit {
                required,
                limit: options.memory_limit_bytes,
            });
        }
        Some(match ClearanceIndex::build(positions, radii) {
            Ok(index) => index,
            Err(IndexBuildError::Dimension) => return Err(PoreError::GridTooLarge),
            Err(IndexBuildError::Allocation) => {
                return Err(PoreError::MemoryLimit {
                    required: usize::MAX,
                    limit: options.memory_limit_bytes,
                });
            }
        })
    } else {
        None
    };

    let mut search = PoreSearch {
        positions,
        radii,
        options,
        axis,
        first_basis,
        second_basis,
        transverse_steps,
        clearance_index,
        nearest_atom: 0,
    };
    let mut profile = Vec::new();
    profile
        .try_reserve_exact(options.samples)
        .map_err(|_| PoreError::MemoryLimit {
            required: usize::MAX,
            limit: options.memory_limit_bytes,
        })?;
    for sample in 0..options.samples {
        profile.push(search.sample_slice(sample));
    }
    Ok(profile)
}

struct PoreSearch<'a> {
    positions: &'a [[f32; 3]],
    radii: &'a [f32],
    options: PoreProfileOptions,
    axis: [f32; 3],
    first_basis: [f32; 3],
    second_basis: [f32; 3],
    transverse_steps: usize,
    clearance_index: Option<ClearanceIndex<'a>>,
    nearest_atom: u32,
}

impl PoreSearch<'_> {
    fn sample_slice(&mut self, sample: usize) -> PoreSample {
        let fraction = if self.options.samples == 1 {
            0.0
        } else {
            usize_to_f32(sample) / usize_to_f32(self.options.samples - 1)
        };
        let axial_coordinate =
            self.options.start + fraction * (self.options.end - self.options.start);
        let axis_centre = add(self.options.axis_origin, scale(self.axis, axial_coordinate));
        let mut best_centre = axis_centre;
        let mut best_clearance = f32::NEG_INFINITY;

        for first_index in 0..=self.transverse_steps {
            let first = transverse_coordinate(
                first_index,
                self.transverse_steps,
                self.options.search_radius,
            );
            for second_index in 0..=self.transverse_steps {
                let second = transverse_coordinate(
                    second_index,
                    self.transverse_steps,
                    self.options.search_radius,
                );
                if first.mul_add(first, second * second) > self.options.search_radius.powi(2) {
                    continue;
                }
                let candidate = add(
                    add(axis_centre, scale(self.first_basis, first)),
                    scale(self.second_basis, second),
                );
                let clearance = self.minimum_clearance(candidate) - self.options.probe_radius;
                if clearance > best_clearance {
                    best_clearance = clearance;
                    best_centre = candidate;
                }
            }
        }

        PoreSample {
            axial_coordinate,
            centre: best_centre,
            radius: best_clearance.max(0.0),
        }
    }

    fn minimum_clearance(&mut self, point: [f32; 3]) -> f32 {
        let Some(index) = self.clearance_index.as_ref() else {
            return minimum_clearance(point, self.positions, self.radii);
        };
        let (clearance, nearest_atom) = index.minimum_clearance(point, self.nearest_atom);
        self.nearest_atom = nearest_atom;
        clearance
    }
}

fn validate(
    positions: &[[f32; 3]],
    radii: &[f32],
    options: PoreProfileOptions,
) -> Result<(), PoreError> {
    if positions.len() != radii.len() {
        return Err(PoreError::LengthMismatch {
            positions: positions.len(),
            radii: radii.len(),
        });
    }
    if positions.is_empty() {
        return Err(PoreError::EmptyAtoms);
    }
    if !positions.iter().flatten().all(|value| value.is_finite())
        || !radii
            .iter()
            .all(|radius| radius.is_finite() && *radius >= 0.0)
    {
        return Err(PoreError::InvalidAtoms);
    }
    let axis_norm = norm(options.axis_direction);
    if !options.axis_origin.iter().all(|value| value.is_finite())
        || !axis_norm.is_finite()
        || axis_norm <= 0.0
        || !options.start.is_finite()
        || !options.end.is_finite()
        || options.samples == 0
        || !options.search_radius.is_finite()
        || options.search_radius <= 0.0
        || !options.grid_spacing.is_finite()
        || options.grid_spacing <= 0.0
        || !options.probe_radius.is_finite()
        || options.probe_radius < 0.0
        || options.memory_limit_bytes == 0
    {
        return Err(PoreError::InvalidOptions);
    }
    Ok(())
}

fn transverse_coordinate(index: usize, steps: usize, radius: f32) -> f32 {
    if steps == 0 {
        0.0
    } else {
        -radius + 2.0 * radius * usize_to_f32(index) / usize_to_f32(steps)
    }
}

fn minimum_clearance(point: [f32; 3], positions: &[[f32; 3]], radii: &[f32]) -> f32 {
    positions
        .iter()
        .zip(radii)
        .map(|(position, radius)| norm(subtract(point, *position)) - radius)
        .fold(f32::INFINITY, f32::min)
}

fn perpendicular(axis: [f32; 3]) -> [f32; 3] {
    let absolute = axis.map(f32::abs);
    let reference = if absolute[0] <= absolute[1] && absolute[0] <= absolute[2] {
        [1.0, 0.0, 0.0]
    } else if absolute[1] <= absolute[2] {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    normalize(cross(axis, reference))
}

fn normalize(vector: [f32; 3]) -> [f32; 3] {
    scale(vector, norm(vector).recip())
}

fn norm(vector: [f32; 3]) -> f32 {
    vector
        .into_iter()
        .map(|value| value * value)
        .sum::<f32>()
        .sqrt()
}

fn cross(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

fn add(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|axis| left[axis] + right[axis])
}

fn subtract(left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
    std::array::from_fn(|axis| left[axis] - right[axis])
}

fn scale(vector: [f32; 3], factor: f32) -> [f32; 3] {
    vector.map(|value| value * factor)
}

#[cfg(test)]
#[path = "pore_tests.rs"]
mod tests;
