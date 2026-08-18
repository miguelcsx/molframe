//! Topology-stable interpolation between caller-selected trajectory frames.
//!
//! Every sample costs `O(atoms)` and allocates exactly one output coordinate
//! array. The cubic path is centripetal Catmull-Rom, reducing loops and
//! overshoot around unevenly spaced conformations.

#[cfg(test)]
#[path = "interpolation_tests.rs"]
mod tests;

/// Interpolation applied between the two central frames.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum TrajectoryInterpolation {
    /// Straight Cartesian interpolation between the central pair.
    #[default]
    Linear,
    /// Centripetal Catmull-Rom using all four frames.
    CentripetalCatmullRom,
}

/// Why trajectory coordinates could not be interpolated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq)]
pub enum TrajectoryInterpolationError {
    /// All frames must share one non-empty topology.
    #[error("trajectory interpolation requires four non-empty atom-aligned frames")]
    TopologyMismatch,
    /// Coordinates and the interpolation fraction must be finite.
    #[error("trajectory interpolation inputs must be finite")]
    NonFiniteInput,
    /// The interpolation fraction lies outside the central interval.
    #[error("trajectory interpolation fraction must lie in [0, 1]")]
    FractionOutOfRange,
}

/// Interpolates the central pair of four topology-aligned coordinate frames.
///
/// Linear interpolation reads only `frames[1]` and `frames[2]`; the four-frame
/// signature lets callers change modes without changing their working set.
///
/// # Errors
///
/// Returns [`TrajectoryInterpolationError`] for malformed frames or fraction.
pub fn interpolate_trajectory_frames(
    frames: [&[[f32; 3]]; 4],
    fraction: f32,
    interpolation: TrajectoryInterpolation,
) -> Result<Vec<[f32; 3]>, TrajectoryInterpolationError> {
    validate(frames, fraction)?;
    Ok((0..frames[0].len())
        .map(|atom| match interpolation {
            TrajectoryInterpolation::Linear => linear(frames[1][atom], frames[2][atom], fraction),
            TrajectoryInterpolation::CentripetalCatmullRom => centripetal(
                frames[0][atom],
                frames[1][atom],
                frames[2][atom],
                frames[3][atom],
                fraction,
            ),
        })
        .collect())
}

fn validate(frames: [&[[f32; 3]]; 4], fraction: f32) -> Result<(), TrajectoryInterpolationError> {
    let atom_count = frames[0].len();
    if atom_count == 0 || frames.iter().any(|frame| frame.len() != atom_count) {
        return Err(TrajectoryInterpolationError::TopologyMismatch);
    }
    if !fraction.is_finite()
        || frames
            .iter()
            .flat_map(|frame| frame.iter().flatten())
            .any(|value| !value.is_finite())
    {
        return Err(TrajectoryInterpolationError::NonFiniteInput);
    }
    if !(0.0..=1.0).contains(&fraction) {
        return Err(TrajectoryInterpolationError::FractionOutOfRange);
    }
    Ok(())
}

fn centripetal(p0: [f32; 3], p1: [f32; 3], p2: [f32; 3], p3: [f32; 3], fraction: f32) -> [f32; 3] {
    let t0 = 0.0_f32;
    let t1 = t0 + knot_distance(p0, p1);
    let t2 = t1 + knot_distance(p1, p2);
    let t3 = t2 + knot_distance(p2, p3);
    if t1 - t0 <= f32::EPSILON || t2 - t1 <= f32::EPSILON || t3 - t2 <= f32::EPSILON {
        return linear(p1, p2, fraction);
    }
    let t = t1 + (t2 - t1) * fraction;
    let a1 = timed_linear(p0, p1, t0, t1, t);
    let a2 = timed_linear(p1, p2, t1, t2, t);
    let a3 = timed_linear(p2, p3, t2, t3, t);
    let b1 = timed_linear(a1, a2, t0, t2, t);
    let b2 = timed_linear(a2, a3, t1, t3, t);
    timed_linear(b1, b2, t1, t2, t)
}

fn knot_distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    squared_distance(left, right).sqrt().sqrt()
}

fn timed_linear(
    left: [f32; 3],
    right: [f32; 3],
    left_time: f32,
    right_time: f32,
    time: f32,
) -> [f32; 3] {
    let fraction = (time - left_time) / (right_time - left_time);
    linear(left, right, fraction)
}

fn linear(left: [f32; 3], right: [f32; 3], fraction: f32) -> [f32; 3] {
    std::array::from_fn(|axis| left[axis] + (right[axis] - left[axis]) * fraction)
}

fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f32 {
    (0..3)
        .map(|axis| {
            let delta = right[axis] - left[axis];
            delta * delta
        })
        .sum()
}
