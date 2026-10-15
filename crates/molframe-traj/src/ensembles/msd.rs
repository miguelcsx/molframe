//! Window-averaged mean-squared displacement over caller-un-wrapped coordinates.

use crate::frame_view::{FrameSource, FrameView};
use crate::numeric::f64_from_u64;

/// MSD at one frame lag.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeanSquaredDisplacement {
    /// Difference in frame indices.
    pub lag: usize,
    /// Number of atom-origin observations averaged.
    pub observations: u64,
    /// Mean squared Cartesian displacement.
    pub value: f64,
}

/// Why an MSD series could not be evaluated.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum MsdError {
    /// At least one non-empty frame is required.
    #[error("at least one non-empty frame is required")]
    Empty,
    /// Frames must contain the same number of atoms.
    #[error("all frames must contain the same number of atoms")]
    DimensionMismatch,
    /// A selected atom is outside the frames.
    #[error("selected atom {0} is outside the frames")]
    AtomOutOfBounds(usize),
    /// Coordinates must be finite.
    #[error("coordinates must be finite")]
    NonFinite,
    /// The maximum lag must be smaller than the frame count.
    #[error("maximum lag must be smaller than the frame count")]
    InvalidLag,
    /// The observation count cannot be represented by the result schema.
    #[error("MSD observation count exceeds exact numeric representation")]
    ObservationLimit,
}

/// Computes window-averaged MSD for explicit atom indices and lags `0..=maximum_lag`.
///
/// Coordinates are used exactly as supplied. Periodic trajectories must be
/// un-wrapped through an explicit trajectory transform before this kernel is
/// called; the function never guesses molecular images.
///
/// # Errors
///
/// Returns [`MsdError`] for empty, ragged, non-finite, out-of-range, or invalid input.
pub fn mean_squared_displacement(
    frames: &[Vec<[f32; 3]>],
    atoms: &[usize],
    maximum_lag: usize,
) -> Result<Vec<MeanSquaredDisplacement>, MsdError> {
    mean_squared_displacement_source(frames, atoms, maximum_lag)
}

/// Computes window-averaged MSD directly from borrowed contiguous frames.
///
/// Empty `atoms` selects all atoms without materialising an index vector.
///
/// # Errors
///
/// Returns [`MsdError`] for empty, malformed, non-finite, out-of-range, or
/// invalid input.
pub fn mean_squared_displacement_view(
    frames: FrameView<'_>,
    atoms: &[usize],
    maximum_lag: usize,
) -> Result<Vec<MeanSquaredDisplacement>, MsdError> {
    mean_squared_displacement_source(&frames, atoms, maximum_lag)
}

fn mean_squared_displacement_source<S: FrameSource + ?Sized>(
    frames: &S,
    atoms: &[usize],
    maximum_lag: usize,
) -> Result<Vec<MeanSquaredDisplacement>, MsdError> {
    let atom_count = validate(frames, atoms, maximum_lag)?;
    (0..=maximum_lag)
        .map(|lag| lag_value(frames, atoms, atom_count, lag))
        .collect()
}

fn validate<S: FrameSource + ?Sized>(
    frames: &S,
    atoms: &[usize],
    maximum_lag: usize,
) -> Result<usize, MsdError> {
    if frames.frame_count() == 0 {
        return Err(MsdError::Empty);
    }
    let atom_count = frames.atom_count();
    if atom_count == 0 {
        return Err(MsdError::Empty);
    }
    for index in 0..frames.frame_count() {
        let Some(frame) = frames.frame(index) else {
            return Err(MsdError::DimensionMismatch);
        };
        if frame.len() != atom_count {
            return Err(MsdError::DimensionMismatch);
        }
        if frame.iter().flatten().any(|value| !value.is_finite()) {
            return Err(MsdError::NonFinite);
        }
    }
    if let Some(atom) = atoms.iter().copied().find(|atom| *atom >= atom_count) {
        return Err(MsdError::AtomOutOfBounds(atom));
    }
    if maximum_lag >= frames.frame_count() {
        return Err(MsdError::InvalidLag);
    }
    Ok(atom_count)
}

fn lag_value<S: FrameSource + ?Sized>(
    frames: &S,
    atoms: &[usize],
    atom_count: usize,
    lag: usize,
) -> Result<MeanSquaredDisplacement, MsdError> {
    let mut sum = 0.0;
    let mut observations = 0_u64;
    for origin in 0..frames.frame_count() - lag {
        let left = frames.frame(origin).ok_or(MsdError::DimensionMismatch)?;
        let right = frames
            .frame(origin + lag)
            .ok_or(MsdError::DimensionMismatch)?;
        if atoms.is_empty() {
            for atom in 0..atom_count {
                sum += squared_distance(left[atom], right[atom]);
                observations = observations
                    .checked_add(1)
                    .ok_or(MsdError::ObservationLimit)?;
            }
        } else {
            for &atom in atoms {
                sum += squared_distance(left[atom], right[atom]);
                observations = observations
                    .checked_add(1)
                    .ok_or(MsdError::ObservationLimit)?;
            }
        }
    }
    let divisor = f64_from_u64(observations).ok_or(MsdError::ObservationLimit)?;
    Ok(MeanSquaredDisplacement {
        lag,
        observations,
        value: sum / divisor,
    })
}

fn squared_distance(left: [f32; 3], right: [f32; 3]) -> f64 {
    left.into_iter()
        .zip(right)
        .map(|(left, right)| {
            let delta = f64::from(right - left);
            delta * delta
        })
        .sum()
}

#[cfg(test)]
#[path = "msd_tests.rs"]
mod tests;
