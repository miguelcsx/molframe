//! Positional fluctuation of atoms across an ensemble of frames.
//!
//! The root-mean-square fluctuation of an atom is how far it strays from its
//! own average position over a set of frames that are already in
//! correspondence. It is not a deviation between two structures — that is
//! [`crate::rmsd`] — and it says nothing about whether the frames were aligned
//! first. Aligning the frames is the caller's decision, because whether the
//! overall tumbling of a molecule counts as fluctuation depends on the
//! question being asked.
//!
//! The mean and the summed squared deviation are accumulated together in one
//! pass with the Welford recurrence, so a long trajectory never forms the
//! difference of two large second moments. Cost is `O(frames · atoms)` time and
//! one output buffer of `atoms` values.

use crate::numeric::exact_u64;

/// Why a fluctuation could not be computed.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FluctuationError {
    /// No frames were supplied, so there is nothing to average over.
    NoFrames,
    /// Two frames disagree on how many atoms they hold, so no atom has a stable
    /// identity across the ensemble.
    RaggedFrames,
    /// The frame count could not be represented by the numeric accumulator.
    TooManyFrames,
}

/// Per-atom root-mean-square fluctuation about each atom's mean position.
///
/// Every frame must list the same atoms in the same order; the returned vector
/// holds one value per atom in that order. A single frame yields all zeros,
/// because one observation cannot depart from its own mean.
///
/// Runs in `O(frames · atoms)` time and allocates the result plus a mean buffer.
///
/// # Errors
///
/// Returns [`FluctuationError::NoFrames`] when `frames` is empty, and
/// [`FluctuationError::RaggedFrames`] when the frames differ in length.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::rmsf;
///
/// // One atom that sits at the origin in one frame and two units along x in
/// // the next: it strays one unit either side of its mean.
/// let a = [[0.0, 0.0, 0.0]];
/// let b = [[2.0, 0.0, 0.0]];
/// let values = rmsf(&[&a, &b])?;
/// assert!((values[0] - 1.0).abs() < 1e-9);
/// # Ok::<(), pdbiox_geom::FluctuationError>(())
/// ```
pub fn rmsf(frames: &[&[[f32; 3]]]) -> Result<Vec<f64>, FluctuationError> {
    let Some((first, rest)) = frames.split_first() else {
        return Err(FluctuationError::NoFrames);
    };

    let atom_count = first.len();
    if rest.iter().any(|frame| frame.len() != atom_count) {
        return Err(FluctuationError::RaggedFrames);
    }

    let mut mean = vec![[0.0f64; 3]; atom_count];
    let mut summed_square = vec![0.0f64; atom_count];
    let mut seen = 0u64;

    for frame in frames {
        seen += 1;
        let inverse_seen = exact_u64(seen)
            .ok_or(FluctuationError::TooManyFrames)?
            .recip();
        for (atom, &position) in frame.iter().enumerate() {
            let point = [
                f64::from(position[0]),
                f64::from(position[1]),
                f64::from(position[2]),
            ];
            let before = [
                point[0] - mean[atom][0],
                point[1] - mean[atom][1],
                point[2] - mean[atom][2],
            ];
            mean[atom][0] += before[0] * inverse_seen;
            mean[atom][1] += before[1] * inverse_seen;
            mean[atom][2] += before[2] * inverse_seen;
            let after = [
                point[0] - mean[atom][0],
                point[1] - mean[atom][1],
                point[2] - mean[atom][2],
            ];
            summed_square[atom] +=
                before[0] * after[0] + before[1] * after[1] + before[2] * after[2];
        }
    }

    let inverse_frames = exact_u64(seen)
        .ok_or(FluctuationError::TooManyFrames)?
        .recip();
    for value in &mut summed_square {
        let mean_square = *value * inverse_frames;
        *value = if mean_square <= 0.0 {
            0.0
        } else {
            mean_square.sqrt()
        };
    }

    Ok(summed_square)
}

#[cfg(test)]
#[path = "rmsf_tests.rs"]
mod tests;
