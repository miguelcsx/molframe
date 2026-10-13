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
/// use molframe_geom::rmsf;
///
/// // One atom that sits at the origin in one frame and two units along x in
/// // the next: it strays one unit either side of its mean.
/// let a = [[0.0, 0.0, 0.0]];
/// let b = [[2.0, 0.0, 0.0]];
/// let values = rmsf(&[&a, &b])?;
/// assert!((values[0] - 1.0).abs() < 1e-9);
/// # Ok::<(), molframe_geom::FluctuationError>(())
/// ```
pub fn rmsf(frames: &[&[[f32; 3]]]) -> Result<Vec<f64>, FluctuationError> {
    let Some((first, rest)) = frames.split_first() else {
        return Err(FluctuationError::NoFrames);
    };

    let atom_count = first.len();
    if rest.iter().any(|frame| frame.len() != atom_count) {
        return Err(FluctuationError::RaggedFrames);
    }

    let mut accumulator = RmsfAccumulator::new(atom_count);
    for frame in frames {
        accumulator.push(frame)?;
    }
    accumulator.finish()
}

/// Incremental per-atom RMSF over frames supplied one at a time.
///
/// Holds a mean and a summed squared deviation per atom, so memory is
/// `O(atoms)` for any number of frames. Frames are folded with the same
/// vectorised Welford kernel as [`rmsf`], so both give identical results.
///
/// # Examples
///
/// ```
/// use molframe_geom::RmsfAccumulator;
///
/// let mut accumulator = RmsfAccumulator::new(1);
/// accumulator.push(&[[0.0, 0.0, 0.0]])?;
/// accumulator.push(&[[2.0, 0.0, 0.0]])?;
/// let values = accumulator.finish()?;
/// assert!((values[0] - 1.0).abs() < 1e-9);
/// # Ok::<(), molframe_geom::FluctuationError>(())
/// ```
#[derive(Clone, Debug)]
pub struct RmsfAccumulator {
    mean: Vec<[f64; 3]>,
    summed_square: Vec<f64>,
    seen: u64,
}

impl RmsfAccumulator {
    /// Starts an empty accumulator for `atoms` atoms.
    #[must_use]
    pub fn new(atoms: usize) -> Self {
        Self {
            mean: vec![[0.0; 3]; atoms],
            summed_square: vec![0.0; atoms],
            seen: 0,
        }
    }

    /// Folds one frame into the running mean and deviation.
    ///
    /// # Errors
    ///
    /// Returns [`FluctuationError::RaggedFrames`] when the frame length differs
    /// from the atom count and [`FluctuationError::TooManyFrames`] when the
    /// frame count leaves the exactly representable range.
    pub fn push(&mut self, frame: &[[f32; 3]]) -> Result<(), FluctuationError> {
        if frame.len() != self.mean.len() {
            return Err(FluctuationError::RaggedFrames);
        }
        let seen = self.seen + 1;
        let inverse_seen = exact_u64(seen)
            .ok_or(FluctuationError::TooManyFrames)?
            .recip();
        self.seen = seen;
        crate::simd::rmsf_update(frame, &mut self.mean, &mut self.summed_square, inverse_seen);
        Ok(())
    }

    /// Frames folded so far.
    #[must_use]
    pub const fn frames(&self) -> u64 {
        self.seen
    }

    /// Returns the RMSF of every atom.
    ///
    /// # Errors
    ///
    /// Returns [`FluctuationError::NoFrames`] when no frame was pushed.
    pub fn finish(self) -> Result<Vec<f64>, FluctuationError> {
        if self.seen == 0 {
            return Err(FluctuationError::NoFrames);
        }
        let inverse_frames = exact_u64(self.seen)
            .ok_or(FluctuationError::TooManyFrames)?
            .recip();
        let mut values = self.summed_square;
        for value in &mut values {
            let mean_square = *value * inverse_frames;
            *value = if mean_square <= 0.0 {
                0.0
            } else {
                mean_square.sqrt()
            };
        }
        Ok(values)
    }
}

#[cfg(test)]
#[path = "rmsf_tests.rs"]
mod tests;
