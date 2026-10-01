//! Resolution shells and amplitude normalization (`F` to `E`).
//!
//! The binner places reflections into shells of `1/d²`; the normalizer turns
//! amplitudes into normalized structure factors with Karle's approach: the mean
//! intensity of each shell, corrected for reflection multiplicity, smoothed over
//! neighbouring shells and interpolated between shell midpoints.

use crate::{CellTransform, SymmetrySet};

/// How shell boundaries are spaced.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinMethod {
    /// Shells hold approximately equal numbers of reflections.
    EqualCount,
    /// Equal steps in `1/d`.
    Dstar,
    /// Equal steps in `1/d²`.
    Dstar2,
    /// Equal steps in `1/d³`.
    Dstar3,
}

/// Why shells or normalizers could not be built.
#[derive(Clone, Copy, Debug, thiserror::Error, PartialEq, Eq)]
pub enum ReflectionBinningError {
    /// Fewer than one shell was requested.
    #[error("the number of resolution bins must be positive")]
    NoBins,
    /// No reflection was supplied.
    #[error("resolution bins need at least one reflection")]
    NoReflections,
    /// A spacing is negative or not finite.
    #[error("reciprocal spacings must be finite and non-negative")]
    InvalidSpacing,
    /// Columns of different lengths were supplied.
    #[error("reflection columns have different lengths")]
    LengthMismatch,
    /// The space group has no symmetry operations.
    #[error("normalization needs the complete set of symmetry operations")]
    NoOperations,
    /// More reflections or shells than the arithmetic can count exactly.
    #[error("more than four billion reflections or bins")]
    TooLarge,
    /// A reflection's indices cannot be transformed.
    #[error("a Miller index overflows under a symmetry operation")]
    IndexOverflow,
}

/// Shell boundaries over `1/d²`.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolutionBinner {
    limits: Vec<f64>,
    mids: Vec<f64>,
    min_inverse_d2: f64,
    max_inverse_d2: f64,
}

impl ResolutionBinner {
    /// Builds `bins` shells over the given `1/d²` values.
    ///
    /// Boundaries are first laid out for twice as many sub-shells, so that each
    /// shell has both an upper limit and a midpoint.
    ///
    /// # Errors
    ///
    /// Returns an error for zero bins, no reflections, or a spacing that is
    /// negative or not finite.
    pub fn new(
        method: BinMethod,
        bins: usize,
        inverse_d2: &[f64],
    ) -> Result<Self, ReflectionBinningError> {
        if bins == 0 {
            return Err(ReflectionBinningError::NoBins);
        }
        if inverse_d2.is_empty() {
            return Err(ReflectionBinningError::NoReflections);
        }
        if inverse_d2
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(ReflectionBinningError::InvalidSpacing);
        }
        let too_large = bins
            .checked_mul(2)
            .and_then(|doubled| u32::try_from(doubled).ok());
        if too_large.is_none() || u32::try_from(inverse_d2.len()).is_err() {
            return Err(ReflectionBinningError::TooLarge);
        }
        let sub_shells = bins * 2;
        let (minimum, maximum) = extent(inverse_d2);
        let mut boundaries = vec![0.0; sub_shells];
        match method {
            BinMethod::EqualCount => {
                let mut sorted = inverse_d2.to_vec();
                sorted.sort_by(f64::total_cmp);
                for index in 1..sub_shells {
                    boundaries[index - 1] = sorted[sorted.len() * index / sub_shells];
                }
            }
            BinMethod::Dstar2 => {
                let step = (maximum - minimum) / usize_as_f64(sub_shells);
                for index in 1..sub_shells {
                    boundaries[index - 1] = minimum + usize_as_f64(index) * step;
                }
            }
            BinMethod::Dstar => {
                let (low, high) = (minimum.sqrt(), maximum.sqrt());
                let step = (high - low) / usize_as_f64(sub_shells);
                for index in 1..sub_shells {
                    let value = low + usize_as_f64(index) * step;
                    boundaries[index - 1] = value * value;
                }
            }
            BinMethod::Dstar3 => {
                let (low, high) = (minimum * minimum.sqrt(), maximum * maximum.sqrt());
                let step = (high - low) / usize_as_f64(sub_shells);
                for index in 1..sub_shells {
                    let value = (low + usize_as_f64(index) * step).cbrt();
                    boundaries[index - 1] = value * value;
                }
            }
        }
        boundaries[sub_shells - 1] = f64::INFINITY;
        let mids = (0..bins).map(|index| boundaries[2 * index]).collect();
        let limits = (0..bins).map(|index| boundaries[2 * index + 1]).collect();
        Ok(Self {
            limits,
            mids,
            min_inverse_d2: minimum,
            max_inverse_d2: maximum,
        })
    }

    /// Builds shells over the spacings of `reflections` in `cell`.
    ///
    /// # Errors
    ///
    /// As [`Self::new`].
    pub fn for_reflections(
        method: BinMethod,
        bins: usize,
        cell: &CellTransform,
        reflections: &[[i32; 3]],
    ) -> Result<Self, ReflectionBinningError> {
        let spacings = reflections
            .iter()
            .map(|hkl| cell.reciprocal_spacing_squared(*hkl))
            .collect::<Vec<_>>();
        Self::new(method, bins, &spacings)
    }

    /// Number of shells.
    #[must_use]
    pub fn len(&self) -> usize {
        self.limits.len()
    }

    /// Whether there are no shells; construction never yields such a binner.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.limits.is_empty()
    }

    /// Upper `1/d²` limit of each shell; the last is infinite.
    #[must_use]
    pub fn limits(&self) -> &[f64] {
        &self.limits
    }

    /// Midpoint `1/d²` of each shell.
    #[must_use]
    pub fn midpoints(&self) -> &[f64] {
        &self.mids
    }

    /// The shell holding `inverse_d2`: the first whose upper limit is not below it.
    #[must_use]
    pub fn bin_index(&self, inverse_d2: f64) -> usize {
        self.limits
            .partition_point(|limit| *limit < inverse_d2)
            .min(self.limits.len() - 1)
    }

    /// Shell of every spacing.
    #[must_use]
    pub fn bin_indices(&self, inverse_d2: &[f64]) -> Vec<usize> {
        inverse_d2
            .iter()
            .map(|value| self.bin_index(*value))
            .collect()
    }

    /// Smallest `d` in a shell, in ångström; `None` outside the shells.
    #[must_use]
    pub fn d_min_of_bin(&self, bin: usize) -> Option<f64> {
        let limit = if bin + 1 == self.limits.len() {
            self.max_inverse_d2
        } else {
            *self.limits.get(bin)?
        };
        Some(limit.sqrt().recip())
    }

    /// Largest `d` in a shell, in ångström; `None` outside the shells.
    #[must_use]
    pub fn d_max_of_bin(&self, bin: usize) -> Option<f64> {
        let limit = if bin == 0 {
            self.min_inverse_d2
        } else {
            *self.limits.get(bin - 1)?
        };
        Some(limit.sqrt().recip())
    }
}

/// Per-reflection multipliers that turn amplitudes into normalized amplitudes.
///
/// Each shell's mean of `F²/ε` is smoothed with the kernel `[0.75, 1, 0.75]`
/// over its neighbours, and the root of it is interpolated linearly in `1/d²`
/// between shell midpoints. A reflection's multiplier is `1/(√ε · rms)`, so
/// `E = F · multiplier`. Reflections whose amplitude is not finite get `NaN` and
/// do not contribute to any shell.
///
/// `epsilon` is the number of space-group operations that fix the indices.
///
/// # Errors
///
/// Returns an error for mismatched columns, an empty operation set, or indices
/// that overflow under an operation.
pub fn amplitude_normalizers(
    cell: &CellTransform,
    symmetry: &SymmetrySet,
    reflections: &[[i32; 3]],
    amplitudes: &[f64],
    binner: &ResolutionBinner,
) -> Result<Vec<f64>, ReflectionBinningError> {
    if reflections.len() != amplitudes.len() {
        return Err(ReflectionBinningError::LengthMismatch);
    }
    if symmetry.operations().is_empty() {
        return Err(ReflectionBinningError::NoOperations);
    }
    let spacings = reflections
        .iter()
        .map(|hkl| cell.reciprocal_spacing_squared(*hkl))
        .collect::<Vec<_>>();
    let bins = binner.bin_indices(&spacings);
    let mut multipliers = vec![f64::NAN; reflections.len()];
    let mut counts = vec![0.0_f64; binner.len()];
    let mut sums = vec![0.0_f64; binner.len()];
    for (row, hkl) in reflections.iter().enumerate() {
        let amplitude = amplitudes[row];
        if !amplitude.is_finite() {
            continue;
        }
        let epsilon = fixing_operations(symmetry, *hkl)?;
        let inverse_epsilon = 1.0 / usize_as_f64(epsilon);
        multipliers[row] = inverse_epsilon.sqrt();
        counts[bins[row]] += 1.0;
        sums[bins[row]] += amplitude * amplitude * inverse_epsilon;
    }
    let rms = smoothed_rms(&counts, &sums);
    let mids = binner.midpoints();
    for (row, multiplier) in multipliers.iter_mut().enumerate() {
        let spacing = spacings[row];
        let mut bin = bins[row];
        let mut scale = rms[bin];
        if spacing > mids[0] && spacing < mids[mids.len() - 1] {
            if spacing > mids[bin] {
                bin += 1;
            }
            let (low, high) = (mids[bin - 1], mids[bin]);
            scale = rms[bin - 1] + (spacing - low) * (rms[bin] - rms[bin - 1]) / (high - low);
        }
        *multiplier /= scale;
    }
    Ok(multipliers)
}

fn fixing_operations(
    symmetry: &SymmetrySet,
    hkl: [i32; 3],
) -> Result<usize, ReflectionBinningError> {
    let mut epsilon = 0;
    for operation in symmetry.operations() {
        let mapped = operation
            .apply_to_hkl(hkl)
            .map_err(|_| ReflectionBinningError::IndexOverflow)?;
        if mapped == hkl {
            epsilon += 1;
        }
    }
    Ok(epsilon)
}

fn smoothed_rms(counts: &[f64], sums: &[f64]) -> Vec<f64> {
    const KERNEL: f64 = 0.75;
    let last = counts.len() - 1;
    (0..counts.len())
        .map(|bin| {
            let mut count = counts[bin];
            let mut sum = sums[bin];
            if bin > 0 {
                count += KERNEL * counts[bin - 1];
                sum += KERNEL * sums[bin - 1];
            }
            if bin < last {
                count += KERNEL * counts[bin + 1];
                sum += KERNEL * sums[bin + 1];
            }
            (sum / count).sqrt()
        })
        .collect()
}

fn extent(values: &[f64]) -> (f64, f64) {
    values
        .iter()
        .fold((values[0], values[0]), |(low, high), value| {
            (low.min(*value), high.max(*value))
        })
}

/// A count as `f64`; callers have checked it fits in `u32`, so it is exact.
fn usize_as_f64(value: usize) -> f64 {
    match u32::try_from(value) {
        Ok(value) => f64::from(value),
        Err(_) => f64::from(u32::MAX),
    }
}

#[cfg(test)]
#[path = "reflection_binning_tests.rs"]
mod tests;
