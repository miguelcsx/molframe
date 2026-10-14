//! Detection of the decimal precision a stream of confidence values carries,
//! so storage can be as small as the source allows without changing any value.
//!
//! One pass over the values keeps, for every candidate number of decimals, a
//! flag for "every value so far decodes back to exactly itself" and the
//! largest scaled integer seen. Choosing a layout afterwards costs nothing.

use super::half::{f16_bits_to_f32, f32_to_f16_bits};
use molframe_core::numeric::{f32_eq, f64_eq, f64_to_f32, f64_to_u64, u64_to_f64};

/// Largest number of decimals probed; confidence files carry at most two.
pub(crate) const MAX_DECIMALS: usize = 4;

const POWERS: [u32; MAX_DECIMALS + 1] = [1, 10, 100, 1_000, 10_000];

/// Scaled integer a value maps to at a given denominator, when it exists.
pub(crate) fn scaled(value: f32, denominator: u32) -> Option<u64> {
    let scaled = (f64::from(value) * f64::from(denominator)).round();
    if scaled < 0.0 || scaled > f64::from(u32::MAX) {
        return None;
    }
    Some(f64_to_u64(scaled))
}

/// The single-precision value a scaled integer stands for.
pub(crate) fn unscaled(scaled: u64, denominator: u32) -> f32 {
    f64_to_f32(u64_to_f64(scaled) / f64::from(denominator))
}

/// Accumulates what storage widths can hold a value stream exactly.
#[derive(Clone, Debug)]
pub struct PrecisionProbe {
    count: u64,
    maximum: f32,
    fixed_exact: [bool; MAX_DECIMALS + 1],
    fixed_maximum: [u64; MAX_DECIMALS + 1],
    half_exact: bool,
    single_exact: bool,
    non_finite: bool,
}

impl Default for PrecisionProbe {
    fn default() -> Self {
        Self::new()
    }
}

impl PrecisionProbe {
    /// Starts an empty probe.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            count: 0,
            maximum: 0.0,
            fixed_exact: [true; MAX_DECIMALS + 1],
            fixed_maximum: [0; MAX_DECIMALS + 1],
            half_exact: true,
            single_exact: true,
            non_finite: false,
        }
    }

    /// Folds one value into the probe.
    pub fn observe(&mut self, value: f32) {
        self.count += 1;
        if !value.is_finite() {
            self.non_finite = true;
            return;
        }
        if value > self.maximum {
            self.maximum = value;
        }
        if self.half_exact && !f32_eq(f16_bits_to_f32(f32_to_f16_bits(value)), value) {
            self.half_exact = false;
        }
        for (slot, denominator) in POWERS.iter().enumerate() {
            if !self.fixed_exact[slot] {
                continue;
            }
            match scaled(value, *denominator) {
                Some(integer) if f32_eq(unscaled(integer, *denominator), value) => {
                    self.fixed_maximum[slot] = self.fixed_maximum[slot].max(integer);
                }
                _ => self.fixed_exact[slot] = false,
            }
        }
    }

    /// Folds one double-precision value into the probe.
    ///
    /// Fixed-point exactness is judged in binary64, so a value such as `12.34`
    /// parsed from text is exact at two decimals even though no `f32` equals it.
    pub fn observe_f64(&mut self, value: f64) {
        self.count += 1;
        if !value.is_finite() {
            self.non_finite = true;
            return;
        }
        let single = f64_to_f32(value);
        if single > self.maximum {
            self.maximum = single;
        }
        if !f64_eq(f64::from(single), value) {
            self.single_exact = false;
        }
        if self.half_exact && !f64_eq(f64::from(f16_bits_to_f32(f32_to_f16_bits(single))), value) {
            self.half_exact = false;
        }
        for (slot, denominator) in POWERS.iter().enumerate() {
            if !self.fixed_exact[slot] {
                continue;
            }
            let integer = (value * f64::from(*denominator)).round();
            if integer >= 0.0
                && integer <= f64::from(u32::MAX)
                && f64_eq(integer / f64::from(*denominator), value)
            {
                self.fixed_maximum[slot] = self.fixed_maximum[slot].max(f64_to_u64(integer));
            } else {
                self.fixed_exact[slot] = false;
            }
        }
    }

    /// Whether every value observed through [`Self::observe_f64`] fits an `f32` exactly.
    #[must_use]
    pub const fn single_is_exact(&self) -> bool {
        self.single_exact
    }

    /// Values observed so far.
    #[must_use]
    pub const fn count(&self) -> u64 {
        self.count
    }

    /// Whether a non-finite value was seen.
    #[must_use]
    pub const fn saw_non_finite(&self) -> bool {
        self.non_finite
    }

    /// Largest finite value seen, or zero for an empty probe.
    #[must_use]
    pub const fn maximum(&self) -> f32 {
        self.maximum
    }

    /// Whether every value survives a round trip through binary16.
    #[must_use]
    pub const fn half_is_exact(&self) -> bool {
        self.half_exact
    }

    /// Smallest decimal count at which every value is an exact fixed-point
    /// number whose scaled integer stays below `ceiling`.
    #[must_use]
    pub fn exact_decimals(&self, ceiling: u64) -> Option<u8> {
        (0..=MAX_DECIMALS)
            .find(|slot| self.fixed_exact[*slot] && self.fixed_maximum[*slot] <= ceiling)
            .and_then(|slot| u8::try_from(slot).ok())
    }

    /// Smallest decimal count at which every value is exact, however large.
    #[must_use]
    pub fn exact_decimals_unbounded(&self) -> Option<u8> {
        self.exact_decimals(u64::from(u32::MAX))
    }
}

/// Ten to the power of a decimal count.
pub(crate) fn denominator(decimals: u8) -> u32 {
    POWERS[usize::from(decimals).min(MAX_DECIMALS)]
}

#[cfg(test)]
#[path = "precision_tests.rs"]
mod tests;
