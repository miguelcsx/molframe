//! The deterministic scalar sequence every synthetic fixture is drawn from.
//!
//! Scale fixtures must be reproducible from a seed alone, so that a failure at
//! ten thousand copies can be reproduced without storing ten thousand copies.
//! That rules out a seeded random-number generator whose algorithm is an
//! implementation detail of another crate, because the sequence would then
//! change with a dependency update and the fixture would stop being the fixture.
//!
//! `SplitMix64` is used instead: four lines, fixed constants, and the same
//! sequence on every platform and every release. It is not a cryptographic
//! generator and nothing here needs one — the requirement is reproducibility,
//! not unpredictability.

use crate::numeric::u64_to_f64;

/// The golden-ratio increment `SplitMix64` advances its state by.
const INCREMENT: u64 = 0x9e37_79b9_7f4a_7c15;
/// First avalanche multiplier.
const FIRST_MULTIPLIER: u64 = 0xbf58_476d_1ce4_e5b9;
/// Second avalanche multiplier.
const SECOND_MULTIPLIER: u64 = 0x94d0_49bb_1331_11eb;

/// A reproducible scalar sequence.
///
/// Two `Seed`s created from the same value produce the same sequence, on any
/// platform and in any release.
///
/// # Examples
///
/// ```
/// use molframe_bench::Seed;
///
/// let mut left = Seed::new(7);
/// let mut right = Seed::new(7);
/// assert_eq!(left.next_unit(), right.next_unit());
/// ```
#[derive(Clone, Debug)]
pub struct Seed {
    state: u64,
}

impl Seed {
    /// Starts a sequence from `seed`.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Advances the sequence and returns the next 64-bit value.
    pub fn next_bits(&mut self) -> u64 {
        self.state = self.state.wrapping_add(INCREMENT);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(FIRST_MULTIPLIER);
        value = (value ^ (value >> 27)).wrapping_mul(SECOND_MULTIPLIER);
        value ^ (value >> 31)
    }

    /// Advances the sequence and returns a value in `[0, 1)`.
    ///
    /// The top fifty-three bits are used, which is exactly what an `f64`
    /// mantissa holds, so the division is exact and no bit is wasted.
    pub fn next_unit(&mut self) -> f64 {
        const MANTISSA_BITS: u32 = 53;
        let bits = self.next_bits() >> (u64::BITS - MANTISSA_BITS);
        u64_to_f64(bits) / u64_to_f64(1u64 << MANTISSA_BITS)
    }

    /// Advances the sequence and returns a value in `[-half_width, half_width)`.
    pub fn next_signed(&mut self, half_width: f64) -> f64 {
        (self.next_unit() - 0.5) * 2.0 * half_width
    }

    /// A sequence derived from this one, independent of how far this one runs.
    ///
    /// Copies of a tile draw from derived sequences rather than from a shared
    /// one, so copy `n` is reproducible without generating copies `0..n`.
    #[must_use]
    pub fn derive(&self, ordinal: u64) -> Self {
        let mut derived = Self::new(self.state ^ ordinal.wrapping_mul(INCREMENT));
        let _warm = derived.next_bits();
        derived
    }
}

#[cfg(test)]
#[path = "seed_tests.rs"]
mod tests;
