//! A fixed multiply-shift hasher for keys that are already integers.
//!
//! `SipHash`, the standard-library default, is built to resist adversarial key
//! choice. Nothing here hashes attacker-chosen keys: the keys are atom indices,
//! symbol identifiers and index pairs the engine assigned itself. Paying a
//! keyed pseudo-random function for them costs roughly a nanosecond per lookup,
//! which at a billion contacts is a second of pure overhead.
//!
//! The multiplier and shifts are compile-time constants rather than process
//! seeded, so iteration order over a map built from the same keys is identical
//! across runs and machines — required by NFR-307.

use std::hash::{BuildHasherDefault, Hasher};

/// The golden-ratio odd multiplier, as used by Fibonacci hashing.
///
/// Multiplying by an odd constant is a bijection on `u64`, so no two distinct
/// absorbed sequences collide through the multiply itself.
const MULTIPLIER: u64 = 0x9E37_79B9_7F4A_7C15;

/// The final xor-shift distance, chosen to fold high entropy into low bits.
///
/// The multiply concentrates entropy in the high bits, while a hash table
/// selects its bucket from the low ones; without this fold, keys differing only
/// in their low bits would land in the same bucket.
const FOLD: u32 = 29;

/// A `Hasher` for integer keys: one multiply and one xor-shift per word.
#[derive(Debug, Clone, Copy, Default)]
pub struct IdentityHasher {
    /// The running hash, updated once per absorbed word.
    state: u64,
}

impl IdentityHasher {
    /// Folds one word into the state.
    ///
    /// Order dependent, so a two-field key such as `(u32, u32)` does not
    /// collide with its own reverse.
    fn absorb(&mut self, value: u64) {
        self.state = (self.state ^ value).wrapping_mul(MULTIPLIER);
        self.state ^= self.state >> FOLD;
    }
}

impl Hasher for IdentityHasher {
    fn finish(&self) -> u64 {
        self.state
    }

    /// Absorbs arbitrary bytes eight at a time, then the remainder.
    ///
    /// Present so the type satisfies `Hasher` for any key; the integer paths
    /// below are what the engine actually uses.
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for chunk in &mut chunks {
            let mut word = [0u8; 8];
            word.copy_from_slice(chunk);
            self.absorb(u64::from_le_bytes(word));
        }

        let remainder = chunks.remainder();
        if !remainder.is_empty() {
            let mut word = [0u8; 8];
            match word.get_mut(..remainder.len()) {
                Some(head) => head.copy_from_slice(remainder),
                None => return,
            }
            self.absorb(u64::from_le_bytes(word));
        }
    }

    fn write_u8(&mut self, value: u8) {
        self.absorb(u64::from(value));
    }

    fn write_u16(&mut self, value: u16) {
        self.absorb(u64::from(value));
    }

    fn write_u32(&mut self, value: u32) {
        self.absorb(u64::from(value));
    }

    fn write_u64(&mut self, value: u64) {
        self.absorb(value);
    }

    fn write_usize(&mut self, value: usize) {
        self.absorb(value as u64);
    }
}

/// The `BuildHasher` to give `HashMap`/`HashSet` over integer keys.
pub type IdentityBuildHasher = BuildHasherDefault<IdentityHasher>;

/// A `HashMap` keyed by integers, hashed by [`IdentityHasher`].
pub type IdentityHashMap<K, V> = std::collections::HashMap<K, V, IdentityBuildHasher>;

/// A `HashSet` of integers, hashed by [`IdentityHasher`].
pub type IdentityHashSet<K> = std::collections::HashSet<K, IdentityBuildHasher>;

#[cfg(test)]
#[path = "mix_tests.rs"]
mod tests;
