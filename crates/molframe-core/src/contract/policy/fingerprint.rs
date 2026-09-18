//! Telling one policy or input from another.
//!
//! The standard hasher is seeded afresh each process, which is right for hash
//! tables and wrong here: a fingerprint that changes between runs identifies the
//! run rather than the thing. This one does not.

use std::fmt;
use std::hash::Hasher;

/// FNV-1a's specified 64-bit initial state.
const FNV1A64_OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
/// FNV-1a's specified 64-bit prime multiplier.
const FNV1A64_PRIME: u64 = 0x0000_0100_0000_01b3;

/// A stable content fingerprint.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct Fingerprint(u64);

impl Fingerprint {
    /// Wraps a finished hash.
    pub(super) const fn from_hash(hash: u64) -> Self {
        Self(hash)
    }

    /// Fingerprints a block of bytes.
    #[must_use]
    pub fn of(bytes: &[u8]) -> Self {
        let mut hasher = Fnv1a::new();
        hasher.write(bytes);
        Self(hasher.finish())
    }

    /// The raw value.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

impl fmt::Display for Fingerprint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "fnv1a64:{:016x}", self.0)
    }
}

impl fmt::Debug for Fingerprint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self}")
    }
}

/// A hasher whose output does not depend on how the program was started.
///
/// The standard hasher is randomly seeded per process, which is right for
/// hash tables and wrong here: a fingerprint that changes between runs cannot
/// identify anything.
pub(super) struct Fnv1a(u64);

impl Fnv1a {
    pub(super) const fn new() -> Self {
        Self(FNV1A64_OFFSET_BASIS)
    }
}

impl Hasher for Fnv1a {
    fn finish(&self) -> u64 {
        self.0
    }

    fn write(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.0 ^= u64::from(*byte);
            // FNV-1a is defined modulo 2^64, so wrapping is part of the
            // algorithm rather than an overflow fallback.
            self.0 = self.0.wrapping_mul(FNV1A64_PRIME);
        }
    }
}
