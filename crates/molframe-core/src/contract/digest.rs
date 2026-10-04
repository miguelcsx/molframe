//! A cryptographic digest of the bytes an analysis read.
//!
//! [`Fingerprint`](super::Fingerprint) is a fast 64-bit hash, enough to tell one
//! policy or cache entry from another and nothing more: it is not collision
//! resistant, so it cannot back the claim that these are exactly the bytes an
//! earlier result was computed from. That claim is what re-execution makes, so it
//! rests on this, SHA-256, and the two types are kept apart so that one is never
//! mistaken for the other.

use sha2::{Digest, Sha256};
use std::fmt;

/// The SHA-256 digest of a block of bytes.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct ContentDigest([u8; 32]);

impl ContentDigest {
    /// Digests a block of bytes.
    #[must_use]
    pub fn of(bytes: &[u8]) -> Self {
        Self(Sha256::digest(bytes).into())
    }

    /// The 32 digest bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    /// Reads a digest from its 64 lowercase hexadecimal digits.
    #[must_use]
    pub fn from_hex(text: &str) -> Option<Self> {
        let text = match text.strip_prefix("sha256:") {
            Some(digits) => digits,
            None => text,
        };
        if text.len() != 64 || !text.bytes().all(|byte| byte.is_ascii_hexdigit()) {
            return None;
        }
        let mut bytes = [0_u8; 32];
        for (slot, pair) in bytes.iter_mut().zip(text.as_bytes().chunks(2)) {
            let digits = std::str::from_utf8(pair).ok()?;
            *slot = u8::from_str_radix(digits, 16).ok()?;
        }
        Some(Self(bytes))
    }
}

impl fmt::Display for ContentDigest {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("sha256:")?;
        for byte in self.0 {
            write!(formatter, "{byte:02x}")?;
        }
        Ok(())
    }
}

impl fmt::Debug for ContentDigest {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self}")
    }
}

#[cfg(test)]
#[path = "digest_tests.rs"]
mod tests;
