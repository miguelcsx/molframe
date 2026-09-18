//! Hostile-input ceilings shared by readers and decoders.

use crate::diagnostic::{Code, Diagnostic};

const DEFAULT_DECOMPRESSED_BYTES: u64 = u64::MAX;
const DEFAULT_COMPRESSION_RATIO: u64 = 1_000;
const DEFAULT_NESTING_DEPTH: u32 = 64;

/// Ceilings that turn hostile input into diagnostics rather than exhaustion.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Limits {
    /// Largest accepted size after decompression, in bytes.
    ///
    /// Unbounded by default. The ceiling exists to stop a small compressed
    /// input expanding without limit, and that job is done by
    /// [`Self::compression_ratio`], which bounds expansion against the bytes
    /// actually supplied. An absolute ceiling additionally refuses large
    /// *uncompressed* input, where there is nothing to defend against: the
    /// file is already that size, and reading it cannot multiply it.
    ///
    /// A caller who wants an absolute ceiling sets one.
    pub decompressed_bytes: u64,
    /// Largest accepted ratio of decompressed to compressed size.
    pub compression_ratio: u64,
    /// Largest accepted declared row count for one category.
    pub rows_per_category: u64,
    /// Deepest accepted nesting.
    pub nesting_depth: u32,
    /// Largest accepted identifier dictionary.
    pub dictionary_entries: u32,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            decompressed_bytes: DEFAULT_DECOMPRESSED_BYTES,
            compression_ratio: DEFAULT_COMPRESSION_RATIO,
            rows_per_category: u64::MAX,
            nesting_depth: DEFAULT_NESTING_DEPTH,
            dictionary_entries: u32::MAX - 1,
        }
    }
}

impl Limits {
    /// Builds the registered finding for a violated ceiling.
    #[must_use]
    pub fn exceeded(limit: &'static str, value: impl std::fmt::Display) -> Diagnostic {
        Diagnostic::new(Code::E1901)
            .with_context("limit", limit)
            .with_context("value", value.to_string())
    }
}
