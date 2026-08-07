//! Hostile-input ceilings shared by readers and decoders.

use crate::diagnostic::{Code, Diagnostic};

const DEFAULT_DECOMPRESSED_BYTES: u64 = 4_u64 << 30;
const DEFAULT_COMPRESSION_RATIO: u64 = 1_000;
const DEFAULT_ROWS_PER_CATEGORY: u64 = 100_000_000;
const DEFAULT_NESTING_DEPTH: u32 = 64;
const DEFAULT_DICTIONARY_ENTRIES: u32 = 1_000_000;

/// Ceilings that turn hostile input into diagnostics rather than exhaustion.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Limits {
    /// Largest accepted size after decompression, in bytes.
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
            rows_per_category: DEFAULT_ROWS_PER_CATEGORY,
            nesting_depth: DEFAULT_NESTING_DEPTH,
            dictionary_entries: DEFAULT_DICTIONARY_ENTRIES,
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
