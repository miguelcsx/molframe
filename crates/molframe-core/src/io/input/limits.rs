//! Hostile-input ceilings shared by readers and decoders.

use crate::diagnostic::{Code, Diagnostic};

/// 16 GiB: the largest published mmCIF files are under 2 GiB uncompressed, so
/// the ceiling sits an order of magnitude above any real entry while still
/// refusing an input that could exhaust a workstation.
const DEFAULT_DECOMPRESSED_BYTES: u64 = 16 << 30;
/// 2^28 rows: far beyond the largest deposited structure (tens of millions of
/// atoms) and consistent with the byte ceiling above.
const DEFAULT_ROWS_PER_CATEGORY: u64 = 1 << 28;
/// 2^24 distinct identifiers: real structures intern thousands.
const DEFAULT_DICTIONARY_ENTRIES: u32 = 1 << 24;
const DEFAULT_COMPRESSION_RATIO: u64 = 1_000;
const DEFAULT_NESTING_DEPTH: u32 = 64;

/// Ceilings that turn hostile input into diagnostics rather than exhaustion.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Limits {
    /// Largest accepted size after decompression, in bytes.
    ///
    /// Defaults to 16 GiB. The ratio ceiling ([`Self::compression_ratio`])
    /// bounds expansion relative to the bytes supplied; this absolute ceiling
    /// additionally stops a hostile or corrupt input from exhausting memory
    /// before any parse. A caller reading a genuinely larger file raises it, or
    /// uses [`Limits::unbounded`] to opt out of every absolute ceiling.
    pub decompressed_bytes: u64,
    /// Largest accepted ratio of decompressed to compressed size.
    pub compression_ratio: u64,
    /// Largest accepted declared row count for one category.
    ///
    /// Defaults to 2^28 rows; [`Limits::unbounded`] removes the ceiling.
    pub rows_per_category: u64,
    /// Deepest accepted nesting.
    pub nesting_depth: u32,
    /// Largest accepted identifier dictionary.
    ///
    /// Defaults to 2^24 entries; [`Limits::unbounded`] removes the ceiling.
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
    /// Limits with no absolute ceiling, for trusted input of any size.
    ///
    /// The compression-ratio and nesting-depth ceilings still apply: they
    /// protect against malformed input rather than large input.
    #[must_use]
    pub const fn unbounded() -> Self {
        Self {
            decompressed_bytes: u64::MAX,
            compression_ratio: DEFAULT_COMPRESSION_RATIO,
            rows_per_category: u64::MAX,
            nesting_depth: DEFAULT_NESTING_DEPTH,
            dictionary_entries: u32::MAX - 1,
        }
    }

    /// Builds the registered finding for a violated ceiling.
    #[must_use]
    pub fn exceeded(limit: &'static str, value: impl std::fmt::Display) -> Diagnostic {
        Diagnostic::new(Code::E1901)
            .with_context("limit", limit)
            .with_context("value", value.to_string())
    }
}
