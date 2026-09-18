//! Exact and deterministic bucketed indexes over sequence sets.

use super::SeedPattern;

/// Standard 32-bit FNV-1a offset basis, widened to the platform index type.
const FNV1A_OFFSET_BASIS_32: usize = 0x811c_9dc5;
/// Standard 32-bit FNV-1a prime; multiplication wraps by algorithm definition.
const FNV1A_PRIME_32: usize = 0x0100_0193;

/// Storage layout for an immutable k-mer index.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KmerStorage {
    /// One sorted map over exact seed bytes.
    Exact,
    /// A fixed number of hash buckets with exact collision resolution.
    Bucketed {
        /// Number of buckets allocated.
        buckets: usize,
    },
}

/// Explicit k-mer extraction and storage policy.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KmerTableOptions {
    /// Contiguous word length when no spaced pattern is supplied.
    pub k: usize,
    /// Optional spaced-seed mask; its weight becomes the indexed word length.
    pub pattern: Option<SeedPattern>,
    /// Index storage layout.
    pub storage: KmerStorage,
}

/// One occurrence of an exact indexed word.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct KmerHit {
    /// Input sequence position.
    pub sequence: usize,
    /// Zero-based window start.
    pub position: usize,
}

/// Invalid extraction or bucket configuration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum KmerTableError {
    /// Contiguous word length is zero.
    EmptyWord,
    /// Bucketed storage requested zero buckets.
    EmptyBuckets,
}

impl std::fmt::Display for KmerTableError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::EmptyWord => formatter.write_str("k-mer length must be positive"),
            Self::EmptyBuckets => formatter.write_str("bucketed k-mer table needs buckets"),
        }
    }
}

impl std::error::Error for KmerTableError {}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Entry {
    word: Vec<u8>,
    hits: Vec<KmerHit>,
}

/// Immutable exact-match index with deterministic iteration and collision handling.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct KmerTable {
    word_length: usize,
    buckets: Vec<Vec<Entry>>,
}

impl KmerTable {
    /// Indexes every sequence under an explicit contiguous or spaced policy.
    ///
    /// # Errors
    ///
    /// Returns [`KmerTableError`] for a zero word length or zero bucket count.
    pub fn build(sequences: &[&[u8]], options: &KmerTableOptions) -> Result<Self, KmerTableError> {
        let word_length = match &options.pattern {
            Some(pattern) => pattern.weight(),
            None if options.k > 0 => options.k,
            None => return Err(KmerTableError::EmptyWord),
        };
        let bucket_count = match options.storage {
            KmerStorage::Exact => 1,
            KmerStorage::Bucketed { buckets: 0 } => return Err(KmerTableError::EmptyBuckets),
            KmerStorage::Bucketed { buckets } => buckets,
        };
        let mut table = Self {
            word_length,
            buckets: vec![Vec::new(); bucket_count],
        };
        for (sequence_index, sequence) in sequences.iter().enumerate() {
            match &options.pattern {
                Some(pattern) => {
                    for (position, word) in pattern.seeds(sequence) {
                        table.insert(
                            word,
                            KmerHit {
                                sequence: sequence_index,
                                position,
                            },
                        );
                    }
                }
                None => {
                    for (position, word) in sequence.windows(options.k).enumerate() {
                        table.insert(
                            word.to_vec(),
                            KmerHit {
                                sequence: sequence_index,
                                position,
                            },
                        );
                    }
                }
            }
        }
        for bucket in &mut table.buckets {
            bucket.sort_by(|left, right| left.word.cmp(&right.word));
        }
        Ok(table)
    }

    /// Exact occurrences of `word`, ordered by sequence then position.
    #[must_use]
    pub fn query(&self, word: &[u8]) -> &[KmerHit] {
        if word.len() != self.word_length {
            return &[];
        }
        let bucket = self.bucket(word);
        match self.buckets[bucket].binary_search_by(|entry| entry.word.as_slice().cmp(word)) {
            Ok(index) => &self.buckets[bucket][index].hits,
            Err(_) => &[],
        }
    }

    /// Indexed word length after applying a spaced pattern.
    #[must_use]
    pub const fn word_length(&self) -> usize {
        self.word_length
    }

    fn insert(&mut self, word: Vec<u8>, hit: KmerHit) {
        let bucket = self.bucket(&word);
        if let Some(entry) = self.buckets[bucket]
            .iter_mut()
            .find(|entry| entry.word == word)
        {
            entry.hits.push(hit);
        } else {
            self.buckets[bucket].push(Entry {
                word,
                hits: vec![hit],
            });
        }
    }

    fn bucket(&self, word: &[u8]) -> usize {
        let hash = word.iter().fold(FNV1A_OFFSET_BASIS_32, |hash, byte| {
            (hash ^ usize::from(*byte)).wrapping_mul(FNV1A_PRIME_32)
        });
        hash % self.buckets.len()
    }
}

#[cfg(test)]
#[path = "table_tests.rs"]
mod tests;
