//! Exact, spaced, bucketed, and similarity-aware sequence words.

mod counts;
mod similar;
mod spaced;
mod table;

pub use counts::{kmer_counts, minimizers};
pub use similar::{SimilarKmer, SimilarKmerError, SimilarKmerOptions, similar_kmers};
pub use spaced::{SeedPattern, SeedPatternError};
pub use table::{KmerHit, KmerStorage, KmerTable, KmerTableError, KmerTableOptions};
