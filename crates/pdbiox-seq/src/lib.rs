//! Pairwise sequence alignment.
//!
//! Algorithms accept compact byte slices, while [`Sequence`] and [`Alphabet`]
//! provide validated typed storage for DNA, RNA, proteins, nucleotides, and
//! explicit custom alphabets. Alignments use affine gap costs, where opening a
//! gap and extending it are charged separately, because a single long insertion
//! is biologically cheaper than many scattered ones and a linear cost cannot
//! express that.
//!
//! Every alignment is deterministic: ties in the traceback are broken by a
//! fixed preference for a match over a gap, and a gap in the second sequence
//! over a gap in the first, so the same inputs always return the same path.

#![forbid(unsafe_code)]
// The alignment matrices are conventionally named `ix` and `iy`; keeping those
// names matches the algorithm they implement more than a lint's spelling rule.

pub mod a2m;
pub mod a3m;
pub mod align;
pub mod alphabet;
pub mod clustal;
pub mod fasta;
pub mod fastq;
pub mod format;
pub mod kmer;
pub mod matrix;
pub mod msa;
pub mod nj;
pub mod phylip;
pub mod scoring;
pub mod seed;
pub mod stockholm;
pub mod tree;
pub mod upgma;

pub use a2m::{A2mError, a2m_match_columns, parse_a2m, write_a2m};
pub use a3m::{a3m_match_columns, parse_a3m, write_a3m};
pub use align::{
    AlignError, Alignment, AlignmentMode, Column, RegionAlignError, RegionError, RegionOptions,
    align_region, global, global_banded, global_matrix, global_score, local, local_matrix,
    semi_global, semi_global_matrix,
};
pub use alphabet::{
    Alphabet, AlphabetError, Code, CustomAlphabet, DNA, DnaAlphabet, NUCLEOTIDE,
    NucleotideAlphabet, PROTEIN, ProteinAlphabet, RNA, RnaAlphabet, Sequence, SequenceError,
};
pub use clustal::{parse_clustal, write_clustal};
pub use fasta::{FastaRecord, parse_fasta, write_fasta};
pub use fastq::{FastqError, FastqRecord, parse_fastq, write_fastq};
pub use format::{
    SequenceDocument, SequenceFormat, SequenceFormatError, SequenceReader, SequenceWriter,
    read_sequence, write_sequence,
};
pub use kmer::{
    KmerHit, KmerStorage, KmerTable, KmerTableError, KmerTableOptions, SeedPattern,
    SeedPatternError, SimilarKmer, SimilarKmerError, SimilarKmerOptions, kmer_counts, minimizers,
    similar_kmers,
};
pub use matrix::{
    MatrixError, MatrixIdentity, MatrixProfile, SubstitutionMatrix, blosum62, load_matrix,
};
pub use msa::{MsaError, MsaOptions, msa, progressive_msa};
pub use nj::neighbor_joining;
pub use phylip::{parse_phylip, write_phylip};
pub use scoring::{Score, Scoring};
pub use seed::seed_and_extend;
pub use stockholm::{parse_stockholm, write_stockholm};
pub use tree::{LadderDirection, NewickError, RerootError, TraversalOrder, Tree};
pub use upgma::upgma;
