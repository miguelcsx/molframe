//! Types shared by sequence format readers and writers.

use crate::{A2mError, FastaRecord, FastqError, FastqRecord, NewickError, Tree};
use std::fmt::{Display, Formatter};

/// A supported sequence/alignment/tree serialization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequenceFormat {
    /// FASTA records.
    Fasta,
    /// FASTQ reads and qualities.
    Fastq,
    /// Stockholm multiple alignment.
    Stockholm,
    /// Clustal multiple alignment.
    Clustal,
    /// Relaxed sequential PHYLIP alignment.
    Phylip,
    /// A2M match/insertion alignment.
    A2m,
    /// A3M match/insertion alignment.
    A3m,
    /// Newick tree.
    Newick,
}

/// Parsed content without lossy coercion between record kinds.
#[derive(Clone, Debug, PartialEq)]
pub enum SequenceDocument {
    /// FASTA-shaped sequence or alignment records.
    Records(Vec<FastaRecord>),
    /// FASTQ records with encoded quality bytes.
    Fastq(Vec<FastqRecord>),
    /// One rooted tree.
    Tree(Tree),
}

/// Why common format dispatch refused input or output.
#[derive(Clone, Debug, PartialEq)]
pub enum SequenceFormatError {
    /// FASTQ-specific validation failed.
    Fastq(FastqError),
    /// A2M/A3M validation failed.
    A2m(A2mError),
    /// Newick parsing failed.
    Newick(NewickError),
    /// Format framing or declared counts are invalid.
    Invalid {
        /// Selected format.
        format: SequenceFormat,
        /// Stable explanation.
        reason: &'static str,
    },
    /// The document variant is incompatible with the requested writer.
    WrongDocument {
        /// Selected format.
        format: SequenceFormat,
    },
}

impl Display for SequenceFormatError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Fastq(error) => Display::fmt(error, formatter),
            Self::A2m(error) => Display::fmt(error, formatter),
            Self::Newick(error) => write!(formatter, "invalid Newick: {error:?}"),
            Self::Invalid { format, reason } => write!(formatter, "invalid {format:?}: {reason}"),
            Self::WrongDocument { format } => {
                write!(formatter, "document is incompatible with {format:?}")
            }
        }
    }
}

impl std::error::Error for SequenceFormatError {}

impl From<FastqError> for SequenceFormatError {
    fn from(value: FastqError) -> Self {
        Self::Fastq(value)
    }
}

impl From<A2mError> for SequenceFormatError {
    fn from(value: A2mError) -> Self {
        Self::A2m(value)
    }
}

impl From<NewickError> for SequenceFormatError {
    fn from(value: NewickError) -> Self {
        Self::Newick(value)
    }
}
