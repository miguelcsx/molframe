//! Stable diagnostic codes for every error the sequence kernels return.
//!
//! Alignment and numeric-limit failures name the code of their kind; every
//! failure about the content of sequences, alphabets, trees, scoring profiles
//! or file formats shares the data-validity code, and the message carries which.

use crate::{
    A2mError, AlignError, AlphabetError, FastqError, KmerTableError, MatrixError, MsaError,
    NewickError, RegionAlignError, RegionError, RerootError, SequenceError, SequenceFormatError,
    SimilarKmerError,
};
use molframe_core::{Code, diagnostic_from};

diagnostic_from!(AlphabetError, |_error| Code::E2101);
diagnostic_from!(SequenceError, |_error| Code::E2101);
diagnostic_from!(SequenceFormatError, |_error| Code::E2101);
diagnostic_from!(A2mError, |_error| Code::E2101);
diagnostic_from!(FastqError, |_error| Code::E2101);
diagnostic_from!(NewickError, |_error| Code::E2101);
diagnostic_from!(RerootError, |_error| Code::E2101);
diagnostic_from!(MatrixError, |_error| Code::E2101);
diagnostic_from!(KmerTableError, |_error| Code::E5101);
diagnostic_from!(RegionError, |_error| Code::E5101);
diagnostic_from!(AlignError, |error| match error {
    AlignError::NumericOverflow => Code::E1903,
    AlignError::NoAlignmentPath => Code::E5105,
});
diagnostic_from!(RegionAlignError, |error| match error {
    RegionAlignError::Region(_) => Code::E5101,
    RegionAlignError::Alignment(AlignError::NumericOverflow) => Code::E1903,
    RegionAlignError::Alignment(AlignError::NoAlignmentPath) => Code::E5105,
});
diagnostic_from!(MsaError, |error| match error {
    MsaError::InvalidGapScore | MsaError::InvalidMemoryLimit { .. } => Code::E5101,
    MsaError::DimensionOverflow | MsaError::NumericOverflow => Code::E1903,
    MsaError::MemoryLimit { .. } => Code::E7001,
});
diagnostic_from!(SimilarKmerError, |error| match error {
    SimilarKmerError::EmptyInput | SimilarKmerError::DuplicateSymbol => Code::E5101,
    SimilarKmerError::CandidateLimit { .. } => Code::E1901,
    SimilarKmerError::NumericOverflow => Code::E1903,
});

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
