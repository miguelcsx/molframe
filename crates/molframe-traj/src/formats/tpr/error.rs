//! TPR decoding failures.

use thiserror::Error;

/// Why a GROMACS TPR topology could not be decoded.
#[derive(Clone, Debug, Error, PartialEq, Eq)]
pub enum TprError {
    /// The file ended before a declared value was complete.
    #[error("TPR input is truncated at byte {offset}")]
    Truncated {
        /// Byte position at which data was required.
        offset: usize,
    },
    /// A length or count is negative or exceeds the input.
    #[error("invalid TPR {field} value {value} at byte {offset}")]
    InvalidCount {
        /// Count vocabulary.
        field: &'static str,
        /// Invalid signed value.
        value: i64,
        /// Byte position.
        offset: usize,
    },
    /// Header marker does not identify a TPR file.
    #[error("input does not begin with a GROMACS VERSION header")]
    InvalidMagic,
    /// Precision is neither single nor double.
    #[error("unsupported TPR real width {0} bytes")]
    UnsupportedPrecision(i32),
    /// The format version is outside the versions this reader validates.
    #[error("unsupported TPR format version {0}")]
    UnsupportedVersion(i32),
    /// A string is not valid UTF-8.
    #[error("TPR {field} is not valid UTF-8")]
    InvalidText {
        /// String vocabulary.
        field: &'static str,
    },
    /// An integer index is outside its table.
    #[error("TPR {field} index {index} is outside {length} entries")]
    InvalidIndex {
        /// Index vocabulary.
        field: &'static str,
        /// Invalid position.
        index: usize,
        /// Table length.
        length: usize,
    },
    /// The file carries no topology section.
    #[error("TPR file has no topology section")]
    MissingTopology,
    /// A force-field function introduced after this version was encountered.
    #[error("unsupported TPR force-field function {0}")]
    UnsupportedFunction(i32),
    /// The 2020 beta serializer cannot be distinguished safely.
    #[error("GROMACS 2020 beta TPR serializer is not supported")]
    UnsupportedBetaSerializer,
    /// Arithmetic on declared sizes overflowed.
    #[error("TPR declared sizes overflow the platform index range")]
    SizeOverflow,
}
