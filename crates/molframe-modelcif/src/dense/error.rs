//! Failures of dense confidence storage, import and extraction.

use std::fmt;

/// Why a dense confidence matrix or score vector could not be built.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ConfidenceError {
    /// The matrix would exceed the memory budget in the chosen storage.
    Budget {
        /// Bytes the matrix needs.
        required: u64,
        /// Bytes the policy allows.
        limit: u64,
    },
    /// The allocator refused the matrix buffer.
    Allocation,
    /// A value is not finite.
    NonFinite,
    /// A cell index lies outside the matrix.
    OutOfRange {
        /// Row index.
        row: usize,
        /// Column index.
        column: usize,
        /// Matrix side.
        size: usize,
    },
    /// The matrix side disagrees with the residue count of the structure.
    SizeMismatch {
        /// Residues in the structure or score vector.
        residues: usize,
        /// Side of the matrix.
        matrix: usize,
    },
    /// The matrix has cells the source never stated.
    Incomplete {
        /// Cells a full matrix has.
        expected: u64,
        /// Cells the source supplied.
        found: u64,
    },
    /// A pair row names a residue the structure does not have.
    UnknownResidue {
        /// Label chain identifier.
        chain: Box<str>,
        /// Label sequence identifier.
        sequence: i64,
    },
    /// The same pair or residue is stated twice with different values.
    Conflict {
        /// Row of the repeated cell.
        row: usize,
        /// Column of the repeated cell.
        column: usize,
    },
    /// The JSON document is not a recognised predicted-aligned-error layout.
    Json {
        /// Byte offset of the first problem.
        offset: usize,
        /// What was expected there.
        expected: &'static str,
    },
    /// Rows of the matrix have different lengths.
    Ragged {
        /// Row whose length differs.
        row: usize,
        /// Length that row should have.
        expected: usize,
        /// Length it has.
        found: usize,
    },
}

impl fmt::Display for ConfidenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Budget { required, limit } => write!(
                formatter,
                "matrix needs {required} bytes, exceeding the {limit}-byte budget"
            ),
            Self::Allocation => formatter.write_str("matrix allocation was refused"),
            Self::NonFinite => formatter.write_str("confidence value is not finite"),
            Self::OutOfRange { row, column, size } => {
                write!(
                    formatter,
                    "cell ({row}, {column}) is outside a {size} x {size} matrix"
                )
            }
            Self::SizeMismatch { residues, matrix } => write!(
                formatter,
                "matrix side {matrix} does not match {residues} residues"
            ),
            Self::Incomplete { expected, found } => {
                write!(formatter, "matrix has {found} of {expected} cells")
            }
            Self::UnknownResidue { chain, sequence } => {
                write!(
                    formatter,
                    "residue {chain}:{sequence} is not in the structure"
                )
            }
            Self::Conflict { row, column } => {
                write!(
                    formatter,
                    "cell ({row}, {column}) is stated twice with different values"
                )
            }
            Self::Json { offset, expected } => {
                write!(formatter, "JSON offset {offset}: expected {expected}")
            }
            Self::Ragged {
                row,
                expected,
                found,
            } => write!(
                formatter,
                "row {row} has {found} entries where {expected} were expected"
            ),
        }
    }
}

impl std::error::Error for ConfidenceError {}
