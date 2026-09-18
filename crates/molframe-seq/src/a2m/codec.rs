//! A2M alignment validation and match-column projection.

use crate::FastaRecord;
use crate::fasta::{parse_fasta, write_fasta};

/// Why FASTA-shaped A2M text is not a coherent alignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum A2mError {
    /// A byte has no A2M column meaning.
    InvalidSymbol {
        /// Zero-based record position.
        record: usize,
        /// Zero-based byte position.
        offset: usize,
    },
    /// Records disagree on their number of match-state columns.
    MatchLength {
        /// Zero-based record position.
        record: usize,
        /// Match-column count established by the first record.
        expected: usize,
        /// Match-column count in this record.
        actual: usize,
    },
}

impl std::fmt::Display for A2mError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidSymbol { record, offset } => {
                write!(
                    formatter,
                    "A2M record {record} has invalid symbol at {offset}"
                )
            }
            Self::MatchLength {
                record,
                expected,
                actual,
            } => write!(
                formatter,
                "A2M record {record} has {actual} match columns, expected {expected}"
            ),
        }
    }
}

impl std::error::Error for A2mError {}

/// Parses A2M and verifies a shared match-state coordinate system.
///
/// Uppercase letters and `-` occupy match columns; lowercase letters and `.`
/// occupy insertion columns. Case and insertion padding remain unchanged.
///
/// # Errors
///
/// Returns [`A2mError`] for invalid bytes or inconsistent match-column counts.
pub fn parse_a2m(text: &str) -> Result<Vec<FastaRecord>, A2mError> {
    let records = parse_fasta(text);
    let mut expected = None;
    for (record, entry) in records.iter().enumerate() {
        validate(&entry.sequence, record)?;
        let actual = a2m_match_columns(&entry.sequence).len();
        match expected {
            Some(expected) if expected != actual => {
                return Err(A2mError::MatchLength {
                    record,
                    expected,
                    actual,
                });
            }
            None => expected = Some(actual),
            _ => {}
        }
    }
    Ok(records)
}

/// Projects an A2M row onto uppercase match states and deletions.
#[must_use]
pub fn a2m_match_columns(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .copied()
        .filter(|symbol| symbol.is_ascii_uppercase() || *symbol == b'-')
        .collect()
}

/// Writes A2M after validating symbols and match-state coordinates.
///
/// # Errors
///
/// Returns [`A2mError`] under the same conditions as [`parse_a2m`].
pub fn write_a2m(records: &[FastaRecord]) -> Result<String, A2mError> {
    let text = write_fasta(records);
    parse_a2m(&text)?;
    Ok(text)
}

fn validate(sequence: &[u8], record: usize) -> Result<(), A2mError> {
    match sequence
        .iter()
        .position(|symbol| !symbol.is_ascii_alphabetic() && !matches!(symbol, b'-' | b'.'))
    {
        Some(offset) => Err(A2mError::InvalidSymbol { record, offset }),
        None => Ok(()),
    }
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
