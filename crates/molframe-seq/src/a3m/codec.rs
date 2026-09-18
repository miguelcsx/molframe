//! Reading the A3M alignment format and reducing it to match columns.
//!
//! A3M is FASTA-shaped, but its letter case carries meaning: upper-case residues
//! and `-` gaps are the columns that align to the query, while lower-case letters
//! are insertions relative to the query and `.` marks an insertion-column gap. So
//! two A3M sequences can differ in length yet share a match-column length once
//! the insertions are removed. Parsing reads the raw records; the match-column
//! reduction drops the insertions so the sequences line up.

use crate::a2m::{A2mError, parse_a2m};
use crate::fasta::{FastaRecord, write_fasta};

/// Parses A3M text into raw records, preserving case, gaps, and insertion dots.
///
/// A3M shares FASTA's header/sequence structure, so the records come back the
/// same way; apply [`a3m_match_columns`] to a record's sequence to get its
/// aligned columns.
///
/// # Errors
///
/// Returns [`A2mError`] for malformed FASTA structure, invalid A3M symbols, or
/// unequal match-column counts.
pub fn parse_a3m(text: &str) -> Result<Vec<FastaRecord>, A2mError> {
    parse_a2m(text)
}

/// Writes A3M after validating symbols and shared match-state coordinates.
///
/// # Errors
///
/// Returns [`A2mError`] for invalid symbols or unequal match-column counts.
pub fn write_a3m(records: &[FastaRecord]) -> Result<String, A2mError> {
    let text = write_fasta(records);
    parse_a2m(&text)?;
    Ok(text)
}

/// Reduces an A3M sequence to its match columns.
///
/// Keeps upper-case residues and `-` gaps — the columns that align to the query —
/// and drops lower-case insertions and their `.` placeholders. Every record of a
/// well-formed A3M alignment reduces to the same length.
#[must_use]
pub fn a3m_match_columns(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .copied()
        .filter(|byte| byte.is_ascii_uppercase() || *byte == b'-')
        .collect()
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
