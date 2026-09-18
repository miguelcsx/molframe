//! Reading and writing FASTA sequence files.
//!
//! FASTA is the plainest sequence format: a header line beginning with `>`
//! naming the sequence, then its residues over one or more lines. Parsing keeps
//! the records in file order and strips the line breaks from each sequence;
//! writing wraps long sequences to a fixed width, so a parse followed by a write
//! and another parse returns the same records.

/// One FASTA record: its identifier, free-text description, and residues.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastaRecord {
    /// The identifier — the first whitespace-delimited token of the header.
    pub id: String,
    /// The rest of the header line, if any.
    pub description: String,
    /// The residue bytes, with all whitespace removed.
    pub sequence: Vec<u8>,
}

/// The line width used when writing sequences.
const WRAP: usize = 60;

/// Parses FASTA text into records, preserving their order.
///
/// Anything before the first header is ignored, and whitespace within a
/// sequence is dropped. Malformed input does not fail; it simply yields whatever
/// records it can.
#[must_use]
pub fn parse_fasta(text: &str) -> Vec<FastaRecord> {
    let mut records = Vec::new();
    let mut current: Option<FastaRecord> = None;
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(record) = current.take() {
                records.push(record);
            }
            current = Some(new_record(header.trim()));
        } else if let Some(record) = current.as_mut() {
            for &byte in line.as_bytes() {
                if !byte.is_ascii_whitespace() {
                    record.sequence.push(byte);
                }
            }
        }
    }
    if let Some(record) = current.take() {
        records.push(record);
    }
    records
}

/// Writes records as FASTA text with sequences wrapped to a fixed width.
#[must_use]
pub fn write_fasta(records: &[FastaRecord]) -> String {
    let mut out = String::new();
    for record in records {
        out.push('>');
        out.push_str(&record.id);
        if !record.description.is_empty() {
            out.push(' ');
            out.push_str(&record.description);
        }
        out.push('\n');
        for chunk in record.sequence.chunks(WRAP) {
            out.push_str(&String::from_utf8_lossy(chunk));
            out.push('\n');
        }
    }
    out
}

/// Splits a header into its identifier and the remaining description.
fn new_record(header: &str) -> FastaRecord {
    let (id, description) = match header.split_once(char::is_whitespace) {
        Some((id, rest)) => (id.to_string(), rest.trim().to_string()),
        None => (header.to_string(), String::new()),
    };
    FastaRecord {
        id,
        description,
        sequence: Vec::new(),
    }
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
