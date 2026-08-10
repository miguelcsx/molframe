//! Reading and writing the Stockholm alignment format.
//!
//! Stockholm records a multiple alignment as `name sequence` lines, optionally
//! split into interleaved blocks that each carry a slice of every sequence, with
//! markup lines beginning with `#` and a `//` terminator. Parsing concatenates
//! the blocks for each name in first-seen order and drops the markup, keeping the
//! gapped sequences; writing emits a single block, so a parse of a written
//! alignment returns the same records.

use crate::fasta::FastaRecord;

/// Parses Stockholm text into alignment records in first-seen order.
///
/// Markup lines (`#…`) are ignored and parsing stops at the `//` terminator.
/// Whitespace within a sequence is dropped, and a name appearing in several
/// interleaved blocks has its pieces concatenated.
#[must_use]
pub fn parse_stockholm(text: &str) -> Vec<FastaRecord> {
    let mut records: Vec<FastaRecord> = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with("//") {
            break;
        }
        let Some((name, rest)) = line.split_once(char::is_whitespace) else {
            continue;
        };
        let residues = rest.bytes().filter(|byte| !byte.is_ascii_whitespace());
        match records.iter_mut().find(|record| record.id == name) {
            Some(record) => record.sequence.extend(residues),
            None => records.push(FastaRecord {
                id: name.to_string(),
                description: String::new(),
                sequence: residues.collect(),
            }),
        }
    }
    records
}

/// Writes alignment records as a single-block Stockholm file.
#[must_use]
pub fn write_stockholm(records: &[FastaRecord]) -> String {
    let mut out = String::from("# STOCKHOLM 1.0\n");
    for record in records {
        out.push_str(&record.id);
        out.push(' ');
        out.push_str(&String::from_utf8_lossy(&record.sequence));
        out.push('\n');
    }
    out.push_str("//\n");
    out
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
