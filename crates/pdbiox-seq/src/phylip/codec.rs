//! Reading and writing the PHYLIP sequence format.
//!
//! PHYLIP opens with a count and a length, then one line per sequence: a name
//! followed by its residues. This reads the relaxed sequential form, where the
//! name is the first whitespace-delimited token and the rest of the line is the
//! sequence, and writes it the same way. The header's numbers are recomputed on
//! write, so a parse of a written file returns the same records.

use crate::fasta::FastaRecord;

/// Parses relaxed sequential PHYLIP into records in file order.
///
/// The first non-blank line is the count/length header and is skipped; each
/// remaining line contributes one record whose name is its first token.
#[must_use]
pub fn parse_phylip(text: &str) -> Vec<FastaRecord> {
    let mut records = Vec::new();
    let mut seen_header = false;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if !seen_header {
            seen_header = true;
            continue;
        }
        let Some((name, rest)) = line.split_once(char::is_whitespace) else {
            continue;
        };
        records.push(FastaRecord {
            id: name.to_string(),
            description: String::new(),
            sequence: rest
                .bytes()
                .filter(|byte| !byte.is_ascii_whitespace())
                .collect(),
        });
    }
    records
}

/// Writes records as relaxed sequential PHYLIP with a recomputed header.
#[must_use]
pub fn write_phylip(records: &[FastaRecord]) -> String {
    let longest = match records.iter().map(|record| record.sequence.len()).max() {
        Some(length) => length,
        None => 0,
    };
    let mut out = format!("{} {}\n", records.len(), longest);
    for record in records {
        out.push_str(&record.id);
        out.push(' ');
        out.push_str(&String::from_utf8_lossy(&record.sequence));
        out.push('\n');
    }
    out
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
