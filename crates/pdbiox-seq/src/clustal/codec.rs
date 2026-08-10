//! Reading and writing the Clustal alignment format.
//!
//! Clustal writes an alignment in blocks: a header line, then repeated blocks of
//! `name sequence` lines followed by a conservation line of spaces and asterisks.
//! Parsing skips the header and the conservation lines — which start with
//! whitespace rather than a name — and concatenates each name's blocks in
//! first-seen order. Writing emits one block, so a written alignment parses back
//! to the same records.

use crate::fasta::FastaRecord;

/// Parses Clustal text into alignment records in first-seen order.
///
/// The header line and the conservation lines are ignored; an optional trailing
/// residue count on a sequence line is ignored too.
#[must_use]
pub fn parse_clustal(text: &str) -> Vec<FastaRecord> {
    let mut records: Vec<FastaRecord> = Vec::new();
    for line in text.lines() {
        if line.starts_with("CLUSTAL") || line.trim().is_empty() {
            continue;
        }
        if line.starts_with(|character: char| character.is_whitespace()) {
            continue;
        }
        let mut parts = line.split_whitespace();
        let (Some(name), Some(sequence)) = (parts.next(), parts.next()) else {
            continue;
        };
        let residues = sequence.bytes();
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

/// Writes alignment records as a single-block Clustal file.
#[must_use]
pub fn write_clustal(records: &[FastaRecord]) -> String {
    let mut out = String::from("CLUSTAL W (pdbiox) multiple sequence alignment\n\n");
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
