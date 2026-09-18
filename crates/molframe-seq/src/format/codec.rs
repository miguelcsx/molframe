//! Dispatch that delegates to each native Rust parser and writer.

use super::{SequenceDocument, SequenceFormat, SequenceFormatError};
use crate::{
    FastaRecord, Tree, parse_a2m, parse_a3m, parse_clustal, parse_fasta, parse_fastq, parse_phylip,
    parse_stockholm, write_a2m, write_a3m, write_clustal, write_fasta, write_fastq, write_phylip,
    write_stockholm,
};

/// Common reader contract for sequence serializations.
pub trait SequenceReader {
    /// Parses and validates the selected serialization.
    ///
    /// # Errors
    ///
    /// Returns [`SequenceFormatError`] rather than a partial document.
    fn read(&self, text: &str) -> Result<SequenceDocument, SequenceFormatError>;
}

/// Common writer contract for sequence serializations.
pub trait SequenceWriter {
    /// Validates and serializes one compatible document.
    ///
    /// # Errors
    ///
    /// Returns [`SequenceFormatError`] rather than partial text.
    fn write(&self, document: &SequenceDocument) -> Result<String, SequenceFormatError>;
}

/// Reads through the common dispatch surface.
///
/// # Errors
///
/// Returns [`SequenceFormatError`] rather than a partial document.
pub fn read_sequence(
    text: &str,
    format: SequenceFormat,
) -> Result<SequenceDocument, SequenceFormatError> {
    format.read(text)
}

/// Writes through the common dispatch surface.
///
/// # Errors
///
/// Returns [`SequenceFormatError`] rather than partial text.
pub fn write_sequence(
    document: &SequenceDocument,
    format: SequenceFormat,
) -> Result<String, SequenceFormatError> {
    format.write(document)
}

impl SequenceReader for SequenceFormat {
    fn read(&self, text: &str) -> Result<SequenceDocument, SequenceFormatError> {
        match self {
            Self::Fasta => read_fasta(text),
            Self::Fastq => Ok(SequenceDocument::Fastq(parse_fastq(text)?)),
            Self::Stockholm => read_stockholm(text),
            Self::Clustal => read_clustal(text),
            Self::Phylip => read_phylip(text),
            Self::A2m => Ok(SequenceDocument::Records(parse_a2m(text)?)),
            Self::A3m => Ok(SequenceDocument::Records(parse_a3m(text)?)),
            Self::Newick => Ok(SequenceDocument::Tree(Tree::from_newick(text)?)),
        }
    }
}

impl SequenceWriter for SequenceFormat {
    fn write(&self, document: &SequenceDocument) -> Result<String, SequenceFormatError> {
        match (self, document) {
            (Self::Fasta, SequenceDocument::Records(records)) => Ok(write_fasta(records)),
            (Self::Fastq, SequenceDocument::Fastq(records)) => Ok(write_fastq(records)?),
            (Self::Stockholm, SequenceDocument::Records(records)) => {
                validate_alignment(records, *self)?;
                Ok(write_stockholm(records))
            }
            (Self::Clustal, SequenceDocument::Records(records)) => {
                validate_alignment(records, *self)?;
                Ok(write_clustal(records))
            }
            (Self::Phylip, SequenceDocument::Records(records)) => {
                validate_alignment(records, *self)?;
                Ok(write_phylip(records))
            }
            (Self::A2m, SequenceDocument::Records(records)) => Ok(write_a2m(records)?),
            (Self::A3m, SequenceDocument::Records(records)) => Ok(write_a3m(records)?),
            (Self::Newick, SequenceDocument::Tree(tree)) => Ok(tree.to_newick()),
            _ => Err(SequenceFormatError::WrongDocument { format: *self }),
        }
    }
}

fn read_fasta(text: &str) -> Result<SequenceDocument, SequenceFormatError> {
    let mut saw_header = false;
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if header.trim().is_empty() {
                return Err(invalid(SequenceFormat::Fasta, "empty record identifier"));
            }
            saw_header = true;
        } else if !saw_header {
            return Err(invalid(
                SequenceFormat::Fasta,
                "missing leading FASTA header",
            ));
        }
    }
    if !saw_header {
        return Err(invalid(
            SequenceFormat::Fasta,
            "missing leading FASTA header",
        ));
    }
    let records = parse_fasta(text);
    if records.iter().any(|record| record.id.is_empty()) {
        return Err(invalid(SequenceFormat::Fasta, "empty record identifier"));
    }
    Ok(SequenceDocument::Records(records))
}

fn read_stockholm(text: &str) -> Result<SequenceDocument, SequenceFormatError> {
    if text.lines().next() != Some("# STOCKHOLM 1.0")
        || !text.lines().any(|line| line.trim() == "//")
    {
        return Err(invalid(
            SequenceFormat::Stockholm,
            "missing header or terminator",
        ));
    }
    alignment_document(parse_stockholm(text), SequenceFormat::Stockholm)
}

fn read_clustal(text: &str) -> Result<SequenceDocument, SequenceFormatError> {
    if !text
        .lines()
        .next()
        .is_some_and(|line| line.to_ascii_uppercase().starts_with("CLUSTAL"))
    {
        return Err(invalid(SequenceFormat::Clustal, "missing CLUSTAL header"));
    }
    alignment_document(parse_clustal(text), SequenceFormat::Clustal)
}

fn read_phylip(text: &str) -> Result<SequenceDocument, SequenceFormatError> {
    let Some(header) = text.lines().find(|line| !line.trim().is_empty()) else {
        return Err(invalid(
            SequenceFormat::Phylip,
            "missing count/length header",
        ));
    };
    let mut fields = header.split_whitespace();
    let count = fields.next().and_then(|value| value.parse::<usize>().ok());
    let length = fields.next().and_then(|value| value.parse::<usize>().ok());
    let (Some(count), Some(length), None) = (count, length, fields.next()) else {
        return Err(invalid(
            SequenceFormat::Phylip,
            "invalid count/length header",
        ));
    };
    let records = parse_phylip(text);
    if records.len() != count || records.iter().any(|record| record.sequence.len() != length) {
        return Err(invalid(
            SequenceFormat::Phylip,
            "declared dimensions do not match records",
        ));
    }
    alignment_document(records, SequenceFormat::Phylip)
}

fn alignment_document(
    records: Vec<FastaRecord>,
    format: SequenceFormat,
) -> Result<SequenceDocument, SequenceFormatError> {
    validate_alignment(&records, format)?;
    Ok(SequenceDocument::Records(records))
}

fn validate_alignment(
    records: &[FastaRecord],
    format: SequenceFormat,
) -> Result<(), SequenceFormatError> {
    let Some(first) = records.first() else {
        return Err(invalid(format, "alignment has no records"));
    };
    if first.id.is_empty() || records.iter().any(|record| record.id.is_empty()) {
        return Err(invalid(format, "alignment has an empty identifier"));
    }
    if records
        .iter()
        .any(|record| record.sequence.len() != first.sequence.len())
    {
        return Err(invalid(format, "alignment rows have different lengths"));
    }
    Ok(())
}

const fn invalid(format: SequenceFormat, reason: &'static str) -> SequenceFormatError {
    SequenceFormatError::Invalid { format, reason }
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
