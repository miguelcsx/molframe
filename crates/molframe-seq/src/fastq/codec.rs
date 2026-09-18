//! Strict FASTQ reading and writing with byte-for-byte quality preservation.

/// One FASTQ sequence and its encoded quality scores.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastqRecord {
    /// First whitespace-delimited header token.
    pub id: String,
    /// Remaining header text.
    pub description: String,
    /// Sequence bytes.
    pub sequence: Vec<u8>,
    /// ASCII-encoded qualities, one per sequence byte.
    pub quality: Vec<u8>,
}

/// Why FASTQ text cannot be represented losslessly.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum FastqError {
    /// A record did not begin with an `@` header.
    MissingHeader {
        /// One-based source line.
        line: usize,
    },
    /// A record did not contain a `+` separator.
    MissingSeparator {
        /// One-based source line.
        line: usize,
    },
    /// Header identifier is empty.
    EmptyIdentifier {
        /// One-based source line or output record position.
        line: usize,
    },
    /// Sequence and quality byte counts differ.
    LengthMismatch {
        /// Zero-based record position.
        record: usize,
        /// Sequence byte count.
        sequence: usize,
        /// Quality byte count.
        quality: usize,
    },
    /// A sequence or quality contains a forbidden control byte.
    InvalidByte {
        /// Zero-based record position.
        record: usize,
        /// Zero-based byte position within the field.
        offset: usize,
    },
}

impl std::fmt::Display for FastqError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingHeader { line } => {
                write!(formatter, "FASTQ header missing at line {line}")
            }
            Self::MissingSeparator { line } => {
                write!(formatter, "FASTQ separator missing at line {line}")
            }
            Self::EmptyIdentifier { line } => {
                write!(formatter, "FASTQ identifier empty at line {line}")
            }
            Self::LengthMismatch {
                record,
                sequence,
                quality,
            } => write!(
                formatter,
                "FASTQ record {record} has {sequence} sequence and {quality} quality bytes"
            ),
            Self::InvalidByte { record, offset } => {
                write!(
                    formatter,
                    "FASTQ record {record} has invalid byte at {offset}"
                )
            }
        }
    }
}

impl std::error::Error for FastqError {}

/// Parses wrapped or four-line FASTQ records without guessing truncated data.
///
/// Sequence lines continue until `+`; quality lines continue until exactly the
/// sequence length. An optional identifier after `+` must equal the header ID.
///
/// # Errors
///
/// Returns [`FastqError`] for malformed framing, invalid bytes or unequal
/// sequence and quality lengths.
pub fn parse_fastq(text: &str) -> Result<Vec<FastqRecord>, FastqError> {
    let lines = text.lines().collect::<Vec<_>>();
    let mut cursor = 0usize;
    let mut records = Vec::new();
    while cursor < lines.len() {
        let header_line = cursor + 1;
        let Some(header) = lines[cursor].strip_prefix('@') else {
            return Err(FastqError::MissingHeader { line: header_line });
        };
        let (id, description) = split_header(header, header_line)?;
        cursor += 1;
        let mut sequence = Vec::new();
        while cursor < lines.len() && !lines[cursor].starts_with('+') {
            sequence.extend_from_slice(lines[cursor].as_bytes());
            cursor += 1;
        }
        let Some(separator) = lines.get(cursor) else {
            return Err(FastqError::MissingSeparator { line: cursor + 1 });
        };
        let separator_id = separator[1..].trim();
        if !separator_id.is_empty() && separator_id != id {
            return Err(FastqError::MissingSeparator { line: cursor + 1 });
        }
        cursor += 1;
        let mut quality = Vec::with_capacity(sequence.len());
        while cursor < lines.len() && quality.len() < sequence.len() {
            quality.extend_from_slice(lines[cursor].as_bytes());
            cursor += 1;
        }
        let record = records.len();
        validate_bytes(&sequence, record)?;
        validate_bytes(&quality, record)?;
        if sequence.len() != quality.len() {
            return Err(FastqError::LengthMismatch {
                record,
                sequence: sequence.len(),
                quality: quality.len(),
            });
        }
        records.push(FastqRecord {
            id,
            description,
            sequence,
            quality,
        });
    }
    Ok(records)
}

/// Writes canonical four-line FASTQ after validating every record.
///
/// # Errors
///
/// Returns [`FastqError`] for empty identifiers, forbidden control bytes or
/// unequal sequence and quality lengths.
pub fn write_fastq(records: &[FastqRecord]) -> Result<String, FastqError> {
    let mut output = String::new();
    for (record, entry) in records.iter().enumerate() {
        if entry.id.is_empty() || entry.id.bytes().any(|byte| byte.is_ascii_whitespace()) {
            return Err(FastqError::EmptyIdentifier { line: record + 1 });
        }
        validate_bytes(&entry.sequence, record)?;
        validate_bytes(&entry.quality, record)?;
        if entry.sequence.len() != entry.quality.len() {
            return Err(FastqError::LengthMismatch {
                record,
                sequence: entry.sequence.len(),
                quality: entry.quality.len(),
            });
        }
        output.push('@');
        output.push_str(&entry.id);
        if !entry.description.is_empty() {
            output.push(' ');
            output.push_str(&entry.description);
        }
        output.push('\n');
        push_ascii(&mut output, &entry.sequence);
        output.push_str("\n+\n");
        push_ascii(&mut output, &entry.quality);
        output.push('\n');
    }
    Ok(output)
}

fn split_header(header: &str, line: usize) -> Result<(String, String), FastqError> {
    let header = header.trim();
    if header.is_empty() {
        return Err(FastqError::EmptyIdentifier { line });
    }
    Ok(match header.split_once(char::is_whitespace) {
        Some((id, description)) => (id.to_owned(), description.trim().to_owned()),
        None => (header.to_owned(), String::new()),
    })
}

fn validate_bytes(bytes: &[u8], record: usize) -> Result<(), FastqError> {
    match bytes.iter().position(u8::is_ascii_control) {
        Some(offset) => Err(FastqError::InvalidByte { record, offset }),
        None => Ok(()),
    }
}

fn push_ascii(output: &mut String, bytes: &[u8]) {
    output.extend(bytes.iter().map(|byte| char::from(*byte)));
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
