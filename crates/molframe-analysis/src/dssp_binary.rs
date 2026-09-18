//! Explicit adapter for the external `mkdssp` executable.

use molframe_cif::Document;
use molframe_core::io::InputBuffer;
use std::path::Path;
use std::process::Command;

/// One contiguous secondary-structure segment emitted by DSSP 4 mmCIF output.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DsspSegment {
    /// mmCIF `_struct_conf.conf_type_id` value.
    pub kind: Box<str>,
    /// Beginning label chain.
    pub begin_chain: Box<str>,
    /// Beginning label sequence number.
    pub begin_sequence: i64,
    /// Ending label chain.
    pub end_chain: Box<str>,
    /// Ending label sequence number.
    pub end_sequence: i64,
}

/// Failure to execute or interpret the external DSSP program.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum DsspBinaryError {
    /// The process could not be started or waited for.
    #[error("could not execute DSSP: {0}")]
    Io(#[from] std::io::Error),
    /// DSSP exited unsuccessfully.
    #[error("DSSP exited unsuccessfully: {message}")]
    Exit {
        /// Standard error decoded lossily.
        message: String,
    },
    /// DSSP output was not valid mmCIF.
    #[error("DSSP output was not valid mmCIF")]
    Parse(Vec<molframe_core::Diagnostic>),
}

/// Runs a caller-selected DSSP executable on a caller-selected structure file.
///
/// No executable lookup or network access occurs until this function is called.
/// DSSP 4 writes annotated mmCIF to standard output when given one input path;
/// its `_struct_conf` rows are returned as typed segments.
///
/// # Errors
///
/// Returns [`DsspBinaryError`] when the process cannot run, exits unsuccessfully,
/// or produces malformed mmCIF.
pub fn run_dssp(
    executable: impl AsRef<Path>,
    input: impl AsRef<Path>,
) -> Result<Vec<DsspSegment>, DsspBinaryError> {
    let output = Command::new(executable.as_ref())
        .arg(input.as_ref())
        .output()?;
    if !output.status.success() {
        return Err(DsspBinaryError::Exit {
            message: String::from_utf8_lossy(&output.stderr).trim().to_owned(),
        });
    }
    parse_dssp_output(&output.stdout)
}

/// Parses DSSP 4's annotated mmCIF output without executing a process.
///
/// # Errors
///
/// Returns [`DsspBinaryError::Parse`] for malformed mmCIF.
pub fn parse_dssp_output(output: &[u8]) -> Result<Vec<DsspSegment>, DsspBinaryError> {
    let input = InputBuffer::from_bytes(output.to_vec());
    let (document, _) = molframe_cif::parse(&input).map_err(DsspBinaryError::Parse)?;
    Ok(segments(&document))
}

fn segments(document: &Document) -> Vec<DsspSegment> {
    let Some(category) = document
        .first_block()
        .and_then(|block| block.category("struct_conf"))
    else {
        return Vec::new();
    };
    (0..category.row_count())
        .filter_map(|row| {
            Some(DsspSegment {
                kind: category
                    .identifier("conf_type_id", row)?
                    .into_owned()
                    .into(),
                begin_chain: category
                    .identifier("beg_label_asym_id", row)?
                    .into_owned()
                    .into(),
                begin_sequence: category.value("beg_label_seq_id", row)?.as_integer()?,
                end_chain: category
                    .identifier("end_label_asym_id", row)?
                    .into_owned()
                    .into(),
                end_sequence: category.value("end_label_seq_id", row)?.as_integer()?,
            })
        })
        .collect()
}

#[cfg(test)]
#[path = "dssp_binary_tests.rs"]
mod tests;
