//! The entry point: bytes in, structure and findings out.

use crate::lower::lower;
use crate::parser::parse;
use pdbiox_core::io::{Format, InputBuffer, ReadOptions, ReadResult, Reader};

/// The reader for PDBx/mmCIF.
#[derive(Clone, Copy, Debug, Default)]
pub struct CifReader;

impl Reader for CifReader {
    const FORMAT: Format = Format::Mmcif;

    fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
        read(input, options)
    }
}

/// Reads a structure.
///
/// Findings from both steps are returned together: what was wrong with the file
/// as text, and what had to be decided while interpreting it.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    let (document, mut findings) = parse(input)?;
    let (structure, interpretation) = lower(&document, options)?;
    findings.extend(interpretation);
    Ok((structure, findings))
}

#[cfg(test)]
#[path = "reader_tests.rs"]
mod tests;
