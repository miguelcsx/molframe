//! The entry point: bytes in, structure and findings out.

use crate::document::Document;
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
    read_with_document(input, options).map(|(_, structure, findings)| (structure, findings))
}

/// Reads both the lossless document and its normalised structure in one parse.
///
/// This is the integration seam for domain crates that interpret additional
/// categories without reparsing the input.
///
/// # Errors
///
/// Returns the findings that stopped parsing or lowering.
pub fn read_with_document(
    input: &InputBuffer,
    options: &ReadOptions,
) -> Result<
    (
        Document,
        pdbiox_core::Structure,
        Vec<pdbiox_core::Diagnostic>,
    ),
    Vec<pdbiox_core::Diagnostic>,
> {
    let (document, mut findings) = parse(input)?;
    let (structure, interpretation) = lower(&document, options)?;
    findings.extend(interpretation);
    options
        .finish(structure, findings)
        .map(|(structure, findings)| (document, structure, findings))
}

#[cfg(test)]
#[path = "structure_tests.rs"]
mod tests;
