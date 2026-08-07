//! `BinaryCIF` container reading followed by shared mmCIF lowering.

use crate::BinaryDocument;
use pdbiox_cif::Document;
use pdbiox_core::diagnostic::Diagnostic;
use pdbiox_core::io::{Format, InputBuffer, Limits, ReadOptions, ReadResult, Reader};

/// Reader for `BinaryCIF`.
#[derive(Clone, Copy, Debug, Default)]
pub struct BcifReader;

impl Reader for BcifReader {
    const FORMAT: Format = Format::BinaryCif;

    fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
        read(input, options)
    }
}

/// Opens a `BinaryCIF` document without decoding its columns.
///
/// # Errors
///
/// Returns a registered container or size diagnostic.
pub fn read_document(input: &InputBuffer, limits: Limits) -> Result<BinaryDocument, Diagnostic> {
    BinaryDocument::parse(input.as_bytes(), limits)
}

/// Reads a `BinaryCIF` structure through the shared `Document` lowering path.
///
/// # Errors
///
/// Returns ordered container, codec, or interpretation diagnostics.
pub fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    read_with_document(input, options).map(|(_, structure, findings)| (structure, findings))
}

/// Reads both the decoded shared document and its normalised structure once.
///
/// # Errors
///
/// Returns ordered container, codec, or interpretation diagnostics.
pub fn read_with_document(
    input: &InputBuffer,
    options: &ReadOptions,
) -> Result<(Document, pdbiox_core::Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
    let binary = read_document(input, options.limits).map_err(|finding| vec![finding])?;
    let document: Document = binary.to_document().map_err(|finding| vec![finding])?;
    let (structure, findings) = pdbiox_cif::lower(&document, options)?;
    options
        .finish(structure, findings)
        .map(|(structure, findings)| (document, structure, findings))
}

#[cfg(test)]
#[path = "reader_tests.rs"]
mod tests;
