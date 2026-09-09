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
    super::direct::read(input, options)
}

/// Reads a structure while retaining selected non-coordinate categories.
///
/// `keep_category` is applied in addition to the bounded set required by core
/// structure lowering. The returned document preserves the selected categories'
/// values and order but is not a lossless representation of the source. Use
/// [`read_with_document`] when every category must survive a round trip.
///
/// # Errors
///
/// Returns the findings that stopped parsing or lowering.
#[doc(hidden)]
pub fn read_with_metadata(
    input: &InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
) -> Result<
    (
        Document,
        pdbiox_core::Structure,
        Vec<pdbiox_core::Diagnostic>,
    ),
    Vec<pdbiox_core::Diagnostic>,
> {
    super::direct::read_with_metadata(input, options, keep_category)
}

/// Result of a combined structure read and external event projection.
#[doc(hidden)]
pub type ProjectedReadResult<T> = Result<
    (
        Document,
        pdbiox_core::Structure,
        T,
        Vec<pdbiox_core::Diagnostic>,
    ),
    Vec<pdbiox_core::Diagnostic>,
>;

/// Reads a structure and an external event projection in one lexer traversal.
///
/// The projection receives only categories it accepts and never requires a
/// lossless document for those values. Selected metadata for other domain
/// extensions remains available in the returned document.
///
/// # Errors
///
/// Returns findings that stopped CIF parsing or structure lowering.
#[doc(hidden)]
pub fn read_with_projection<S>(
    input: &InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
    projection: S,
) -> ProjectedReadResult<S::Output>
where
    S: crate::CifEventSink,
{
    super::direct::read_with_projection(input, options, keep_category, projection)
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
