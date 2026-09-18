//! `BinaryCIF` container and structure reading.

use crate::BinaryDocument;
use molframe_cif::Document;
use molframe_core::diagnostic::Diagnostic;
use molframe_core::io::{Format, InputBuffer, Limits, ReadOptions, ReadResult, Reader};

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

/// Reads a `BinaryCIF` structure without materialising a coordinate DOM.
///
/// # Errors
///
/// Returns ordered container, codec, or interpretation diagnostics.
pub fn read(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    super::direct::read(input, options)
}

/// Reads a structure while retaining selected non-coordinate categories.
///
/// The returned document is a bounded metadata projection, not a lossless
/// representation. Use [`read_with_document`] when every category must survive
/// a round trip.
///
/// # Errors
///
/// Returns ordered container, codec, or interpretation diagnostics.
#[doc(hidden)]
pub fn read_with_metadata(
    input: &InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
) -> Result<(Document, molframe_core::Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
    super::direct::read_with_metadata(input, options, keep_category)
}

/// Result of a direct structure read and compact non-coordinate projection.
#[doc(hidden)]
pub type ProjectedReadResult<T> = Result<
    (
        Document,
        molframe_core::Structure,
        T,
        Vec<molframe_core::Diagnostic>,
    ),
    Vec<molframe_core::Diagnostic>,
>;

/// Reads a structure and an external compact projection from one container decode.
///
/// Projected categories are decoded one column at a time and are not retained
/// in the returned metadata document unless `keep_category` also requests
/// them. Coordinate projection remains owned by the direct structure reader.
///
/// # Errors
///
/// Returns ordered container, codec, or interpretation diagnostics.
#[doc(hidden)]
pub fn read_with_projection<S>(
    input: &InputBuffer,
    options: &ReadOptions,
    keep_category: fn(&str) -> bool,
    projection: S,
) -> ProjectedReadResult<S::Output>
where
    S: molframe_cif::CifEventSink,
{
    super::direct::read_with_projection(input, options, keep_category, projection)
}

/// Reads both the decoded shared document and its normalised structure once.
///
/// # Errors
///
/// Returns ordered container, codec, or interpretation diagnostics.
pub fn read_with_document(
    input: &InputBuffer,
    options: &ReadOptions,
) -> Result<(Document, molframe_core::Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
    let binary = read_document(input, options.limits).map_err(|finding| vec![finding])?;
    let document: Document = binary.to_document().map_err(|finding| vec![finding])?;
    let (structure, findings) = molframe_cif::lower(&document, options)?;
    options
        .finish(structure, findings)
        .map(|(structure, findings)| (document, structure, findings))
}

#[cfg(test)]
#[path = "structure_tests.rs"]
mod tests;
