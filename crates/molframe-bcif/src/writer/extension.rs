//! Structure writing with extra categories supplied as a parsed document.

use super::document::encode_category;
use super::structure::write_with_extra;
use molframe_cif::{CifWriteOptions, Document};
use molframe_core::diagnostic::Diagnostic;
use molframe_core::io::DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES;
use molframe_core::structure::Structure;

/// Writes a structure and appends every category of `extension` after its coordinates.
///
/// The extension document carries categories the canonical projection does not
/// know, such as the `ma_*` tables of a model archive file. Each is encoded
/// with the same smallest-lossless column codecs as [`super::write_document`].
///
/// # Errors
///
/// Returns canonical projection or `BinaryCIF` encoding diagnostics.
pub fn write_structure_with_extension(
    structure: &Structure,
    options: &CifWriteOptions,
    extension: &Document,
) -> Result<Vec<u8>, Vec<Diagnostic>> {
    let mut extra = Vec::new();
    for block in extension.blocks() {
        for category in block.categories() {
            extra.push(encode_category(category).map_err(|error| vec![error])?);
        }
    }
    let mut output = Vec::new();
    write_with_extra(
        structure,
        options,
        DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES,
        extra,
        &mut output,
    )?;
    Ok(output)
}
