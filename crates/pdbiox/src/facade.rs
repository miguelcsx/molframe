//! Reading and writing, dispatched to whichever format crate is linked.
//!
//! The dispatch is the only thing here. Deciding which format a file is belongs
//! to the core, reading it belongs to the format crate, and this is the seam
//! between them.

use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::{Format, InputBuffer, Limits, ReadOptions, ReadResult};
use pdbiox_core::structure::Structure;
use std::path::Path;

/// Reads a structure, discarding what was wrong with the file.
///
/// Convenient, and less honest than [`read_with_diagnostics`]: a file that
/// parses is not the same as a file that is right.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read(path: impl AsRef<Path>) -> Result<Structure, Vec<Diagnostic>> {
    read_with_diagnostics(path).map(|(structure, _)| structure)
}

/// Reads a structure and everything that was wrong with the file it came from.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_with_diagnostics(
    path: impl AsRef<Path>,
) -> Result<(Structure, Vec<Diagnostic>), Vec<Diagnostic>> {
    read_with_options(path, &ReadOptions::new())
}

/// Reads a structure under explicit options.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_with_options(path: impl AsRef<Path>, options: &ReadOptions) -> ReadResult {
    let path = path.as_ref();
    let input = InputBuffer::open(path, options.limits).map_err(|finding| vec![finding])?;
    let name = path.file_name().and_then(|name| name.to_str());
    read_buffer(&input, name, options)
}

/// Reads a structure from bytes already in hand.
///
/// `name` is only used to settle the format when the content does not.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_bytes(bytes: Vec<u8>, name: Option<&str>, options: &ReadOptions) -> ReadResult {
    let input = InputBuffer::from_bytes(bytes);
    read_buffer(&input, name, options)
}

fn read_buffer(input: &InputBuffer, name: Option<&str>, options: &ReadOptions) -> ReadResult {
    let format = Format::detect(options.format, input, name).map_err(|finding| vec![finding])?;
    match format {
        #[cfg(feature = "mmcif")]
        Format::Mmcif => pdbiox_cif::read(input, options),
        #[cfg(feature = "pdb")]
        Format::Pdb => pdbiox_pdb::read(input, options),
        other => Err(vec![unsupported(other)]),
    }
}

/// Writes a structure in the legacy fixed-column format, or explains why it
/// cannot be written.
///
/// # Errors
///
/// Returns every capacity the structure exceeds.
#[cfg(feature = "pdb")]
pub fn write_pdb(
    structure: &Structure,
    options: &pdbiox_pdb::PdbOptions,
) -> Result<String, Vec<Diagnostic>> {
    pdbiox_pdb::write(structure, options)
}

/// The finding raised for a format this build cannot read.
fn unsupported(format: Format) -> Diagnostic {
    Diagnostic::new(Code::E1001)
        .with_message("this build does not include a reader for the detected format")
        .with_context("format", format.name())
        .with_context("crate", crate_for(format))
}

fn crate_for(format: Format) -> &'static str {
    match format {
        Format::Mmcif => "pdbiox-cif",
        Format::Pdb => "pdbiox-pdb",
        // The format list grows as crates land, and a build that has not been
        // rebuilt against the newer list should say so rather than fail to
        // compile against it.
        _ => "(not yet assigned)",
    }
}

/// Reads a document, preserving everything the file held.
///
/// # Errors
///
/// Returns the findings that stopped the read.
#[cfg(feature = "mmcif")]
pub fn read_document(path: impl AsRef<Path>) -> Result<pdbiox_cif::Document, Vec<Diagnostic>> {
    let path = path.as_ref();
    let input = InputBuffer::open(path, Limits::default()).map_err(|finding| vec![finding])?;
    pdbiox_cif::parse(&input).map(|(document, _)| document)
}

/// Writes a structure as a valid, self-consistent mmCIF file.
#[cfg(feature = "mmcif")]
#[must_use]
pub fn write_mmcif(structure: &Structure) -> String {
    pdbiox_cif::write_canonical(structure)
}

/// The ceilings a read runs under by default.
#[must_use]
pub fn default_limits() -> Limits {
    Limits::default()
}

#[cfg(test)]
#[path = "facade_tests.rs"]
mod tests;
