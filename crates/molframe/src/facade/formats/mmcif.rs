//! mmCIF reading and writing, including the streaming write path `bcif` and
//! `modelcif` share.

use crate::facade::unsupported;
use molframe_core::diagnostic::{Code, Diagnostic, Findings};
use molframe_core::io::{Format, InputBuffer, Limits};
use molframe_core::structure::Structure;
use std::io;
use std::io::Write;
use std::path::Path;

/// Reads a document, preserving everything the file held.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_document(path: impl AsRef<Path>) -> Result<molframe_cif::Document, Findings> {
    let path = path.as_ref();
    let input = InputBuffer::open(path, Limits::default()).map_err(Findings::from)?;
    let name = path.file_name().and_then(|name| name.to_str());
    let format = Format::detect(Format::Auto, &input, name).map_err(Findings::from)?;
    match format {
        Format::Mmcif => molframe_cif::parse(&input)
            .map(|(document, _)| document)
            .map_err(Findings::from),
        Format::Pdbml => molframe_cif::parse_pdbml_document(input.as_bytes()).map_err(|error| {
            Findings::from(Diagnostic::new(Code::E1102).with_message(error.to_string()))
        }),
        #[cfg(feature = "bcif")]
        Format::BinaryCif => molframe_bcif::read_document(&input, Limits::default())
            .and_then(|document| document.to_document())
            .map_err(Findings::from),
        other => Err(unsupported(other).into()),
    }
}

/// Renders a structure as valid, self-consistent mmCIF in memory.
///
/// # Errors
///
/// Returns the same diagnostics the file verbs return.
pub fn write_mmcif(structure: &Structure) -> Result<String, Findings> {
    write_mmcif_with_options(structure, &molframe_cif::CifWriteOptions::new())
}

/// Renders mmCIF in memory with explicit identifier decisions.
///
/// # Errors
///
/// Returns the same diagnostics the file verbs return.
pub fn write_mmcif_with_options(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> Result<String, Findings> {
    let mut output = Vec::with_capacity(structure.atom_count() as usize * 100);
    write_mmcif_to_with_options(structure, options, &mut output)
        .map_err(|error| cif_write_findings(&error))?;
    String::from_utf8(output).map_err(|error| {
        cif_write_findings(&molframe_cif::CifWriteToError::Output(io::Error::new(
            io::ErrorKind::InvalidData,
            error,
        )))
    })
}

/// Streams canonical mmCIF with explicit identifier decisions.
///
/// The in-memory writers above are the public ladder; this streaming form
/// serves the file verbs and stays private so the facade exposes one write
/// convention per family.
///
/// # Errors
///
/// Returns canonical projection or destination errors.
pub(crate) fn write_mmcif_to_with_options<W: Write>(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
    output: &mut W,
) -> Result<(), molframe_cif::CifWriteToError> {
    #[cfg(feature = "modelcif")]
    if let Some(model) = structure
        .extensions()
        .get::<molframe_modelcif::ModelCif>(molframe_modelcif::MODEL_CIF_EXTENSION)
    {
        return molframe_modelcif::write_canonical_to(structure, model, options, output);
    }
    molframe_cif::write_canonical_to(structure, options, output)
}

pub(crate) fn cif_write_findings(error: &molframe_cif::CifWriteToError) -> Findings {
    match error {
        molframe_cif::CifWriteToError::Projection(error) => Diagnostic::new(Code::E4105)
            .with_message("structure cannot be projected to canonical CIF")
            .with_context("reason", error.to_string()),
        molframe_cif::CifWriteToError::Output(error) => Diagnostic::new(Code::E7901)
            .with_message("canonical CIF output failed")
            .with_context("reason", error.to_string()),
    }
    .into()
}
