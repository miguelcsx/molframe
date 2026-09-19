//! The per-format convenience verbs: the readers and writers a caller names
//! by format, which the facade's dispatch does not need in line.
//!
//! The file ladder (`read`, `write`, `write_with_options`) lives in the parent and
//! streams through the same kernels these wrap.

use molframe_core::diagnostic::{Code, Diagnostic, Findings};
use molframe_core::io::{Format, InputBuffer, Limits};
use molframe_core::structure::Structure;
#[cfg(feature = "mmcif")]
use std::io;
#[cfg(feature = "mmcif")]
use std::io::Write;
use std::path::Path;

use super::unsupported;
#[cfg(feature = "chem")]
use molframe_core::contract::DictionaryVersion;
#[cfg(feature = "geom")]
use molframe_core::selection::AtomSelection;
#[cfg(feature = "geom")]
use molframe_geom::Rigid;

/// Writes a structure in the legacy fixed-column format, or explains why it
/// cannot be written.
///
/// # Errors
///
/// Returns every capacity the structure exceeds.
#[cfg(feature = "pdb")]
pub fn write_pdb(
    structure: &Structure,
    options: &molframe_pdb::PdbOptions,
) -> Result<String, Findings> {
    molframe_pdb::write(structure, options).map_err(Findings::from)
}

/// Reads a document, preserving everything the file held.
///
/// # Errors
///
/// Returns the findings that stopped the read.
#[cfg(feature = "mmcif")]
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

/// Reads a Chemical Component Dictionary file as a versioned provider.
///
/// This explicit entry point avoids guessing whether a `.cif` file is a
/// structure or a component dictionary. Parsing and component lowering remain
/// in their owning Rust crates.
///
/// # Errors
///
/// Returns file, CIF syntax, or component-definition diagnostics.
#[cfg(feature = "chem")]
pub fn read_component_dictionary(
    path: impl AsRef<Path>,
    version: DictionaryVersion,
) -> Result<(molframe_chem::CifProvider, Vec<Diagnostic>), Findings> {
    let input = InputBuffer::open(path, Limits::default()).map_err(Findings::from)?;
    molframe_chem::read_ccd(&input, version).map_err(Findings::from)
}

/// Renders a structure as valid, self-consistent mmCIF in memory.
///
/// # Errors
///
/// Returns the same diagnostics the file verbs return.
#[cfg(feature = "mmcif")]
pub fn write_mmcif(structure: &Structure) -> Result<String, Findings> {
    write_mmcif_with_options(structure, &molframe_cif::CifWriteOptions::new())
}

/// Renders mmCIF in memory with explicit identifier decisions.
///
/// # Errors
///
/// Returns the same diagnostics the file verbs return.
#[cfg(feature = "mmcif")]
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
#[cfg(feature = "mmcif")]
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

/// Renders deterministic `BinaryCIF` bytes in memory.
///
/// # Errors
///
/// Returns a diagnostic if a projected column cannot be represented.
#[cfg(feature = "bcif")]
pub fn write_bcif(structure: &Structure) -> Result<Vec<u8>, Findings> {
    molframe_bcif::write_structure(structure).map_err(Findings::from)
}

/// Renders deterministic `BinaryCIF` in memory with explicit identifier decisions.
///
/// # Errors
///
/// Returns a diagnostic if canonical preflight or binary encoding fails.
#[cfg(feature = "bcif")]
pub fn write_bcif_with_options(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> Result<Vec<u8>, Findings> {
    molframe_bcif::write_structure_with_options(structure, options).map_err(Findings::from)
}

#[cfg(feature = "mmcif")]
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

/// Applies one rigid transform to selected atoms in every dense model.
///
/// The coordinate transaction owns storage and generation tracking; the
/// geometric crate owns the transform arithmetic. This facade only joins the
/// two capabilities and publishes the resulting immutable snapshot.
///
/// # Errors
///
/// Returns a diagnostic when a selected atom is outside the topology, the
/// structure is a ragged ensemble, or the transformed snapshot is invalid.
#[cfg(feature = "geom")]
pub fn transform(
    structure: &Structure,
    selection: &AtomSelection,
    rigid: &Rigid,
) -> Result<Structure, Findings> {
    if structure.ragged_models().is_some() {
        return Err(Diagnostic::new(Code::E6008).into());
    }
    if let Some(atom) = selection
        .iter()
        .find(|atom| *atom >= structure.atom_count())
    {
        return Err(Diagnostic::new(Code::E6009)
            .with_context("atom", atom.to_string())
            .into());
    }

    let mut editor = structure.edit();
    editor
        .transform(selection, |position| rigid.apply(position))
        .map_err(Findings::from)?;
    editor.commit().map_err(Findings::from)
}
