//! Reading and writing, dispatched to whichever format crate is linked.
//!
//! The dispatch is the only thing here. Deciding which format a file is belongs
//! to the core, reading it belongs to the format crate, and this is the seam
//! between them.

use molframe_core::diagnostic::{Code, Diagnostic, Findings};
use molframe_core::io::{Format, InputBuffer, Limits, OutputOptions, OutputSink, ReadOptions};
use molframe_core::structure::Structure;
#[cfg(feature = "mmcif")]
use std::io;
use std::io::Write;
use std::path::Path;

#[cfg(feature = "mmcif")]
mod cif_family;
#[cfg(feature = "mmcif")]
mod extensions;
// The enum holds one variant per format crate, so with none of them linked it
// would be an empty type whose `next_batch` match has no arms to reach. The
// module is gated by the same four features that gate its variants, rather
// than publishing an enum nothing can construct.
#[cfg(any(
    feature = "mmcif",
    feature = "pdb",
    feature = "bcif",
    feature = "modelcif"
))]
mod structure_batches;

#[cfg(any(
    feature = "mmcif",
    feature = "pdb",
    feature = "bcif",
    feature = "modelcif"
))]
pub use structure_batches::{StructureBatchReader, open_structure_batches};

#[cfg(feature = "bcif")]
use cif_family::read_bcif_buffer;
#[cfg(feature = "mmcif")]
use cif_family::read_mmcif_buffer;

#[cfg(feature = "chem")]
use molframe_core::contract::DictionaryVersion;

#[cfg(feature = "geom")]
use molframe_core::selection::AtomSelection;
#[cfg(feature = "geom")]
use molframe_geom::Rigid;

/// Reads a structure, discarding what was wrong with the file.
///
/// Convenient, and less honest than [`read_with_diagnostics`]: a file that
/// parses is not the same as a file that is right.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read(path: impl AsRef<Path>) -> Result<Structure, Findings> {
    read_with_diagnostics(path).map(|(structure, _)| structure)
}

/// Reads a structure and everything that was wrong with the file it came from.
///
/// The findings that came with a structure that parsed are the pair's second
/// element. The findings that stopped a read are the error.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_with_diagnostics(
    path: impl AsRef<Path>,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    read_with_options(path, &ReadOptions::new())
}

/// Reads a structure under explicit options.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_with_options(
    path: impl AsRef<Path>,
    options: &ReadOptions,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    let path = path.as_ref();
    let input = InputBuffer::open(path, options.limits).map_err(Findings::from)?;
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
pub fn read_bytes(
    bytes: Vec<u8>,
    name: Option<&str>,
    options: &ReadOptions,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    let input = InputBuffer::from_bytes(bytes);
    read_buffer(&input, name, options)
}

/// Writes a structure in the format named by the destination.
///
/// A `.gz` or `.zst` suffix applies deterministic compression after the format
/// suffix has selected the writer.
///
/// # Errors
///
/// Returns a diagnostic when the name does not select a linked writer, the
/// structure cannot be represented, compression is unavailable, or the
/// destination cannot be written.
pub fn write(path: impl AsRef<Path>, structure: &Structure) -> Result<(), Findings> {
    write_with_options(path, structure, OutputOptions::default())
}

/// Streams a structure under an explicit output-memory policy.
///
/// The default policy is 100 MB and no policy above 500 MB is accepted. The
/// file sink itself retains at most 64 KiB; format writers receive rows or
/// columns incrementally and never construct a complete file buffer.
///
/// # Errors
///
/// Returns format, projection, memory-policy, compression, or destination
/// diagnostics.
pub fn write_with_options(
    path: impl AsRef<Path>,
    structure: &Structure,
    options: OutputOptions,
) -> Result<(), Findings> {
    let path = path.as_ref();
    let name = path.file_name().and_then(|name| name.to_str());
    let Some(format) = name.and_then(Format::from_name) else {
        return Err(Findings::from(
            Diagnostic::new(Code::E1001).with_context("name", path.display().to_string()),
        ));
    };
    let mut output = OutputSink::create(path, options).map_err(Findings::from)?;
    write_stream(&mut output, structure, format, options.memory_limit_bytes)?;
    output.finish().map_err(Findings::from)
}

fn write_stream<W: Write>(
    output: &mut W,
    structure: &Structure,
    format: Format,
    memory_limit_bytes: usize,
) -> Result<(), Findings> {
    #[cfg(not(feature = "bcif"))]
    let _ = memory_limit_bytes;
    #[cfg(not(any(feature = "mmcif", feature = "bcif", feature = "pdb")))]
    let _ = (output, structure);
    match format {
        #[cfg(feature = "mmcif")]
        Format::Mmcif => {
            write_mmcif_to_with_options(structure, &molframe_cif::CifWriteOptions::new(), output)
                .map_err(|error| cif_write_findings(&error))
        }
        #[cfg(feature = "bcif")]
        Format::BinaryCif => molframe_bcif::write_structure_to_with_memory_limit(
            structure,
            &molframe_cif::CifWriteOptions::new(),
            memory_limit_bytes,
            output,
        )
        .map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pdb => molframe_pdb::write_to(structure, &molframe_pdb::PdbOptions::new(), output)
            .map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Mmtf => molframe_pdb::write_mmtf_to(structure, output).map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pqr => {
            molframe_pdb::write_pqr_to(structure, &molframe_pdb::PdbOptions::new(), output)
                .map_err(Findings::from)
        }
        #[cfg(feature = "pdb")]
        Format::Pdbqt => {
            molframe_pdb::write_pdbqt_to(structure, &molframe_pdb::PdbOptions::new(), output)
                .map_err(Findings::from)
        }
        other => Err(unsupported_writer(other).into()),
    }
}

fn read_buffer(
    input: &InputBuffer,
    name: Option<&str>,
    options: &ReadOptions,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    let format = Format::detect(options.format, input, name).map_err(Findings::from)?;
    match format {
        #[cfg(feature = "mmcif")]
        Format::Mmcif => read_mmcif_buffer(input, options),
        #[cfg(feature = "mmcif")]
        Format::Pdbml => read_pdbml_buffer(input, options),
        #[cfg(feature = "bcif")]
        Format::BinaryCif => read_bcif_buffer(input, options),
        #[cfg(feature = "pdb")]
        Format::Pdb => molframe_pdb::read(input, options).map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Mmtf => molframe_pdb::read_mmtf(input, options).map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pqr => molframe_pdb::read_pqr(input, options).map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pdbqt => molframe_pdb::read_pdbqt(input, options).map_err(Findings::from),
        other => Err(unsupported(other).into()),
    }
}

#[cfg(feature = "mmcif")]
fn read_pdbml_buffer(
    input: &InputBuffer,
    options: &ReadOptions,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    // `PdbmlReadError` renders its `Pdbml` variant by delegating to the inner
    // error, so the catch-all arm below produces the same finding the variant
    // arm would have.
    match molframe_cif::read_pdbml(input.as_bytes(), options) {
        Ok((document, structure, findings)) => {
            extensions::attach_materialized_cif_metadata(&document, structure, findings, options)
                .map_err(Findings::from)
        }
        Err(molframe_cif::PdbmlReadError::Findings(findings)) => Err(findings.into()),
        Err(error) => Err(Findings::from(
            Diagnostic::new(Code::E1102)
                .with_message("PDBML/XML could not be decoded")
                .with_context("decoder", error.to_string()),
        )),
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
    options: &molframe_pdb::PdbOptions,
) -> Result<String, Findings> {
    molframe_pdb::write(structure, options).map_err(Findings::from)
}

/// The finding raised for a format this build cannot read.
fn unsupported(format: Format) -> Diagnostic {
    Diagnostic::new(Code::E1001)
        .with_message("this build does not include a reader for the detected format")
        .with_context("format", format.name())
        .with_context("crate", crate_for(format))
}

fn unsupported_writer(format: Format) -> Diagnostic {
    Diagnostic::new(Code::E1001)
        .with_message("this format has no bounded incremental writer")
        .with_context("format", format.name())
        .with_context("crate", crate_for(format))
}

fn crate_for(format: Format) -> &'static str {
    match format {
        Format::Mmcif | Format::Pdbml => "molframe-cif",
        Format::BinaryCif => "molframe-bcif",
        Format::Pdb | Format::Pqr | Format::Pdbqt | Format::Mmtf => "molframe-pdb",
        // `Format` is non-exhaustive across crate versions. An unknown variant
        // is unsupported by this compiled facade rather than assigned a guessed
        // owner.
        _ => "unlinked-format",
    }
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
/// Returns the canonical preflight or in-memory destination error.
#[cfg(feature = "mmcif")]
pub fn write_mmcif(structure: &Structure) -> Result<String, molframe_cif::CifWriteToError> {
    write_mmcif_with_options(structure, &molframe_cif::CifWriteOptions::new())
}

/// Renders mmCIF in memory with explicit identifier decisions.
///
/// # Errors
///
/// Returns the canonical preflight or in-memory destination error.
#[cfg(feature = "mmcif")]
pub fn write_mmcif_with_options(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> Result<String, molframe_cif::CifWriteToError> {
    let mut output = Vec::with_capacity(structure.atom_count() as usize * 100);
    write_mmcif_to_with_options(structure, options, &mut output)?;
    String::from_utf8(output).map_err(|error| {
        molframe_cif::CifWriteToError::Output(io::Error::new(io::ErrorKind::InvalidData, error))
    })
}

/// Streams canonical mmCIF to a byte destination.
///
/// # Errors
///
/// Returns canonical projection or destination errors.
#[cfg(feature = "mmcif")]
pub fn write_mmcif_to<W: Write>(
    structure: &Structure,
    output: &mut W,
) -> Result<(), molframe_cif::CifWriteToError> {
    write_mmcif_to_with_options(structure, &molframe_cif::CifWriteOptions::new(), output)
}

/// Streams canonical mmCIF with explicit identifier decisions.
///
/// # Errors
///
/// Returns canonical projection or destination errors.
#[cfg(feature = "mmcif")]
pub fn write_mmcif_to_with_options<W: Write>(
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
fn cif_write_findings(error: &molframe_cif::CifWriteToError) -> Findings {
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

#[cfg(test)]
#[path = "facade/facade_tests.rs"]
mod tests;
