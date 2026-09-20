//! Reading and writing, dispatched to whichever format crate is linked.
//!
//! The dispatch is the only thing here. Deciding which format a file is belongs
//! to the core, reading it belongs to the format crate, and this is the seam
//! between them.

#[cfg(feature = "mmcif")]
use formats::{cif_write_findings, write_mmcif_to_with_options};
use molframe_core::diagnostic::{Code, Diagnostic, Findings};
use molframe_core::io::{Format, InputBuffer, OutputOptions, OutputSink, ReadOptions};
use molframe_core::structure::Structure as CoreStructure;
use std::io::Write;
use std::path::Path;

use crate::structure::Structure;

mod formats;

#[cfg(feature = "chemistry")]
pub use formats::read_component_dictionary;
#[cfg(feature = "geometry")]
pub use formats::transform;
#[cfg(feature = "pdb")]
pub use formats::write_pdb;
#[cfg(feature = "mmcif")]
pub use formats::{read_document, write_mmcif, write_mmcif_with_options};
#[cfg(feature = "bcif")]
pub use formats::{write_bcif, write_bcif_with_options};

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
    write_with_options(path, structure, &WriteOptions::canonical())
}

/// Streams a structure under explicit projection and memory decisions.
///
/// The output policy defaults to 100 MB and no policy above 500 MB is
/// accepted. The file sink itself retains at most 64 KiB; format writers
/// receive rows or columns incrementally and never construct a complete file
/// buffer.
///
/// # Errors
///
/// Returns format, projection, memory-policy, compression, or destination
/// diagnostics.
pub fn write_with_options(
    path: impl AsRef<Path>,
    structure: &Structure,
    options: &WriteOptions,
) -> Result<(), Findings> {
    let path = path.as_ref();
    let name = path.file_name().and_then(|name| name.to_str());
    let Some(format) = name.and_then(Format::from_name) else {
        return Err(Findings::from(
            Diagnostic::new(Code::E1001).with_context("name", path.display().to_string()),
        ));
    };
    let mut output = OutputSink::create(path, options.output).map_err(Findings::from)?;
    write_stream(&mut output, structure, format, options)?;
    output.finish().map_err(Findings::from)
}

/// The decisions a write makes: which identifiers the projection invents, per
/// format family, plus the output working-memory policy.
///
/// [`WriteOptions::canonical`] is the default projection; the `with_*`
/// builders hand the per-family decision types to callers that need explicit
/// identifiers.
#[derive(Clone, Debug)]
pub struct WriteOptions {
    output: OutputOptions,
    #[cfg(feature = "mmcif")]
    cif: molframe_cif::CifWriteOptions,
    #[cfg(feature = "pdb")]
    pdb: molframe_pdb::PdbOptions,
}

impl WriteOptions {
    /// Canonical projection with the default 100 MB output-workspace ceiling.
    #[must_use]
    pub fn canonical() -> Self {
        Self {
            output: OutputOptions::default(),
            #[cfg(feature = "mmcif")]
            cif: molframe_cif::CifWriteOptions::new(),
            #[cfg(feature = "pdb")]
            pdb: molframe_pdb::PdbOptions::new(),
        }
    }

    /// Replaces the output working-memory policy.
    #[must_use]
    pub fn with_output(mut self, output: OutputOptions) -> Self {
        self.output = output;
        self
    }

    /// Replaces the CIF-family projection decisions.
    #[cfg(feature = "mmcif")]
    #[must_use]
    pub fn with_cif(mut self, cif: molframe_cif::CifWriteOptions) -> Self {
        self.cif = cif;
        self
    }

    /// Replaces the PDB-family projection decisions.
    #[cfg(feature = "pdb")]
    #[must_use]
    pub fn with_pdb(mut self, pdb: molframe_pdb::PdbOptions) -> Self {
        self.pdb = pdb;
        self
    }

    /// The output working-memory policy.
    #[must_use]
    pub const fn output(&self) -> &OutputOptions {
        &self.output
    }

    /// The CIF-family projection decisions.
    #[cfg(feature = "mmcif")]
    #[must_use]
    pub const fn cif(&self) -> &molframe_cif::CifWriteOptions {
        &self.cif
    }

    /// The PDB-family projection decisions.
    #[cfg(feature = "pdb")]
    #[must_use]
    pub const fn pdb(&self) -> &molframe_pdb::PdbOptions {
        &self.pdb
    }
}

fn write_stream<W: Write>(
    output: &mut W,
    structure: &Structure,
    format: Format,
    options: &WriteOptions,
) -> Result<(), Findings> {
    #[cfg(not(any(feature = "mmcif", feature = "bcif", feature = "pdb")))]
    let _ = (output, structure, options);
    match format {
        #[cfg(feature = "mmcif")]
        Format::Mmcif => write_mmcif_to_with_options(structure, options.cif(), output)
            .map_err(|error| cif_write_findings(&error)),
        #[cfg(feature = "bcif")]
        Format::BinaryCif => molframe_bcif::write_structure_to_with_memory_limit(
            structure.engine(),
            options.cif(),
            options.output().memory_limit_bytes,
            output,
        )
        .map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pdb => molframe_pdb::write_to(structure.engine(), options.pdb(), output)
            .map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Mmtf => {
            molframe_pdb::write_mmtf_to(structure.engine(), output).map_err(Findings::from)
        }
        #[cfg(feature = "pdb")]
        Format::Pqr => molframe_pdb::write_pqr_to(structure.engine(), options.pdb(), output)
            .map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pdbqt => molframe_pdb::write_pdbqt_to(structure.engine(), options.pdb(), output)
            .map_err(Findings::from),
        other => Err(unsupported_writer(other).into()),
    }
}

/// Reads a structure from an already-acquired buffer, without a further copy.
///
/// This is the zero-copy entry point `read_bytes` builds on: `InputBuffer` is
/// already `Arc`-backed and cheap to clone, so a caller holding one (an FFI
/// boundary, a memory-mapped file) can read through it directly instead of
/// paying `read_bytes`'s `Vec<u8>` copy.
///
/// `name` is only used to settle the format when the content does not.
///
/// # Errors
///
/// Returns the findings that stopped the read.
pub fn read_buffer(
    input: &InputBuffer,
    name: Option<&str>,
    options: &ReadOptions,
) -> Result<(Structure, Vec<Diagnostic>), Findings> {
    dispatch_read(input, name, options)
        .map(|(structure, diagnostics)| (structure.into(), diagnostics))
}

fn dispatch_read(
    input: &InputBuffer,
    name: Option<&str>,
    options: &ReadOptions,
) -> Result<(CoreStructure, Vec<Diagnostic>), Findings> {
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
) -> Result<(CoreStructure, Vec<Diagnostic>), Findings> {
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

#[cfg(test)]
#[path = "facade/facade_tests.rs"]
mod tests;
