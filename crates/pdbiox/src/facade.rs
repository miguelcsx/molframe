//! Reading and writing, dispatched to whichever format crate is linked.
//!
//! The dispatch is the only thing here. Deciding which format a file is belongs
//! to the core, reading it belongs to the format crate, and this is the seam
//! between them.

use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::{Format, InputBuffer, Limits, ReadOptions, ReadResult, write_output};
use pdbiox_core::structure::Structure;
use std::path::Path;

#[cfg(feature = "chem")]
use pdbiox_core::contract::DictionaryVersion;

#[cfg(feature = "geom")]
use pdbiox_core::selection::AtomSelection;
#[cfg(feature = "geom")]
use pdbiox_geom::Rigid;

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
pub fn write(path: impl AsRef<Path>, structure: &Structure) -> Result<(), Vec<Diagnostic>> {
    let path = path.as_ref();
    let name = path.file_name().and_then(|name| name.to_str());
    let Some(format) = name.and_then(Format::from_name) else {
        return Err(vec![
            Diagnostic::new(Code::E1001).with_context("name", path.display().to_string()),
        ]);
    };
    let rendered = render(structure, format)?;
    write_output(path, &rendered).map_err(|finding| vec![finding])
}

fn render(structure: &Structure, format: Format) -> Result<Vec<u8>, Vec<Diagnostic>> {
    match format {
        #[cfg(feature = "mmcif")]
        Format::Mmcif => {
            #[cfg(feature = "modelcif")]
            if let Some(model) = structure
                .extensions()
                .get::<pdbiox_modelcif::ModelCif>(pdbiox_modelcif::MODEL_CIF_EXTENSION)
            {
                return pdbiox_modelcif::write_canonical(structure, model)
                    .map(String::into_bytes)
                    .map_err(|error| cif_write_findings(&error));
            }
            pdbiox_cif::write_canonical(structure)
                .map(String::into_bytes)
                .map_err(|error| cif_write_findings(&error))
        }
        #[cfg(feature = "mmcif")]
        Format::Pdbml => write_pdbml_structure(structure),
        #[cfg(feature = "bcif")]
        Format::BinaryCif => pdbiox_bcif::write_structure(structure),
        #[cfg(feature = "pdb")]
        Format::Pdb => {
            pdbiox_pdb::write(structure, &pdbiox_pdb::PdbOptions::new()).map(String::into_bytes)
        }
        #[cfg(feature = "pdb")]
        Format::Mmtf => pdbiox_pdb::write_mmtf(structure),
        #[cfg(feature = "pdb")]
        Format::Pqr => {
            pdbiox_pdb::write_pqr(structure, &pdbiox_pdb::PdbOptions::new()).map(String::into_bytes)
        }
        #[cfg(feature = "pdb")]
        Format::Pdbqt => pdbiox_pdb::write_pdbqt(structure, &pdbiox_pdb::PdbOptions::new())
            .map(String::into_bytes),
        other => Err(vec![unsupported(other)]),
    }
}

fn read_buffer(input: &InputBuffer, name: Option<&str>, options: &ReadOptions) -> ReadResult {
    let format = Format::detect(options.format, input, name).map_err(|finding| vec![finding])?;
    match format {
        #[cfg(feature = "mmcif")]
        Format::Mmcif => read_mmcif_buffer(input, options),
        #[cfg(feature = "mmcif")]
        Format::Pdbml => read_pdbml_buffer(input, options),
        #[cfg(feature = "bcif")]
        Format::BinaryCif => read_bcif_buffer(input, options),
        #[cfg(feature = "pdb")]
        Format::Pdb => pdbiox_pdb::read(input, options),
        #[cfg(feature = "pdb")]
        Format::Mmtf => pdbiox_pdb::read_mmtf(input, options),
        #[cfg(feature = "pdb")]
        Format::Pqr => pdbiox_pdb::read_pqr(input, options),
        #[cfg(feature = "pdb")]
        Format::Pdbqt => pdbiox_pdb::read_pdbqt(input, options),
        other => Err(vec![unsupported(other)]),
    }
}

#[cfg(feature = "mmcif")]
fn read_pdbml_buffer(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    match pdbiox_cif::read_pdbml(input.as_bytes(), options) {
        Ok((document, structure, findings)) => {
            attach_cif_metadata(&document, structure, findings, options)
        }
        Err(pdbiox_cif::PdbmlReadError::Findings(findings)) => Err(findings),
        Err(pdbiox_cif::PdbmlReadError::Pdbml(error)) => Err(vec![
            Diagnostic::new(Code::E1102)
                .with_message("PDBML/XML could not be decoded")
                .with_context("decoder", error.to_string()),
        ]),
        Err(error) => Err(vec![
            Diagnostic::new(Code::E1102)
                .with_message("PDBML/XML could not be decoded")
                .with_context("decoder", error.to_string()),
        ]),
    }
}

#[cfg(feature = "mmcif")]
fn write_pdbml_structure(structure: &Structure) -> Result<Vec<u8>, Vec<Diagnostic>> {
    let canonical =
        pdbiox_cif::write_canonical(structure).map_err(|error| cif_write_findings(&error))?;
    let input = InputBuffer::from_bytes(canonical.into_bytes());
    let (document, _) = pdbiox_cif::parse(&input)?;
    pdbiox_cif::write_pdbml(&document)
        .map(String::into_bytes)
        .map_err(|error| {
            vec![
                Diagnostic::new(Code::E1102)
                    .with_message("PDBML/XML could not be encoded")
                    .with_context("encoder", error.to_string()),
            ]
        })
}

#[cfg(feature = "bcif")]
fn read_bcif_buffer(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    let (document, structure, findings) = pdbiox_bcif::read_with_document(input, options)?;
    if options.only_atomic_coords {
        return Ok((structure, findings));
    }
    attach_cif_metadata(&document, structure, findings, options)
}

#[cfg(feature = "mmcif")]
fn read_mmcif_buffer(input: &InputBuffer, options: &ReadOptions) -> ReadResult {
    let (document, structure, findings) = pdbiox_cif::read_with_document(input, options)?;
    if options.only_atomic_coords {
        return Ok((structure, findings));
    }
    attach_cif_metadata(&document, structure, findings, options)
}

#[cfg(any(feature = "mmcif", feature = "bcif"))]
fn attach_cif_metadata(
    document: &pdbiox_cif::Document,
    structure: Structure,
    findings: Vec<Diagnostic>,
    options: &ReadOptions,
) -> ReadResult {
    let (structure, findings) = attach_optional_extensions(document, structure, findings);
    options.finish(structure, findings)
}

#[cfg(any(feature = "modelcif", feature = "xtal"))]
fn attach_optional_extensions(
    document: &pdbiox_cif::Document,
    mut structure: Structure,
    mut findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    #[cfg(feature = "modelcif")]
    {
        let (model, model_findings) = pdbiox_modelcif::lower(document);
        findings.extend(model_findings);
        if !model.is_empty() {
            structure = structure.with_extension(pdbiox_modelcif::MODEL_CIF_EXTENSION, model);
        }
    }
    #[cfg(feature = "xtal")]
    {
        match pdbiox_xtal::lower_assemblies(document) {
            Ok(assemblies) if !assemblies.is_empty() => {
                structure = structure.with_extension(pdbiox_xtal::ASSEMBLIES_EXTENSION, assemblies);
            }
            Ok(_) => {}
            Err(assembly_findings) => findings.extend(assembly_findings),
        }
        match pdbiox_xtal::lower_ncs(document) {
            Ok(ncs) if !ncs.is_empty() => {
                structure = structure.with_extension(pdbiox_xtal::NCS_EXTENSION, ncs);
            }
            Ok(_) => {}
            Err(ncs_findings) => findings.extend(ncs_findings),
        }
        match pdbiox_xtal::lower_symmetry(document) {
            Ok(symmetry) if !symmetry.is_empty() => {
                structure = structure.with_extension(pdbiox_xtal::SYMMETRY_EXTENSION, symmetry);
            }
            Ok(_) => {}
            Err(symmetry_findings) => findings.extend(symmetry_findings),
        }
    }
    (structure, findings)
}

#[cfg(not(any(feature = "modelcif", feature = "xtal")))]
fn attach_optional_extensions(
    _document: &pdbiox_cif::Document,
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    (structure, findings)
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
        Format::Mmcif | Format::Pdbml => "pdbiox-cif",
        Format::BinaryCif => "pdbiox-bcif",
        Format::Pdb | Format::Pqr | Format::Pdbqt | Format::Mmtf => "pdbiox-pdb",
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
pub fn read_document(path: impl AsRef<Path>) -> Result<pdbiox_cif::Document, Vec<Diagnostic>> {
    let path = path.as_ref();
    let input = InputBuffer::open(path, Limits::default()).map_err(|finding| vec![finding])?;
    let name = path.file_name().and_then(|name| name.to_str());
    let format = Format::detect(Format::Auto, &input, name).map_err(|finding| vec![finding])?;
    match format {
        Format::Mmcif => pdbiox_cif::parse(&input).map(|(document, _)| document),
        Format::Pdbml => pdbiox_cif::parse_pdbml_document(input.as_bytes())
            .map_err(|error| vec![Diagnostic::new(Code::E1102).with_message(error.to_string())]),
        #[cfg(feature = "bcif")]
        Format::BinaryCif => pdbiox_bcif::read_document(&input, Limits::default())
            .and_then(|document| document.to_document())
            .map_err(|finding| vec![finding]),
        other => Err(vec![unsupported(other)]),
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
) -> Result<(pdbiox_chem::CifProvider, Vec<Diagnostic>), Vec<Diagnostic>> {
    let input = InputBuffer::open(path, Limits::default()).map_err(|finding| vec![finding])?;
    pdbiox_chem::read_ccd(&input, version)
}

/// Writes a structure as a valid, self-consistent mmCIF file.
///
/// # Errors
///
/// Returns the canonical preflight error without returning partial output.
#[cfg(feature = "mmcif")]
pub fn write_mmcif(structure: &Structure) -> Result<String, pdbiox_cif::CifWriteError> {
    write_mmcif_with_options(structure, &pdbiox_cif::CifWriteOptions::new())
}

/// Writes a structure as mmCIF with explicit identifier decisions.
///
/// # Errors
///
/// Returns the canonical preflight error without returning partial output.
#[cfg(feature = "mmcif")]
pub fn write_mmcif_with_options(
    structure: &Structure,
    options: &pdbiox_cif::CifWriteOptions,
) -> Result<String, pdbiox_cif::CifWriteError> {
    #[cfg(feature = "modelcif")]
    if let Some(model) = structure
        .extensions()
        .get::<pdbiox_modelcif::ModelCif>(pdbiox_modelcif::MODEL_CIF_EXTENSION)
    {
        return pdbiox_modelcif::write_canonical_with_options(structure, model, options);
    }
    pdbiox_cif::write_canonical_with_options(structure, options)
}

/// Writes a structure as deterministic `BinaryCIF` bytes.
///
/// # Errors
///
/// Returns a diagnostic if a projected column cannot be represented.
#[cfg(feature = "bcif")]
pub fn write_bcif(structure: &Structure) -> Result<Vec<u8>, Vec<Diagnostic>> {
    pdbiox_bcif::write_structure(structure)
}

/// Writes deterministic `BinaryCIF` with explicit canonical identifier decisions.
///
/// # Errors
///
/// Returns a diagnostic if canonical preflight or binary encoding fails.
#[cfg(feature = "bcif")]
pub fn write_bcif_with_options(
    structure: &Structure,
    options: &pdbiox_cif::CifWriteOptions,
) -> Result<Vec<u8>, Vec<Diagnostic>> {
    pdbiox_bcif::write_structure_with_options(structure, options)
}

/// The ceilings a read runs under by default.
#[must_use]
pub fn default_limits() -> Limits {
    Limits::default()
}

#[cfg(feature = "mmcif")]
fn cif_write_findings(error: &pdbiox_cif::CifWriteError) -> Vec<Diagnostic> {
    vec![
        Diagnostic::new(Code::E4105)
            .with_message("structure cannot be projected to canonical CIF")
            .with_context("reason", error.to_string()),
    ]
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
) -> Result<Structure, Vec<Diagnostic>> {
    if structure.ragged_models().is_some() {
        return Err(vec![Diagnostic::new(Code::E6008)]);
    }
    if let Some(atom) = selection
        .iter()
        .find(|atom| *atom >= structure.atom_count())
    {
        return Err(vec![
            Diagnostic::new(Code::E6009).with_context("atom", atom.to_string()),
        ]);
    }

    let mut editor = structure.edit();
    editor
        .transform(selection, |position| rigid.apply(position))
        .map_err(|finding| vec![finding])?;
    editor.commit()
}

#[cfg(test)]
#[path = "facade_tests.rs"]
mod tests;
