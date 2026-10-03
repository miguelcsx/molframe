//! Format detection and the owner each format reads through.

#[cfg(feature = "bcif")]
use super::cif_family::read_bcif_buffer;
#[cfg(feature = "mmcif")]
use super::cif_family::read_mmcif_buffer;
#[cfg(feature = "mmcif")]
use super::extensions;
#[cfg(feature = "pdb")]
use super::pdb_symmetry;
#[cfg(feature = "chemistry")]
use super::small_molecule;
use molframe_core::diagnostic::{Code, Diagnostic, Findings};
use molframe_core::io::{Format, InputBuffer, ReadOptions};
use molframe_core::structure::Structure as CoreStructure;

pub(super) fn dispatch_read(
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
        Format::Pdb => molframe_pdb::read(input, options)
            .map(|(structure, findings)| pdb_symmetry::attach_pdb_symmetry(structure, findings))
            .map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Mmtf => molframe_pdb::read_mmtf(input, options).map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pqr => molframe_pdb::read_pqr(input, options).map_err(Findings::from),
        #[cfg(feature = "pdb")]
        Format::Pdbqt => molframe_pdb::read_pdbqt(input, options).map_err(Findings::from),
        #[cfg(feature = "chemistry")]
        Format::Sdf => small_molecule::read_sdf(input),
        #[cfg(feature = "chemistry")]
        Format::Mol2 => small_molecule::read_mol2(input),
        #[cfg(all(feature = "chemistry", feature = "mmcif"))]
        Format::SmallCif => small_molecule::read_small_cif(input),
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
pub(super) fn unsupported(format: Format) -> Diagnostic {
    Diagnostic::new(Code::E1001)
        .with_message("this build does not include a reader for the detected format")
        .with_context("format", format.name())
        .with_context("crate", crate_for(format))
}

pub(super) fn unsupported_writer(format: Format) -> Diagnostic {
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
        Format::Sdf | Format::Mol2 | Format::SmallCif => "molframe-chem",
        // `Format` is non-exhaustive across crate versions. An unknown variant
        // is unsupported by this compiled facade rather than assigned a guessed
        // owner.
        _ => "unlinked-format",
    }
}
