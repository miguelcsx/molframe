//! `BinaryCIF` write verbs.

use molframe_core::diagnostic::Findings;

use crate::structure::Structure;

/// Renders deterministic `BinaryCIF` bytes in memory.
///
/// # Errors
///
/// Returns a diagnostic if a projected column cannot be represented.
pub fn write_bcif(structure: &Structure) -> Result<Vec<u8>, Findings> {
    write_bcif_with_options(structure, &molframe_cif::CifWriteOptions::new())
}

/// Renders deterministic `BinaryCIF` in memory with explicit identifier decisions.
///
/// # Errors
///
/// Returns a diagnostic if canonical preflight or binary encoding fails.
pub fn write_bcif_with_options(
    structure: &Structure,
    options: &molframe_cif::CifWriteOptions,
) -> Result<Vec<u8>, Findings> {
    let options = super::complete_cif_options(structure.engine(), options);
    #[cfg(feature = "modelcif")]
    if let Some(model) = structure
        .engine()
        .extensions()
        .get::<molframe_modelcif::ModelCif>(molframe_modelcif::MODEL_CIF_EXTENSION)
    {
        return write_with_model(structure, model, &options);
    }
    molframe_bcif::write_structure_with_options(structure.engine(), &options)
        .map_err(Findings::from)
}

/// Appends the retained `ModelCIF` categories, encoded as binary columns.
#[cfg(feature = "modelcif")]
fn write_with_model(
    structure: &Structure,
    model: &molframe_modelcif::ModelCif,
    options: &molframe_cif::CifWriteOptions,
) -> Result<Vec<u8>, Findings> {
    use molframe_core::diagnostic::{Code, Diagnostic};
    let failure = |message: String| {
        Findings::from(
            Diagnostic::new(Code::E1901)
                .with_message(message)
                .in_category("ma_*"),
        )
    };
    let mut text = b"data_extension\n".to_vec();
    molframe_modelcif::write_extension_to(model, &mut text)
        .map_err(|error| failure(error.to_string()))?;
    let input = molframe_core::io::InputBuffer::from_bytes(text);
    let (document, _) = molframe_cif::parse(&input).map_err(|findings| {
        failure(format!(
            "ModelCIF categories did not re-parse: {findings:?}"
        ))
    })?;
    molframe_bcif::write_structure_with_extension(structure.engine(), options, &document)
        .map_err(Findings::from)
}
