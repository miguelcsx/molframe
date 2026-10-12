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
    molframe_bcif::write_structure_with_options(structure.engine(), &options)
        .map_err(Findings::from)
}
