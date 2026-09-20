//! The PDB write verb.

use molframe_core::diagnostic::Findings;

use crate::structure::Structure;

/// Writes a structure in the legacy fixed-column format, or explains why it
/// cannot be written.
///
/// # Errors
///
/// Returns every capacity the structure exceeds.
pub fn write_pdb(
    structure: &Structure,
    options: &molframe_pdb::PdbOptions,
) -> Result<String, Findings> {
    molframe_pdb::write(structure.engine(), options).map_err(Findings::from)
}
