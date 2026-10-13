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
    molframe_pdb::write(
        structure.engine(),
        &super::complete_pdb_options(structure.engine(), options),
    )
    .map_err(Findings::from)
}

/// Writes a structure like [`write_pdb`] and also reports what the format
/// could not carry.
///
/// Each loss is one aggregated `W4101` finding with a count and an example, so
/// a lossy projection from mmCIF is never silent.
///
/// # Errors
///
/// Returns every capacity the structure exceeds.
pub fn write_pdb_with_findings(
    structure: &Structure,
    options: &molframe_pdb::PdbOptions,
) -> Result<(String, Findings), Findings> {
    let completed = super::complete_pdb_options(structure.engine(), options);
    let text = molframe_pdb::write(structure.engine(), &completed).map_err(Findings::from)?;
    let losses = molframe_pdb::projection_loss(structure.engine(), &completed);
    Ok((text, Findings::from(losses)))
}
