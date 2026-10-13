//! The PQR and PDBQT write verbs that report what they cannot carry.

use molframe_core::diagnostic::{Diagnostic, Findings};

use crate::structure::Structure;

/// Writes PQR and reports what the format could not carry.
///
/// # Errors
///
/// Returns the capacity or missing-annotation findings that stop the write.
pub fn write_pqr_with_findings(
    structure: &Structure,
    options: &molframe_pdb::PdbOptions,
) -> Result<(String, Findings), Findings> {
    let text = molframe_pdb::write_pqr(structure.engine(), options).map_err(Findings::from)?;
    let losses = molframe_pdb::pqr_loss(structure.engine(), options);
    Ok((text, with_extensions(structure, losses, "PQR")))
}

/// Writes a rigid PDBQT model and reports what the format could not carry.
///
/// # Errors
///
/// Returns the capacity or missing-annotation findings that stop the write.
pub fn write_pdbqt_with_findings(
    structure: &Structure,
    options: &molframe_pdb::PdbOptions,
) -> Result<(String, Findings), Findings> {
    let text = molframe_pdb::write_pdbqt(structure.engine(), options).map_err(Findings::from)?;
    let losses = molframe_pdb::pdbqt_loss(structure.engine(), options);
    Ok((text, with_extensions(structure, losses, "PDBQT")))
}

fn with_extensions(structure: &Structure, mut losses: Vec<Diagnostic>, target: &str) -> Findings {
    losses.extend(super::loss::uncarried_extensions(
        structure.engine(),
        &[],
        target,
    ));
    Findings::from(losses)
}
