//! Symmetry for legacy PDB files, which name their space group only in `CRYST1`.

use molframe_core::{Diagnostic, Structure};

/// Resolves the space-group name of a PDB `CRYST1` record into the same symmetry
/// extension a CIF read attaches, so a legacy file is as usable crystallographically.
///
/// An unknown or absent symbol leaves the structure unchanged: a name this
/// catalogue does not know is not evidence of any symmetry.
#[cfg(feature = "crystal")]
pub(super) fn attach_pdb_symmetry(
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    let Some(symbol) = structure.data().entry.space_group.clone() else {
        return (structure, findings);
    };
    match molframe_xtal::space_group_by_hermann_mauguin(&symbol) {
        Ok(setting) => (
            structure.with_extension(molframe_xtal::SYMMETRY_EXTENSION, setting.symmetry_set()),
            findings,
        ),
        Err(_) => (structure, findings),
    }
}

#[cfg(not(feature = "crystal"))]
pub(super) fn attach_pdb_symmetry(
    structure: Structure,
    findings: Vec<Diagnostic>,
) -> (Structure, Vec<Diagnostic>) {
    (structure, findings)
}
