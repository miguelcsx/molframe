//! The PDB write options a structure's own metadata implies.
//!
//! A structure read from mmCIF keeps its space group and its
//! non-crystallographic operators as typed extensions that the PDB crate cannot
//! see. Without them `CRYST1` names `P 1` for every crystal and the `MTRIX`
//! records are gone, so the facade completes the caller's options from both.

use molframe_core::structure::Structure;

/// Options extended with the structure's space group and `MTRIX` records.
///
/// A space group or records the caller already chose are never replaced.
pub(crate) fn complete(
    structure: &Structure,
    options: &molframe_pdb::PdbOptions,
) -> molframe_pdb::PdbOptions {
    // Only the crystal metadata below rewrites the options.
    #[cfg_attr(not(feature = "crystal"), allow(unused_mut))]
    let mut completed = options.clone();
    #[cfg(feature = "crystal")]
    {
        let extensions = structure.extensions();
        if structure.data().entry.space_group.is_none()
            && let Some(symbol) = extensions
                .get::<molframe_xtal::SymmetrySet>(molframe_xtal::SYMMETRY_EXTENSION)
                .and_then(|set| set.hermann_mauguin.as_deref())
        {
            completed = completed.with_space_group(symbol);
        }
        if completed.extra_records().is_none()
            && let Some(set) = extensions.get::<molframe_xtal::NcsSet>(molframe_xtal::NCS_EXTENSION)
        {
            let text = mtrix_records(set);
            if !text.is_empty() {
                completed = completed.with_extra_records(text);
            }
        }
    }
    #[cfg(not(feature = "crystal"))]
    let _ = structure;
    completed
}

/// `MTRIX1..3` for the operators whose copies are already in the file.
#[cfg(feature = "crystal")]
fn mtrix_records(set: &molframe_xtal::NcsSet) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    for (position, operator) in set
        .operators()
        .filter(|operator| operator.code == molframe_xtal::NcsCode::Given)
        .enumerate()
    {
        let transform = &operator.transform;
        for row in 0..3 {
            let _ = writeln!(
                out,
                "MTRIX{} {:>3}{:10.6}{:10.6}{:10.6}     {:10.5}    1",
                row + 1,
                position + 1,
                transform.matrix[row][0],
                transform.matrix[row][1],
                transform.matrix[row][2],
                transform.translation[row],
            );
        }
    }
    out
}
