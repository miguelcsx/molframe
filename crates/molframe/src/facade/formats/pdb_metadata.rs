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

/// `MTRIX1..3` for every operator; the last column marks the copies already
/// in the file.
///
/// Numeric identifiers are kept as the serials, so reading the records back
/// names each operator as before; other identifiers are numbered in order.
#[cfg(feature = "crystal")]
fn mtrix_records(set: &molframe_xtal::NcsSet) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    let numeric: Option<Vec<u32>> = set
        .operators()
        .map(|operator| operator.id.parse().ok())
        .collect();
    let mut ordered: Vec<(u32, &molframe_xtal::NcsOperator)> = set
        .operators()
        .enumerate()
        .map(|(position, operator)| {
            let recorded = numeric.as_ref().and_then(|ids| ids.get(position).copied());
            let id = match recorded {
                Some(id) => id,
                None => molframe_core::numeric::usize_to_u32(position + 1),
            };
            (id, operator)
        })
        .collect();
    ordered.sort_by_key(|(serial, _)| *serial);
    for (serial, operator) in ordered {
        let transform = &operator.transform;
        let given = if operator.code == molframe_xtal::NcsCode::Given {
            '1'
        } else {
            ' '
        };
        for row in 0..3 {
            let _ = writeln!(
                out,
                "MTRIX{} {:>3}{:10.6}{:10.6}{:10.6}     {:10.5}    {}",
                row + 1,
                serial,
                transform.matrix[row][0],
                transform.matrix[row][1],
                transform.matrix[row][2],
                transform.translation[row],
                given,
            );
        }
    }
    out
}
