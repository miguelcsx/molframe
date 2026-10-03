//! File-declared per-residue secondary structure from mmCIF categories.

use crate::document::{CifValue, DataBlock};
use molframe_core::diagnostic::Diagnostics;
use molframe_core::structure::StructureData;
use molframe_core::{SecondarySource, SecondaryStructure};

/// Reads helix, sheet and turn ranges in one residue-column allocation.
///
/// Every residue a range names gets its state with source
/// [`SecondarySource::File`]; the rest stay unknown with no source.
pub(super) fn read(
    block: &DataBlock,
    data: &StructureData,
    _findings: &mut Diagnostics,
) -> (Vec<SecondaryStructure>, Vec<SecondarySource>) {
    let mut states = vec![SecondaryStructure::Unknown; data.topology.residues.len()];
    if let Some(category) = block.category("struct_conf") {
        for row in 0..category.row_count() {
            let Some(kind) = conf_state(
                category.text("conf_type_id", row),
                category
                    .value("pdbx_PDB_helix_class", row)
                    .and_then(CifValue::as_integer),
            ) else {
                continue;
            };
            assign_range(
                data,
                &mut states,
                kind,
                category.identifier("beg_label_asym_id", row).as_deref(),
                sequence(category, "beg_label_seq_id", row),
                category.identifier("end_label_asym_id", row).as_deref(),
                sequence(category, "end_label_seq_id", row),
            );
        }
    }
    if let Some(category) = block.category("struct_sheet_range") {
        for row in 0..category.row_count() {
            assign_range(
                data,
                &mut states,
                SecondaryStructure::Strand,
                category.identifier("beg_label_asym_id", row).as_deref(),
                sequence(category, "beg_label_seq_id", row),
                category.identifier("end_label_asym_id", row).as_deref(),
                sequence(category, "end_label_seq_id", row),
            );
        }
    }
    let sources = states
        .iter()
        .map(|state| match state {
            SecondaryStructure::Unknown => SecondarySource::None,
            _ => SecondarySource::File,
        })
        .collect();
    (states, sources)
}

/// The state a `struct_conf` row declares, from its type and helix class.
///
/// The PDB helix class numbers right-handed α as 1, π as 3 and 3₁₀ as 5; a
/// generic `HELX_P` without a class is taken to be α, which is what nearly all
/// deposited generic helices are.
pub(super) fn conf_state(
    conf_type: Option<&str>,
    class: Option<i64>,
) -> Option<SecondaryStructure> {
    let conf_type = conf_type?.to_ascii_uppercase();
    if conf_type.starts_with("TURN") {
        return Some(SecondaryStructure::Turn);
    }
    if !conf_type.starts_with("HELX") {
        return None;
    }
    Some(match conf_type.as_str() {
        "HELX_RH_AL_P" => SecondaryStructure::AlphaHelix,
        "HELX_RH_3T_P" => SecondaryStructure::ThreeTenHelix,
        "HELX_RH_PI_P" => SecondaryStructure::PiHelix,
        "HELX_LH_PP_P" => SecondaryStructure::PolyProline,
        "HELX_P" => match class {
            None | Some(1) => SecondaryStructure::AlphaHelix,
            Some(3) => SecondaryStructure::PiHelix,
            Some(5) => SecondaryStructure::ThreeTenHelix,
            Some(10) => SecondaryStructure::PolyProline,
            Some(_) => SecondaryStructure::OtherHelix,
        },
        _ => SecondaryStructure::OtherHelix,
    })
}

fn sequence(category: &crate::document::Category, item: &str, row: usize) -> Option<i32> {
    category
        .value(item, row)
        .and_then(CifValue::as_integer)
        .and_then(|value| i32::try_from(value).ok())
}

fn assign_range(
    data: &StructureData,
    states: &mut [SecondaryStructure],
    kind: SecondaryStructure,
    begin_chain: Option<&str>,
    begin_sequence: Option<i32>,
    end_chain: Option<&str>,
    end_sequence: Option<i32>,
) {
    let (Some(begin_chain), Some(begin_sequence), Some(end_chain), Some(end_sequence)) =
        (begin_chain, begin_sequence, end_chain, end_sequence)
    else {
        return;
    };
    if begin_chain != end_chain {
        return;
    }
    for chain in data.chains() {
        // Label endpoints never match author identifiers from another namespace.
        if chain.label() != Some(begin_chain) {
            continue;
        }
        let mut first = None;
        let mut last = None;
        let mut ambiguous = false;
        for residue in chain.residues() {
            if residue.label_seq_id() == Some(begin_sequence) {
                ambiguous |= first.replace(residue.index().as_usize()).is_some();
            }
            if residue.label_seq_id() == Some(end_sequence) {
                ambiguous |= last.replace(residue.index().as_usize()).is_some();
            }
        }
        let (Some(first), Some(last)) = (first, last) else {
            continue;
        };
        if !ambiguous
            && first <= last
            && let Some(range) = states.get_mut(first..=last)
        {
            range.fill(kind);
        }
    }
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
