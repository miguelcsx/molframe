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
        "HELX_P" => match class {
            None | Some(1) => SecondaryStructure::AlphaHelix,
            Some(3) => SecondaryStructure::PiHelix,
            Some(5) => SecondaryStructure::ThreeTenHelix,
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
    if begin_chain != end_chain || begin_sequence > end_sequence {
        return;
    }
    for chain in data.chains() {
        if chain.label() != Some(begin_chain) && chain.auth_label() != Some(begin_chain) {
            continue;
        }
        for residue in chain.residues() {
            let Some(sequence) = residue.label_seq_id().or_else(|| residue.auth_seq_id()) else {
                continue;
            };
            if (begin_sequence..=end_sequence).contains(&sequence)
                && let Some(slot) = states.get_mut(residue.index().as_usize())
            {
                *slot = kind;
            }
        }
    }
}

#[cfg(test)]
#[path = "secondary_tests.rs"]
mod tests;
