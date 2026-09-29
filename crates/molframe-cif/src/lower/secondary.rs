//! File-declared per-residue secondary structure from mmCIF categories.

use crate::document::{CifValue, DataBlock};
use molframe_core::SecondaryStructure;
use molframe_core::diagnostic::Diagnostics;
use molframe_core::structure::StructureData;

/// Reads helix, sheet and turn ranges in one residue-column allocation.
pub(super) fn read(
    block: &DataBlock,
    data: &StructureData,
    _findings: &mut Diagnostics,
) -> Vec<SecondaryStructure> {
    let mut states = vec![SecondaryStructure::Unknown; data.topology.residues.len()];
    if let Some(category) = block.category("struct_conf") {
        for row in 0..category.row_count() {
            let kind = match category.text("conf_type_id", row) {
                Some(value) if value.to_ascii_uppercase().contains("HELX") => {
                    SecondaryStructure::Helix
                }
                Some(value) if value.to_ascii_uppercase().contains("TURN") => {
                    SecondaryStructure::Turn
                }
                _ => continue,
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
    states
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
