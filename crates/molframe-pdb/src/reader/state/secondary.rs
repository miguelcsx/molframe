//! Legacy PDB HELIX and SHEET records lowered to residue state.

use crate::fixed;
use crate::header::PdbHeaders;
use molframe_core::SecondaryStructure;
use molframe_core::structure::StructureData;

pub(super) fn read(data: &StructureData, headers: &PdbHeaders) -> Vec<SecondaryStructure> {
    let mut states = vec![SecondaryStructure::Unknown; data.topology.residues.len()];
    for record in headers.named("HELIX") {
        assign(
            data,
            &mut states,
            SecondaryStructure::Helix,
            fixed::text(record.line(), 20, 20),
            fixed::integer(record.line(), 22, 25).and_then(|value| i32::try_from(value).ok()),
            fixed::text(record.line(), 32, 32),
            fixed::integer(record.line(), 34, 37).and_then(|value| i32::try_from(value).ok()),
        );
    }
    for record in headers.named("SHEET") {
        assign(
            data,
            &mut states,
            SecondaryStructure::Strand,
            fixed::text(record.line(), 22, 22),
            fixed::integer(record.line(), 23, 26).and_then(|value| i32::try_from(value).ok()),
            fixed::text(record.line(), 33, 33),
            fixed::integer(record.line(), 34, 37).and_then(|value| i32::try_from(value).ok()),
        );
    }
    states
}

fn assign(
    data: &StructureData,
    states: &mut [SecondaryStructure],
    kind: SecondaryStructure,
    begin_chain: &str,
    begin_sequence: Option<i32>,
    end_chain: &str,
    end_sequence: Option<i32>,
) {
    let (Some(begin_sequence), Some(end_sequence)) = (begin_sequence, end_sequence) else {
        return;
    };
    if begin_chain.is_empty() || begin_chain != end_chain || begin_sequence > end_sequence {
        return;
    }
    for chain in data.chains() {
        if chain.label() != Some(begin_chain) && chain.auth_label() != Some(begin_chain) {
            continue;
        }
        for residue in chain.residues() {
            let Some(sequence) = residue.auth_seq_id().or_else(|| residue.label_seq_id()) else {
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
