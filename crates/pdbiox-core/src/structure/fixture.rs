//! A small structure the tests in this module share.
//!
//! Two chains of one entity, three residues each, four atoms each — enough to
//! exercise every level of the hierarchy without any of the arithmetic being
//! coincidentally right.

use super::data::{CoordinateStore, Structure, StructureData};
use crate::chunk::{AtomRecord, ChunkBuilder};
use crate::column::Presence;
use crate::element::Element;
use crate::index::ResidueIndex;
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::symbol::AltId;
use crate::topology::{ChainRecord, EntityKind, PolymerKind, ResidueRecord};

/// Atom names cycled through each residue.
const ATOM_NAMES: [&str; 4] = ["N", "CA", "C", "O"];

/// Builds the shared fixture.
///
/// Panics only if the identifier dictionary refuses a handful of strings, which
/// cannot happen at this size and would be a defect in the interner.
pub fn sample() -> Structure {
    let mut data = StructureData::empty();
    let dictionary = &mut data.dictionary;

    let Ok(entity_id) = dictionary.intern("1") else {
        return Structure::new(data);
    };
    let Ok(alanine) = dictionary.intern("ALA") else {
        return Structure::new(data);
    };
    let sequence = [alanine; 4];
    let entity = data.topology.entities.push(
        entity_id,
        EntityKind::Polymer,
        OptionalSymbol::NONE,
        &sequence,
    );

    let mut builder = ChunkBuilder::new();
    builder.start_model(0);
    let mut residue_position = 0u32;

    for (chain_number, label) in ["A", "B"].into_iter().enumerate() {
        let Ok(chain_label) = data.dictionary.intern(label) else {
            continue;
        };
        let first_residue = residue_position;

        for offset in 0..3u32 {
            let residue = data.topology.residues.push(
                ResidueRecord {
                    label_comp_id: alanine,
                    auth_comp_id: OptionalSymbol::NONE,
                    label_seq_id: OptionalI32::some(offset.cast_signed() + 1),
                    auth_seq_id: OptionalI32::some(offset.cast_signed() + 100),
                    ins_code: OptionalSymbol::NONE,
                    het: false,
                },
                residue_position * 4..(residue_position + 1) * 4,
            );

            for (slot, name) in ATOM_NAMES.iter().enumerate() {
                let Ok(atom_name) = data.dictionary.intern(name) else {
                    continue;
                };
                builder.push(AtomRecord {
                    position: Some([chain_number as f32, offset as f32, slot as f32]),
                    element: element_of(slot),
                    atom_name,
                    auth_atom_name: OptionalSymbol::NONE,
                    alt_id: AltId::BLANK,
                    residue: ResidueIndex::new(residue.get()),
                    occupancy: (1.0, Presence::Present),
                    b_factor: (20.0 + offset as f32, Presence::Present),
                    formal_charge: (0, Presence::Inapplicable),
                    atom_site_id: residue_position * 4 + slot as u32 + 1,
                });
            }
            residue_position += 1;
        }

        data.topology.chains.push(
            ChainRecord {
                label_asym_id: chain_label,
                auth_asym_id: OptionalSymbol::some(chain_label),
                entity,
                polymer_kind: PolymerKind::Protein,
            },
            first_residue..residue_position,
        );
    }

    data.topology.models.push(1, 0..2);
    let (chunks, coords) = builder.finish();
    data.chunks = chunks;
    data.coords = CoordinateStore::Single(coords);
    Structure::new(data)
}

fn element_of(slot: usize) -> Element {
    match slot {
        0 => Element::NITROGEN,
        3 => Element::OXYGEN,
        _ => Element::CARBON,
    }
}
