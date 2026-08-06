use super::*;
use crate::chunk::{AtomRecord, ChunkBuilder};
use crate::column::Presence;
use crate::element::Element;
use crate::index::{EntityIndex, ResidueIndex};
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::structure::fixture;
use crate::symbol::{AltId, SymbolId};
use crate::topology::{ChainRecord, PolymerKind, ResidueRecord};

fn codes(data: &StructureData) -> Vec<Code> {
    validate(data).iter().map(Diagnostic::code).collect()
}

#[test]
fn a_well_formed_structure_violates_nothing() {
    assert!(validate(fixture::sample().data()).is_empty());
}

#[test]
fn an_empty_structure_violates_nothing() {
    assert!(validate(&StructureData::empty()).is_empty());
}

#[test]
fn a_chain_pointing_at_a_nonexistent_entity_is_reported() {
    let mut data = StructureData::empty();
    data.topology.chains.push(
        ChainRecord {
            label_asym_id: SymbolId::from_raw(0),
            auth_asym_id: OptionalSymbol::NONE,
            entity: EntityIndex::new(7),
            polymer_kind: PolymerKind::None,
        },
        0..0,
    );
    assert!(codes(&data).contains(&Code::E3005));
}

#[test]
fn residues_that_leave_a_gap_in_the_atom_order_are_reported() {
    let mut data = StructureData::empty();
    for range in [0..4u32, 8..12] {
        data.topology.residues.push(
            ResidueRecord {
                label_comp_id: SymbolId::from_raw(0),
                auth_comp_id: OptionalSymbol::NONE,
                label_seq_id: OptionalI32::NONE,
                auth_seq_id: OptionalI32::NONE,
                ins_code: OptionalSymbol::NONE,
                het: false,
            },
            range,
        );
    }
    assert!(codes(&data).contains(&Code::E3004));
}

#[test]
fn models_with_different_atom_counts_cannot_share_a_coordinate_set() {
    let mut data = StructureData::empty();
    let mut builder = ChunkBuilder::new();
    data.topology.residues.push(
        ResidueRecord {
            label_comp_id: SymbolId::from_raw(0),
            auth_comp_id: OptionalSymbol::NONE,
            label_seq_id: OptionalI32::NONE,
            auth_seq_id: OptionalI32::NONE,
            ins_code: OptionalSymbol::NONE,
            het: false,
        },
        0..2,
    );
    for serial in 0..2u32 {
        builder.push(AtomRecord {
            position: Some([0.0; 3]),
            element: Element::CARBON,
            atom_name: SymbolId::from_raw(0),
            auth_atom_name: OptionalSymbol::NONE,
            alt_id: AltId::BLANK,
            residue: ResidueIndex::new(0),
            occupancy: (1.0, Presence::Present),
            b_factor: (0.0, Presence::Present),
            formal_charge: (0, Presence::Inapplicable),
            atom_site_id: serial,
        });
    }
    let (chunks, first) = builder.finish();
    data.chunks = chunks;
    let short = (0..1).map(|_| [0.0; 3]).collect();
    data.coords = CoordinateStore::Dense {
        frames: vec![first, short],
    };

    assert!(codes(&data).contains(&Code::E3010));
}

#[test]
fn an_occupancy_outside_the_permitted_range_is_reported() {
    let mut data = StructureData::empty();
    let mut builder = ChunkBuilder::new();
    data.topology.residues.push(
        ResidueRecord {
            label_comp_id: SymbolId::from_raw(0),
            auth_comp_id: OptionalSymbol::NONE,
            label_seq_id: OptionalI32::NONE,
            auth_seq_id: OptionalI32::NONE,
            ins_code: OptionalSymbol::NONE,
            het: false,
        },
        0..1,
    );
    builder.push(AtomRecord {
        position: Some([0.0; 3]),
        element: Element::CARBON,
        atom_name: SymbolId::from_raw(0),
        auth_atom_name: OptionalSymbol::NONE,
        alt_id: AltId::BLANK,
        residue: ResidueIndex::new(0),
        occupancy: (1.5, Presence::Present),
        b_factor: (0.0, Presence::Present),
        formal_charge: (0, Presence::Inapplicable),
        atom_site_id: 1,
    });
    let (chunks, coords) = builder.finish();
    data.chunks = chunks;
    data.coords = CoordinateStore::Single(coords);

    assert!(codes(&data).contains(&Code::E3009));
}

#[test]
fn every_violation_is_reported_rather_than_only_the_first() {
    let mut data = StructureData::empty();
    for entity in [EntityIndex::new(3), EntityIndex::new(4)] {
        data.topology.chains.push(
            ChainRecord {
                label_asym_id: SymbolId::from_raw(0),
                auth_asym_id: OptionalSymbol::NONE,
                entity,
                polymer_kind: PolymerKind::None,
            },
            0..0,
        );
    }
    assert_eq!(codes(&data).len(), 2);
}
