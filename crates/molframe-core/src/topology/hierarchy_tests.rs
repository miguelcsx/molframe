use super::*;
use crate::index::{ChainIndex, EntityIndex, ModelIndex, ResidueIndex};
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::symbol::SymbolId;
use crate::topology::{ChainRecord, ChainTable, EntityKind, PolymerKind, ResidueRecord};

fn sample() -> Topology {
    let mut topology = Topology::default();
    let entity = topology
        .entities
        .push(
            SymbolId::from_raw(1),
            EntityKind::Polymer,
            OptionalSymbol::NONE,
            &[
                SymbolId::from_raw(10),
                SymbolId::from_raw(11),
                SymbolId::from_raw(12),
            ],
        )
        .expect("small entity table");
    for start in [0u32, 4] {
        topology
            .residues
            .push(
                ResidueRecord {
                    label_comp_id: SymbolId::from_raw(10),
                    auth_comp_id: OptionalSymbol::NONE,
                    label_seq_id: OptionalI32::some(start.cast_signed() / 4 + 1),
                    auth_seq_id: OptionalI32::some(start.cast_signed() / 4 + 1),
                    ins_code: OptionalSymbol::NONE,
                    het: false,
                },
                start..start + 4,
            )
            .expect("small residue table");
    }
    topology
        .chains
        .push(
            ChainRecord {
                label_asym_id: SymbolId::from_raw(2),
                auth_asym_id: OptionalSymbol::NONE,
                entity,
                polymer_kind: PolymerKind::Protein,
            },
            0..2,
        )
        .expect("small chain table");
    topology.models.push(1, 0..1).expect("small model table");
    topology
}

#[test]
fn a_residues_atoms_are_a_range_computed_from_two_loads() {
    let topology = sample();
    assert_eq!(topology.residues.atoms(ResidueIndex::new(0)), Some(0..4));
    assert_eq!(topology.residues.atoms(ResidueIndex::new(1)), Some(4..8));
    assert_eq!(topology.residues.atoms(ResidueIndex::new(2)), None);
}

#[test]
fn walking_down_the_hierarchy_yields_ranges_at_every_level() {
    let topology = sample();
    let chains = topology.models.chains(ModelIndex::new(0));
    assert_eq!(chains, Some(0..1));
    assert_eq!(topology.chains.residues(ChainIndex::new(0)), Some(0..2));
    assert_eq!(topology.atom_count(), 8);
}

#[test]
fn a_deposited_model_number_is_kept_rather_than_replaced_by_its_position() {
    let mut topology = Topology::default();
    topology.models.push(5, 0..0).expect("small model table");
    topology.models.push(7, 0..0).expect("small model table");
    assert_eq!(topology.models.model_num(ModelIndex::new(0)), Some(5));
    assert_eq!(topology.models.model_num(ModelIndex::new(1)), Some(7));
}

#[test]
fn the_residue_containing_an_atom_is_found_by_searching_the_offsets() {
    let topology = sample();
    assert_eq!(topology.residues.containing(0), Some(ResidueIndex::new(0)));
    assert_eq!(topology.residues.containing(3), Some(ResidueIndex::new(0)));
    assert_eq!(topology.residues.containing(4), Some(ResidueIndex::new(1)));
    assert_eq!(topology.residues.containing(7), Some(ResidueIndex::new(1)));
    assert_eq!(topology.residues.containing(8), None);
}

#[test]
fn chains_that_are_copies_of_one_species_are_found_by_reading_a_column() {
    let mut topology = Topology::default();
    let first = topology
        .entities
        .push(
            SymbolId::from_raw(1),
            EntityKind::Polymer,
            OptionalSymbol::NONE,
            &[],
        )
        .expect("small entity table");
    let second = topology
        .entities
        .push(
            SymbolId::from_raw(2),
            EntityKind::Water,
            OptionalSymbol::NONE,
            &[],
        )
        .expect("small entity table");
    for entity in [first, second, first] {
        topology
            .chains
            .push(
                ChainRecord {
                    label_asym_id: SymbolId::from_raw(3),
                    auth_asym_id: OptionalSymbol::NONE,
                    entity,
                    polymer_kind: PolymerKind::None,
                },
                0..0,
            )
            .expect("small chain table");
    }
    let copies: Vec<_> = topology.chains.instances_of(first).collect();
    assert_eq!(copies, [ChainIndex::new(0), ChainIndex::new(2)]);
}

#[test]
fn an_entitys_canonical_sequence_is_what_should_be_there_not_what_was_modelled() {
    let topology = sample();
    let sequence = topology.entities.canonical_sequence(EntityIndex::new(0));
    assert_eq!(sequence.len(), 3, "three residues in the sequence");
    assert_eq!(topology.chains.residues(ChainIndex::new(0)), Some(0..2));
}

#[test]
fn an_entity_is_found_by_the_identifier_the_file_declared() {
    let topology = sample();
    assert_eq!(
        topology.entities.find_by_id(SymbolId::from_raw(1)),
        Some(EntityIndex::new(0))
    );
    assert_eq!(topology.entities.find_by_id(SymbolId::from_raw(99)), None);
}

#[test]
fn an_absent_depositor_label_is_distinguishable_from_one_equal_to_the_label() {
    let mut chains = ChainTable::default();
    chains
        .push(
            ChainRecord {
                label_asym_id: SymbolId::from_raw(1),
                auth_asym_id: OptionalSymbol::NONE,
                entity: EntityIndex::new(0),
                polymer_kind: PolymerKind::None,
            },
            0..0,
        )
        .expect("small chain table");
    chains
        .push(
            ChainRecord {
                label_asym_id: SymbolId::from_raw(1),
                auth_asym_id: OptionalSymbol::some(SymbolId::from_raw(1)),
                entity: EntityIndex::new(0),
                polymer_kind: PolymerKind::None,
            },
            0..0,
        )
        .expect("small chain table");
    assert_eq!(chains.auth_asym_id(ChainIndex::new(0)), None);
    assert_eq!(
        chains.auth_asym_id(ChainIndex::new(1)),
        Some(SymbolId::from_raw(1))
    );
}
