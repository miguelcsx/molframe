use super::*;
use crate::optional::{OptionalI32, OptionalSymbol};
use crate::symbol::SymbolId;
use crate::topology::ResidueRecord;

/// A table of `count` residues, each holding `per` atoms.
fn residues(count: u32, per: u32) -> ResidueTable {
    let mut table = ResidueTable::default();
    for index in 0..count {
        table
            .push(
                ResidueRecord {
                    label_comp_id: SymbolId::from_raw(1),
                    auth_comp_id: OptionalSymbol::NONE,
                    label_seq_id: OptionalI32::some(index.cast_signed()),
                    auth_seq_id: OptionalI32::some(index.cast_signed()),
                    ins_code: OptionalSymbol::NONE,
                    het: false,
                },
                index * per..(index + 1) * per,
            )
            .expect("small residue table");
    }
    table
}

/// The residue of each atom, in order, for the table above.
fn parents_in_order(count: u32, per: u32) -> Vec<u32> {
    (0..count)
        .flat_map(|residue| std::iter::repeat_n(residue, per as usize))
        .collect()
}

#[test]
fn every_mapping_answers_with_the_same_residue_for_every_atom() {
    let (count, per) = (500u32, 9u32);
    let table = residues(count, per);
    let order = parents_in_order(count, per);

    let mappings = [
        ParentMapping::OffsetsOnly,
        ParentMapping::explicit(&order),
        ParentMapping::block_indexed(&order),
    ];
    for mapping in &mappings {
        for atom in 0..count * per {
            let resolved = mapping.resolve(atom, 0, &table);
            assert_eq!(
                resolved,
                Some(ResidueIndex::new(atom / per)),
                "{mapping:?} at atom {atom}"
            );
        }
    }
}

#[test]
fn a_chunk_that_does_not_start_at_the_first_atom_still_resolves() {
    let table = residues(20, 4);
    let order = parents_in_order(20, 4);
    // A chunk covering atoms 40 onwards, so residue 10 onwards.
    let mapping = ParentMapping::block_indexed(&order[40..]);
    assert_eq!(mapping.resolve(0, 40, &table), Some(ResidueIndex::new(10)));
    assert_eq!(mapping.resolve(7, 40, &table), Some(ResidueIndex::new(11)));
}

#[test]
fn the_block_mapping_costs_a_thirty_second_of_the_explicit_one() {
    let order = parents_in_order(100, 10);
    let explicit = ParentMapping::explicit(&order);
    let blocked = ParentMapping::block_indexed(&order);
    assert!(blocked.bytes() * 32 <= explicit.bytes());
    assert_eq!(ParentMapping::OffsetsOnly.bytes(), 0);
}

#[test]
fn an_atom_past_the_last_residue_resolves_to_nothing_rather_than_looping() {
    let table = residues(3, 2);
    let order = parents_in_order(3, 2);
    for mapping in [
        ParentMapping::OffsetsOnly,
        ParentMapping::block_indexed(&order),
    ] {
        assert_eq!(mapping.resolve(99, 0, &table), None);
    }
}
