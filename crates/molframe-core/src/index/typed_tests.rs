use super::*;

#[test]
fn positions_of_different_levels_are_distinct_types() {
    // That this file compiles is the assertion: the two share a
    // representation and cannot be substituted for one another.
    let atom = AtomIndex::new(3);
    let residue = ResidueIndex::new(3);
    assert_eq!(atom.get(), residue.get());
}

#[test]
fn debug_names_the_level_and_display_gives_the_bare_ordinal() {
    assert_eq!(format!("{:?}", ChainIndex::new(2)), "chain(2)");
    assert_eq!(format!("{}", ChainIndex::new(2)), "2");
}

#[test]
fn positions_order_by_ordinal_so_tables_can_be_binary_searched() {
    let mut ids = [
        ResidueIndex::new(5),
        ResidueIndex::new(1),
        ResidueIndex::new(3),
    ];
    ids.sort_unstable();
    assert_eq!(
        ids,
        [
            ResidueIndex::new(1),
            ResidueIndex::new(3),
            ResidueIndex::new(5)
        ]
    );
}

#[test]
fn advancing_past_the_last_representable_position_is_rejected() {
    assert_eq!(AtomIndex::new(u32::MAX).next(), None);
    assert_eq!(AtomIndex::new(0).next(), Some(AtomIndex::new(1)));
}
