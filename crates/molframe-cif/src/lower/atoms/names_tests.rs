use super::*;

#[test]
fn a_residue_name_set_answers_the_same_way_either_side_of_its_threshold() {
    let mut names = ResidueNames::default();
    let alt = AltId::BLANK;

    for raw in 0..200_u32 {
        let name = SymbolId::from_raw(raw);
        assert!(!names.contains(name, alt), "reported a repeat for {raw}");
        names.insert(name, alt);
        assert!(names.contains(name, alt), "lost {raw} after inserting it");
    }

    // Everything inserted is still present once the index has taken over.
    for raw in 0..200_u32 {
        assert!(names.contains(SymbolId::from_raw(raw), alt), "lost {raw}");
    }
    assert!(!names.contains(SymbolId::from_raw(200), alt));
}

#[test]
fn clearing_a_residue_name_set_drops_its_index_too() {
    let mut names = ResidueNames::default();
    let alt = AltId::BLANK;
    for raw in 0..200_u32 {
        names.insert(SymbolId::from_raw(raw), alt);
    }
    names.clear();
    assert!(!names.contains(SymbolId::from_raw(0), alt));
    assert!(!names.contains(SymbolId::from_raw(199), alt));
}

#[test]
fn the_alternate_location_is_part_of_the_identity() {
    let mut names = ResidueNames::default();
    let name = SymbolId::from_raw(7);
    names.insert(name, AltId::BLANK);
    assert!(names.contains(name, AltId::BLANK));
    let Some(labelled) = AltId::labelled(SymbolId::from_raw(3)) else {
        panic!("a labelled alternate location must be constructible")
    };
    assert!(!names.contains(name, labelled));
}
