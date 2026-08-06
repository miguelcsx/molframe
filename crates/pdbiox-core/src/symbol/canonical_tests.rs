use super::*;

#[test]
fn no_string_is_declared_twice_because_that_would_orphan_an_identifier() {
    let mut seen = hashbrown::HashSet::with_capacity(CANONICAL.len());
    for text in CANONICAL {
        assert!(seen.insert(*text), "{text} is declared more than once");
    }
}

#[test]
fn every_declared_string_resolves_back_to_its_own_position() {
    for (ordinal, text) in CANONICAL.iter().enumerate() {
        assert_eq!(ordinal_of(text), Some(ordinal as u32), "{text}");
        assert_eq!(text_of(ordinal as u32), Some(*text));
    }
}

#[test]
fn a_string_that_is_not_declared_has_no_canonical_identifier() {
    assert_eq!(ordinal_of("MY_LAB_LIGAND"), None);
    assert_eq!(ordinal_of(""), None);
    assert_eq!(text_of(CANONICAL.len() as u32), None);
}

#[test]
fn the_backbone_atom_names_every_protein_uses_are_canonical() {
    for name in ["N", "CA", "C", "O", "CB", "OXT"] {
        assert!(ordinal_of(name).is_some(), "{name} should be canonical");
    }
}
