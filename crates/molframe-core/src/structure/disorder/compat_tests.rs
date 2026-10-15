use super::*;
use crate::symbol::Interner;

fn label(text: &str) -> AltId {
    let mut interner = Interner::default();
    let Ok(symbol) = interner.intern(text) else {
        panic!("tiny dictionary");
    };
    match AltId::labelled(symbol) {
        Some(label) => label,
        None => panic!("small identifier fits"),
    }
}

#[test]
fn blank_is_compatible_with_everything_and_different_labels_are_not() {
    let (a, b) = (label("A"), label("B"));
    assert!(altloc_compatible(AltId::BLANK, AltId::BLANK));
    assert!(altloc_compatible(AltId::BLANK, a));
    assert!(altloc_compatible(b, AltId::BLANK));
    assert!(altloc_compatible(a, a));
    assert!(!altloc_compatible(a, b));
}

#[test]
fn pairs_skip_incompatible_conformers() {
    let (a, b) = (label("A"), label("B"));
    let labels = [AltId::BLANK, a, b, a];
    let pairs: Vec<_> = compatible_pairs(&labels).collect();
    assert_eq!(pairs, vec![(0, 1), (0, 2), (0, 3), (1, 3)]);
}
