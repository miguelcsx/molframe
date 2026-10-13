use super::super::fixture::sample;
use crate::index::ChainIndex;

#[test]
fn chains_with_their_own_labels_are_left_alone() {
    assert!(sample().with_distinct_chain_labels().is_none());
}

#[test]
fn a_repeated_label_gains_a_number_and_keeps_its_depositor_label() {
    let structure = sample();
    let mut data = structure.data().clone();
    let Ok(repeat) = data.dictionary.intern("A") else {
        panic!("intern");
    };
    assert!(
        data.topology
            .chains
            .relabel(ChainIndex::new(1), repeat)
            .is_ok()
    );
    let repeated = super::Structure::new(data);
    let Some(distinct) = repeated.with_distinct_chain_labels() else {
        panic!("repeat not detected");
    };
    let labels: Vec<_> = distinct
        .data()
        .chains()
        .map(super::super::handle::ChainRef::label)
        .collect();
    assert_eq!(labels, [Some("A"), Some("A2")]);
    let auth: Vec<_> = distinct
        .data()
        .chains()
        .map(super::super::handle::ChainRef::auth_label)
        .collect();
    let original: Vec<_> = repeated
        .data()
        .chains()
        .map(super::super::handle::ChainRef::auth_label)
        .collect();
    assert_eq!(auth, original);
}
