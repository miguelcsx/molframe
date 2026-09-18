use super::*;
use crate::structure::fixture;

#[test]
fn walking_down_the_hierarchy_reaches_a_named_atom() {
    let structure = fixture::sample();
    let data = structure.data();

    let atom = data
        .model(ModelIndex::new(0))
        .and_then(|model| model.chain("A"))
        .and_then(|chain| chain.residue(100))
        .and_then(|residue| residue.atom("CA"));

    assert_eq!(atom.and_then(AtomRef::name), Some("CA"));
    assert_eq!(atom.and_then(AtomRef::element), Some(Element::CARBON));
}

#[test]
fn a_chain_is_addressable_by_either_namespace() {
    let structure = fixture::sample();
    let Some(model) = structure.data().model(ModelIndex::new(0)) else {
        panic!("expected a model")
    };
    assert!(model.chain("A").is_some());
    assert!(model.chain("B").is_some());
    assert!(model.chain("Z").is_none());
}

#[test]
fn iterating_every_level_visits_each_member_once() {
    let structure = fixture::sample();
    let data = structure.data();

    assert_eq!(data.models().count(), 1);
    assert_eq!(data.chains().count(), 2);
    assert_eq!(data.residues().count(), 6);
    assert_eq!(data.atoms().count(), 24);

    let through_hierarchy: usize = data
        .models()
        .flat_map(ModelRef::chains)
        .flat_map(ChainRef::residues)
        .flat_map(ResidueRef::atoms)
        .count();
    assert_eq!(through_hierarchy, 24);
}

#[test]
fn an_atom_reports_the_residue_it_belongs_to() {
    let structure = fixture::sample();
    let data = structure.data();

    for atom in data.atoms() {
        let residue = atom.residue().map(ResidueRef::index);
        assert_eq!(residue, Some(ResidueIndex::new(atom.index().get() / 4)));
    }
}

#[test]
fn a_deposited_residue_number_survives_rather_than_being_replaced_by_a_position() {
    let structure = fixture::sample();
    let numbers: Vec<_> = structure
        .data()
        .chains()
        .flat_map(ChainRef::residues)
        .filter_map(ResidueRef::auth_seq_id)
        .collect();
    assert_eq!(numbers, [100, 101, 102, 100, 101, 102]);
}

#[test]
fn a_residue_carries_both_namespaces_without_deriving_one_from_the_other() {
    let structure = fixture::sample();
    let Some(residue) = structure.data().residues().next() else {
        panic!("expected a residue")
    };
    assert_eq!(residue.label_seq_id(), Some(1));
    assert_eq!(residue.auth_seq_id(), Some(100));
    assert_eq!(residue.name(), Some("ALA"));
    assert_eq!(residue.ins_code(), None);
    assert!(!residue.is_het());
}

#[test]
fn every_atom_reports_its_recorded_position() {
    let structure = fixture::sample();
    let missing = structure
        .data()
        .atoms()
        .filter(|atom| atom.position().is_none())
        .count();
    assert_eq!(missing, 0);
}
