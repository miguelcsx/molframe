use super::*;
use crate::{ComponentBond, ComponentKind};

fn atom(name: &str, element: Element, charge: i8, aromatic: bool) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge,
        aromatic,
        leaving: false,
        stereo: None,
    }
}

fn bond(first: &str, second: &str, order: BondOrder) -> ComponentBond {
    ComponentBond {
        atom_a: first.into(),
        atom_b: second.into(),
        order,
        aromatic: false,
        stereo: None,
    }
}

/// A glycine-like backbone: amide nitrogen, alpha carbon, carbonyl and oxygen.
fn glycine() -> Arc<Component> {
    Arc::new(Component {
        id: "GLY".into(),
        name: "GLY".into(),
        parent: None,
        one_letter_code: None,
        formula: None,
        kind: ComponentKind::AminoAcid,
        atoms: vec![
            atom("N", Element::NITROGEN, 0, false),
            atom("CA", Element::CARBON, 0, false),
            atom("C", Element::CARBON, 0, false),
            atom("O", Element::OXYGEN, 0, false),
            atom("H", Element::HYDROGEN, 0, false),
        ]
        .into(),
        bonds: vec![
            bond("N", "CA", BondOrder::Single),
            bond("CA", "C", BondOrder::Single),
            bond("C", "O", BondOrder::Double),
            bond("N", "H", BondOrder::Single),
        ]
        .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    })
}

#[test]
fn every_atom_is_found_by_name() {
    let index = ComponentIndex::build(glycine());

    for expected in ["N", "CA", "C", "O", "H"] {
        let found = index.atom(expected).expect("atom is present");
        assert_eq!(found.name.as_ref(), expected);
    }
}

#[test]
fn an_absent_name_resolves_to_nothing_rather_than_a_neighbour() {
    let index = ComponentIndex::build(glycine());

    assert!(index.atom("CB").is_none());
    assert!(index.atom("").is_none());
    assert!(index.position("CA1").is_none());
}

#[test]
fn a_nitrogen_carrying_a_hydrogen_donates() {
    let index = ComponentIndex::build(glycine());

    assert!(index.is_donor("N"));
    assert!(!index.is_donor("CA"), "carbon never donates");
    assert!(!index.is_donor("O"), "this oxygen carries no hydrogen");
}

#[test]
fn a_carbonyl_oxygen_accepts_and_an_amide_nitrogen_does_not() {
    let index = ComponentIndex::build(glycine());

    assert!(index.is_acceptor("O"));
    // N is bonded to CA, which is not doubly bonded to an oxygen, so N is not
    // an amide here and remains an acceptor.
    assert!(index.is_acceptor("N"));
    assert!(!index.is_acceptor("CA"));
}

#[test]
fn a_nitrogen_bonded_to_a_carbonyl_carbon_is_an_amide_and_does_not_accept() {
    let amide = Arc::new(Component {
        id: "AMD".into(),
        name: "AMD".into(),
        parent: None,
        one_letter_code: None,
        formula: None,
        kind: ComponentKind::NonPolymer,
        atoms: vec![
            atom("N", Element::NITROGEN, 0, false),
            atom("C", Element::CARBON, 0, false),
            atom("O", Element::OXYGEN, 0, false),
        ]
        .into(),
        bonds: vec![
            bond("N", "C", BondOrder::Single),
            bond("C", "O", BondOrder::Double),
        ]
        .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    });
    let index = ComponentIndex::build(amide);

    assert!(!index.is_acceptor("N"), "an amide nitrogen does not accept");
    assert!(index.is_acceptor("O"));
}

#[test]
fn a_positively_charged_nitrogen_does_not_accept() {
    let charged = Arc::new(Component {
        id: "CHG".into(),
        name: "CHG".into(),
        parent: None,
        one_letter_code: None,
        formula: None,
        kind: ComponentKind::NonPolymer,
        atoms: vec![atom("NZ", Element::NITROGEN, 1, false)].into(),
        bonds: Vec::new().into(),
        ideal_coordinates: None,
        model_coordinates: None,
    });
    let index = ComponentIndex::build(charged);

    assert!(!index.is_acceptor("NZ"));
}

#[test]
fn a_component_with_no_atoms_indexes_without_panicking() {
    let empty = Arc::new(Component {
        id: "NUL".into(),
        name: "NUL".into(),
        parent: None,
        one_letter_code: None,
        formula: None,
        kind: ComponentKind::NonPolymer,
        atoms: Vec::new().into(),
        bonds: Vec::new().into(),
        ideal_coordinates: None,
        model_coordinates: None,
    });
    let index = ComponentIndex::build(empty);

    assert!(index.atom("N").is_none());
    assert!(!index.is_donor("N"));
    assert!(!index.is_acceptor("N"));
    assert_eq!(index.component().atoms.len(), 0);
}
