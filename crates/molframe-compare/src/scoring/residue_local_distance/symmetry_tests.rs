use super::SymmetryCache;
use molframe_chem::{Component, ComponentAtom, ComponentBond, ComponentKind};
use molframe_core::{BondOrder, Element};

fn atom(name: &str, element: Element) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
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

fn fork(order_one: BondOrder, order_two: BondOrder, leaf: Element) -> Component {
    Component {
        id: "FRK".into(),
        name: "fork".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: vec![
            atom("CB", Element::CARBON),
            atom("A", leaf),
            atom("B", leaf),
        ]
        .into(),
        bonds: vec![bond("CB", "A", order_one), bond("CB", "B", order_two)].into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

#[test]
fn leaves_on_a_planar_centre_are_exchanged_although_only_one_bond_is_double() {
    let component = fork(BondOrder::Double, BondOrder::Single, Element::OXYGEN);
    let found = SymmetryCache::default().permutations(&component, &[0, 1, 2], 16);
    assert_eq!(found.map(|all| all.len()), Some(1));
}

#[test]
fn leaves_on_a_tetrahedral_centre_are_not_exchanged() {
    let component = fork(BondOrder::Single, BondOrder::Single, Element::CARBON);
    let found = SymmetryCache::default().permutations(&component, &[0, 1, 2], 16);
    assert_eq!(found.map(|all| all.len()), Some(0));
}

#[test]
fn exceeding_the_automorphism_limit_reports_no_result() {
    let component = fork(BondOrder::Double, BondOrder::Single, Element::OXYGEN);
    assert!(
        SymmetryCache::default()
            .permutations(&component, &[0, 1, 2], 1)
            .is_none()
    );
}
