use super::side_chain_definition;
use crate::{Component, ComponentAtom, ComponentBond, ComponentKind};
use pdbiox_core::{BondOrder, Element};
use std::sync::Arc;

#[test]
fn arginine_path_reaches_chi5_and_branch_ties_use_atom_names() {
    let component = component(
        &[
            "N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH2", "NH1",
        ],
        &[
            ("N", "CA", BondOrder::Single),
            ("CA", "C", BondOrder::Single),
            ("C", "O", BondOrder::Double),
            ("CA", "CB", BondOrder::Single),
            ("CB", "CG", BondOrder::Single),
            ("CG", "CD", BondOrder::Single),
            ("CD", "NE", BondOrder::Single),
            ("NE", "CZ", BondOrder::Single),
            ("CZ", "NH2", BondOrder::Aromatic),
            ("CZ", "NH1", BondOrder::Aromatic),
        ],
    );
    let Some(definition) = side_chain_definition(&component) else {
        panic!("definition absent")
    };
    assert_eq!(definition.torsion_count(), 5);
    assert_eq!(
        definition.atoms.iter().map(Box::as_ref).collect::<Vec<_>>(),
        ["N", "CA", "CB", "CG", "CD", "NE", "CZ", "NH1"]
    );
}

#[test]
fn aromatic_bond_terminates_after_the_conventional_second_torsion() {
    let component = component(
        &["N", "CA", "CB", "CG", "CD1", "CE1"],
        &[
            ("N", "CA", BondOrder::Single),
            ("CA", "CB", BondOrder::Single),
            ("CB", "CG", BondOrder::Single),
            ("CG", "CD1", BondOrder::Aromatic),
            ("CD1", "CE1", BondOrder::Aromatic),
        ],
    );
    assert_eq!(
        side_chain_definition(&component).map(|value| value.torsion_count()),
        Some(2)
    );
}

fn component(names: &[&str], bonds: &[(&str, &str, BondOrder)]) -> Component {
    Component {
        id: "TEST".into(),
        name: "test".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        formula: None,
        atoms: names
            .iter()
            .map(|name| ComponentAtom {
                name: (*name).into(),
                alternate_name: None,
                element: if name.starts_with('O') {
                    Element::OXYGEN
                } else if name.starts_with('N') {
                    Element::NITROGEN
                } else {
                    Element::CARBON
                },
                charge: 0,
                aromatic: false,
                leaving: false,
                stereo: None,
            })
            .collect(),
        bonds: bonds
            .iter()
            .map(|(left, right, order)| ComponentBond {
                atom_a: (*left).into(),
                atom_b: (*right).into(),
                order: *order,
                aromatic: *order == BondOrder::Aromatic,
                stereo: None,
            })
            .collect(),
        ideal_coordinates: Some(Arc::from([])),
        model_coordinates: None,
    }
}
