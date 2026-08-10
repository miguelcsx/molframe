use super::{SideChainRoles, side_chain_definition};
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
    let side_chain_atoms = side_chain_names(&component);
    let Some(definition) = side_chain_definition(&component, &roles(&side_chain_atoms)) else {
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
    let side_chain_atoms = side_chain_names(&component);
    assert_eq!(
        side_chain_definition(&component, &roles(&side_chain_atoms))
            .map(|value| value.torsion_count()),
        Some(2)
    );
}

#[test]
fn explicit_roles_support_nonstandard_component_atom_identifiers() {
    let component = component(
        &["amide", "centre", "branch", "outer"],
        &[
            ("amide", "centre", BondOrder::Single),
            ("centre", "branch", BondOrder::Single),
            ("branch", "outer", BondOrder::Single),
        ],
    );
    let explicit = ["branch", "outer"];
    let definition = side_chain_definition(
        &component,
        &SideChainRoles {
            nitrogen: "amide",
            alpha_carbon: "centre",
            side_chain_atoms: &explicit,
        },
    );
    let Some(definition) = definition else {
        panic!("definition absent")
    };
    assert_eq!(
        definition.atoms.iter().map(Box::as_ref).collect::<Vec<_>>(),
        ["amide", "centre", "branch", "outer"]
    );
}

fn side_chain_names(component: &Component) -> Vec<&str> {
    component
        .atoms
        .iter()
        .filter(|atom| !matches!(atom.name.as_ref(), "N" | "CA" | "C" | "O"))
        .map(|atom| atom.name.as_ref())
        .collect()
}

fn roles<'a>(side_chain_atoms: &'a [&'a str]) -> SideChainRoles<'a> {
    SideChainRoles {
        nitrogen: "N",
        alpha_carbon: "CA",
        side_chain_atoms,
    }
}

fn component(names: &[&str], bonds: &[(&str, &str, BondOrder)]) -> Component {
    Component {
        id: "TEST".into(),
        name: "test".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: None,
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
