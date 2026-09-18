use super::ligand_symmetry_rmsd;
use molframe_chem::{Component, ComponentAtom, ComponentBond, ComponentKind};
use molframe_core::{BondOrder, Element};

fn symmetric_component() -> Component {
    Component {
        id: "CO2".into(),
        name: "carbon dioxide".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: Some("CO2".into()),
        atoms: vec![
            atom("C", Element::CARBON),
            atom("O1", Element::OXYGEN),
            atom("O2", Element::OXYGEN),
        ]
        .into(),
        bonds: vec![bond("C", "O1"), bond("C", "O2")].into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

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

fn bond(left: &str, right: &str) -> ComponentBond {
    ComponentBond {
        atom_a: left.into(),
        atom_b: right.into(),
        order: BondOrder::Double,
        aromatic: false,
        stereo: None,
    }
}

#[test]
fn swapped_equivalent_atoms_have_zero_ligand_rmsd() {
    let component = symmetric_component();
    let reference = [[0.0, 0.0, 0.0], [-1.2, 0.0, 0.0], [1.2, 0.0, 0.0]];
    let model = [[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [-1.2, 0.0, 0.0]];
    let result = ligand_symmetry_rmsd(&reference, &model, &component, 4)
        .unwrap_or_else(|error| panic!("ligand comparison failed: {error}"));
    assert!(result.rmsd < f64::EPSILON);
    assert_eq!(result.mapping.reference_to_model.as_ref(), &[0, 2, 1]);
}
