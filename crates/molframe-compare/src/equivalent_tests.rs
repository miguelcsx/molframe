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

#[test]
fn a_named_comparison_uses_only_the_atoms_both_sets_have_and_their_symmetry() {
    use super::named_ligand_rmsd;
    let component = methyl_component();
    // The component defines a carbon and three equivalent hydrogens; the
    // coordinate sets name the carbon and only two of the hydrogens, swapped.
    let reference = [
        ("C", [0.0, 0.0, 0.0]),
        ("H1", [1.0, 0.0, 0.0]),
        ("H2", [0.0, 1.0, 0.0]),
    ];
    let model = [
        ("H2", [1.0, 0.0, 0.0]),
        ("C", [0.0, 0.0, 0.0]),
        ("H1", [0.0, 1.0, 0.0]),
        ("Zn", [9.0, 9.0, 9.0]),
    ];
    let result = named_ligand_rmsd(&reference, &model, &component, 8).expect("shared atoms");
    assert!(result.rmsd.abs() < 1e-6, "{}", result.rmsd);
}

#[test]
fn a_named_comparison_with_no_shared_atom_is_an_error() {
    use super::named_ligand_rmsd;
    let component = methyl_component();
    assert!(named_ligand_rmsd(&[("Q", [0.0; 3])], &[("Q", [0.0; 3])], &component, 8).is_err());
}

fn methyl_component() -> Component {
    Component {
        id: "MET".into(),
        name: "methyl".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: vec![
            atom("C", Element::CARBON),
            atom("H1", Element::HYDROGEN),
            atom("H2", Element::HYDROGEN),
            atom("H3", Element::HYDROGEN),
        ]
        .into(),
        bonds: vec![bond("C", "H1"), bond("C", "H2"), bond("C", "H3")].into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}
