use super::{
    PeoeAtom, PeoeAtomType, PeoeBond, PeoeError, PeoeOptions, component_peoe_charges, peoe_charges,
};
use crate::{Component, ComponentAtom, ComponentBond, ComponentKind};
use pdbiox_core::{BondOrder, Element};

fn charges(atom_types: &[PeoeAtomType], bonds: &[(usize, usize)]) -> Vec<f64> {
    let atoms: Vec<_> = atom_types
        .iter()
        .map(|atom_type| PeoeAtom {
            atom_type: *atom_type,
            formal_charge: 0.0,
        })
        .collect();
    let bonds: Vec<_> = bonds
        .iter()
        .map(|&(atom_a, atom_b)| PeoeBond { atom_a, atom_b })
        .collect();
    match peoe_charges(&atoms, &bonds, PeoeOptions::default()) {
        Ok(charges) => charges,
        Err(error) => panic!("valid PEOE input failed: {error}"),
    }
}

#[test]
fn carbonyl_oxygen_receives_negative_charge() {
    let charges = charges(&[PeoeAtomType::CSp2, PeoeAtomType::OSp2], &[(0, 1)]);
    assert!(charges[1] < 0.0 && charges[0] > 0.0);
    assert!((charges[0] + charges[1]).abs() < 1.0e-9);
}

#[test]
fn orbital_state_changes_the_computed_charge() {
    let alkane = charges(&[PeoeAtomType::CSp3, PeoeAtomType::OSp3], &[(0, 1)]);
    let carbonyl = charges(&[PeoeAtomType::CSp2, PeoeAtomType::OSp2], &[(0, 1)]);
    assert!((alkane[0] - carbonyl[0]).abs() > 1.0e-3);
}

#[test]
fn formal_charge_is_conserved() {
    let atoms = [
        PeoeAtom {
            atom_type: PeoeAtomType::NSp2,
            formal_charge: 1.0,
        },
        PeoeAtom {
            atom_type: PeoeAtomType::CSp2,
            formal_charge: 0.0,
        },
    ];
    let result = peoe_charges(
        &atoms,
        &[PeoeBond {
            atom_a: 0,
            atom_b: 1,
        }],
        PeoeOptions::default(),
    );
    let Ok(charges) = result else {
        panic!("charged component should be supported")
    };
    assert!((charges.iter().sum::<f64>() - 1.0).abs() < 1.0e-9);
}

#[test]
fn component_topology_drives_hybridisation_perception() {
    let component = component(
        &[("C", Element::CARBON), ("O", Element::OXYGEN)],
        &[("C", "O", BondOrder::Double)],
    );
    let result = component_peoe_charges(&component, PeoeOptions::default());
    let Ok(charges) = result else {
        panic!("CCD carbonyl should be parameterised")
    };
    assert!(charges[1] < 0.0 && charges[0] > 0.0);
}

#[test]
fn unparameterised_elements_are_explicit_errors() {
    let component = component(&[("FE", Element::IRON)], &[]);
    assert!(matches!(
        component_peoe_charges(&component, PeoeOptions::default()),
        Err(PeoeError::UnsupportedAtom { atom: 0, .. })
    ));
}

#[test]
fn default_parameter_profile_has_stable_provenance_identity() {
    let profile = PeoeOptions::default().profile;
    assert_eq!(profile.name(), "gasteiger-marsili");
    assert_eq!(profile.version(), "1");
}

fn component(atoms: &[(&str, Element)], bonds: &[(&str, &str, BondOrder)]) -> Component {
    Component {
        id: "TEST".into(),
        name: "test component".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: atoms
            .iter()
            .map(|(name, element)| ComponentAtom {
                name: (*name).into(),
                alternate_name: None,
                element: *element,
                charge: 0,
                aromatic: false,
                leaving: false,
                stereo: None,
            })
            .collect::<Vec<_>>()
            .into(),
        bonds: bonds
            .iter()
            .map(|(atom_a, atom_b, order)| ComponentBond {
                atom_a: (*atom_a).into(),
                atom_b: (*atom_b).into(),
                order: *order,
                aromatic: *order == BondOrder::Aromatic,
                stereo: None,
            })
            .collect::<Vec<_>>()
            .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}
