use super::*;
use pdbiox_chem::{
    ComponentAtom, ComponentBond, ComponentKind, MemoryProvider, StereoConfiguration,
};
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::{BondOrder, Element};
use std::sync::Arc;

fn structure() -> Structure {
    let source = "data_g\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
HETATM 1 C A CMP X 1 0 0 0\n\
HETATM 2 C B CMP X 1 2 0 0\n\
HETATM 3 C C CMP X 1 2 2 0\n";
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let Ok((structure, _)) = pdbiox_cif::read(&input, &ReadOptions::new()) else {
        panic!("valid geometry fixture");
    };
    structure
}

fn atom(name: &str) -> ComponentAtom {
    ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element: Element::CARBON,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: Option::<StereoConfiguration>::None,
    }
}

fn provider() -> MemoryProvider {
    MemoryProvider::new(
        DictionaryVersion::new("ccd-test"),
        [Component {
            id: "CMP".into(),
            name: "component".into(),
            kind: ComponentKind::NonPolymer,
            parent: None,
            one_letter_code: None,
            formula: None,
            atoms: Arc::from([atom("A"), atom("B"), atom("C")]),
            bonds: Arc::from([
                ComponentBond {
                    atom_a: "A".into(),
                    atom_b: "B".into(),
                    order: BondOrder::Single,
                    aromatic: false,
                    stereo: None,
                },
                ComponentBond {
                    atom_a: "B".into(),
                    atom_b: "C".into(),
                    order: BondOrder::Single,
                    aromatic: false,
                    stereo: None,
                },
            ]),
            ideal_coordinates: Some(Arc::from([
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.5, 0.866_025_4, 0.0],
            ])),
            model_coordinates: None,
        }],
    )
    .expect("component fixture is unique")
}

#[test]
fn ccd_ideal_bonds_and_angles_share_explicit_coverage() {
    let Ok(report) = reference_geometry(
        &structure(),
        &provider(),
        Namespace::Label,
        ReferenceGeometryOptions {
            maximum_bond_deviation: 0.2,
            maximum_angle_deviation_degrees: 5.0,
        },
    ) else {
        panic!("CCD reference geometry should be assessable");
    };
    assert_eq!((report.intended, report.assessed), (3, 3));
    assert_eq!(report.bonds.len(), 2);
    assert_eq!(report.angles.len(), 1);
}

#[test]
fn explicit_namespace_is_rejected_instead_of_guessed() {
    let result = reference_geometry(
        &structure(),
        &provider(),
        Namespace::Explicit,
        ReferenceGeometryOptions {
            maximum_bond_deviation: 0.2,
            maximum_angle_deviation_degrees: 5.0,
        },
    );
    assert!(result.is_err());
}
