use super::{BasePairOptions, base_pairs};
use crate::HydrogenBondOptions;
use pdbiox_chem::{
    Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider, PolymerLinkPolicy,
    apply_component_chemistry,
};
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::{BondOrder, Element, ExecutionContext};
use pdbiox_spatial::SpatialBackend;
use std::sync::Arc;

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N D ADE A 1 0 0 0\n\
ATOM 2 H H ADE A 1 1 0 0\n\
ATOM 3 O A URA B 1 2.8 0 0\n";

#[test]
fn ccd_identity_and_oriented_bond_define_a_pair() {
    let structure = annotated_structure();
    let Ok(pairs) = base_pairs(
        &structure,
        &provider(b'U'),
        options(),
        &ExecutionContext::default(),
    ) else {
        panic!("valid base-pair analysis");
    };
    assert_eq!(pairs.len(), 1);
    assert_eq!(pairs[0].hydrogen_bond_count, 1);
}

#[test]
fn non_complementary_ccd_codes_do_not_pair() {
    let structure = annotated_structure();
    let Ok(pairs) = base_pairs(
        &structure,
        &provider(b'G'),
        options(),
        &ExecutionContext::default(),
    ) else {
        panic!("valid non-complementary analysis");
    };
    assert!(pairs.is_empty());
}

#[test]
fn supporting_bond_count_is_an_explicit_policy() {
    let structure = annotated_structure();
    let mut options = options();
    options.minimum_hydrogen_bonds = 2;
    let Ok(pairs) = base_pairs(
        &structure,
        &provider(b'U'),
        options,
        &ExecutionContext::default(),
    ) else {
        panic!("valid stringent analysis");
    };
    assert!(pairs.is_empty());
}

fn annotated_structure() -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    let structure = match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    match apply_component_chemistry(&structure, &provider(b'U'), PolymerLinkPolicy::Disabled) {
        Ok(report) => report.structure,
        Err(error) => panic!("chemistry failed: {error}"),
    }
}

fn options() -> BasePairOptions {
    BasePairOptions {
        hydrogen_bonds: HydrogenBondOptions {
            maximum_donor_acceptor_distance: 3.5,
            minimum_angle_degrees: 150.0,
            backend: SpatialBackend::BruteForce,
            periodic: false,
        },
        minimum_hydrogen_bonds: 1,
    }
}

fn provider(second_code: u8) -> MemoryProvider {
    MemoryProvider::new(
        DictionaryVersion::new("test"),
        [donor_component(), acceptor_component(second_code)],
    )
    .expect("component fixtures are unique")
}

fn donor_component() -> Component {
    Component {
        id: "ADE".into(),
        name: "adenine test".into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(b'A'),
        formula: None,
        atoms: Arc::from([atom("D", Element::NITROGEN), atom("H", Element::HYDROGEN)]),
        bonds: Arc::from([ComponentBond {
            atom_a: "D".into(),
            atom_b: "H".into(),
            order: BondOrder::Single,
            aromatic: false,
            stereo: None,
        }]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

fn acceptor_component(code: u8) -> Component {
    Component {
        id: "URA".into(),
        name: "acceptor test".into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(code),
        formula: None,
        atoms: Arc::from([atom("A", Element::OXYGEN)]),
        bonds: Arc::from([]),
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
