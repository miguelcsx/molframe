use super::*;
use crate::{ComponentAtom, ComponentBond, MemoryProvider, StereoConfiguration};
use pdbiox_core::Element;
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::io::{InputBuffer, ReadOptions};

const STRUCTURE: &str = "data_gly\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 N N GLY A 1 0.0 0 0\nATOM 2 C CA GLY A 1 1.4 0 0\n\
ATOM 3 C C GLY A 1 2.8 0 0\nATOM 4 N N GLY A 2 4.1 0 0\n\
ATOM 5 C CA GLY A 2 5.5 0 0\nATOM 6 C C GLY A 2 6.9 0 0\n";

fn component() -> Component {
    let atom = |name: &str, element| ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo: (name == "CA").then_some(StereoConfiguration::S),
    };
    Component {
        id: "GLY".into(),
        name: "GLYCINE".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        formula: Some("C2 H5 N O2".into()),
        atoms: vec![
            atom("N", Element::NITROGEN),
            atom("CA", Element::CARBON),
            atom("C", Element::CARBON),
        ]
        .into(),
        bonds: vec![
            ComponentBond {
                atom_a: "N".into(),
                atom_b: "CA".into(),
                order: BondOrder::Single,
                aromatic: false,
                stereo: None,
            },
            ComponentBond {
                atom_a: "CA".into(),
                atom_b: "C".into(),
                order: BondOrder::Single,
                aromatic: false,
                stereo: None,
            },
        ]
        .into(),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}

#[test]
fn ccd_edges_polymer_links_and_chain_classification_are_applied_together() {
    let input = InputBuffer::from_bytes(STRUCTURE.as_bytes().to_vec());
    let (structure, _) = pdbiox_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider = MemoryProvider::new(DictionaryVersion::new("test-1"), [component()]);
    let report = apply_component_chemistry(&structure, &provider).expect("provider works");
    assert_eq!(report.structure.data().bonds.len(), 5);
    assert_eq!(
        report
            .structure
            .chain(pdbiox_core::ChainIndex::new(0))
            .map(pdbiox_core::structure::ChainRef::polymer_kind),
        Some(PolymerKind::Protein)
    );
    assert_eq!(report.dictionary_version.as_str(), "test-1");
    assert!(report.findings.is_empty());
}

#[test]
fn an_unknown_component_is_reported_once_and_yields_a_known_empty_graph() {
    let input = InputBuffer::from_bytes(STRUCTURE.as_bytes().to_vec());
    let (structure, _) = pdbiox_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider = MemoryProvider::new(DictionaryVersion::new("empty"), []);
    let report = apply_component_chemistry(&structure, &provider).expect("provider works");
    assert!(report.structure.data().bonds.is_available());
    assert_eq!(
        report
            .findings
            .iter()
            .filter(|finding| finding.code() == Code::W3201)
            .count(),
        1
    );
}
