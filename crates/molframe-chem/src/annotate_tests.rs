use super::*;
use crate::{ComponentAtom, ComponentBond, MemoryProvider, StereoConfiguration};
use molframe_core::Element;
use molframe_core::contract::DictionaryVersion;
use molframe_core::io::{InputBuffer, ReadOptions};

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
        charge: if name == "N" { -1 } else { 0 },
        aromatic: name == "N",
        leaving: false,
        stereo: (name == "CA").then_some(StereoConfiguration::S),
    };
    Component {
        id: "GLY".into(),
        name: "GLYCINE".into(),
        kind: ComponentKind::AminoAcid,
        parent: None,
        one_letter_code: Some(b'G'),
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
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider = MemoryProvider::new(DictionaryVersion::new("test-1"), [component()])
        .expect("component fixture is unique");
    let report = apply_component_chemistry(
        &structure,
        &provider,
        PolymerLinkPolicy::explicit(
            2.1,
            vec![PolymerLinkRule::new(
                ComponentKind::AminoAcid,
                "C",
                ComponentKind::AminoAcid,
                "N",
            )],
        ),
    )
    .expect("provider works");
    assert_eq!(report.structure.data().bonds.len(), 5);
    assert_eq!(
        report
            .structure
            .chain(molframe_core::ChainIndex::new(0))
            .map(molframe_core::structure::ChainRef::polymer_kind),
        Some(PolymerKind::Protein)
    );
    assert_eq!(report.dictionary_version.as_str(), "test-1");
    assert!(report.findings.is_empty());
    let Some(molframe_core::AtomAnnotation::Boolean(aromatic)) = report
        .structure
        .annotations()
        .get(molframe_core::AROMATIC_ATOM_ANNOTATION)
    else {
        panic!("aromatic annotation absent")
    };
    assert_eq!(
        aromatic.get(0),
        Some((true, molframe_core::Presence::Present))
    );
    let Some(molframe_core::AtomAnnotation::Integer(charges)) = report
        .structure
        .annotations()
        .get(molframe_core::FORMAL_CHARGE_ANNOTATION)
    else {
        panic!("formal charge annotation absent")
    };
    assert_eq!(charges.get(0), Some((-1, molframe_core::Presence::Present)));
}

#[test]
fn modelled_stereo_is_anchored_to_the_component_reference_geometry() {
    let source = "data_stereo\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
HETATM 1 C CTR LIG A 1 0 0 0\n\
HETATM 2 N A LIG A 1 1 0 0\n\
HETATM 3 O B LIG A 1 0 1 0\n\
HETATM 4 S C LIG A 1 0 0 -1\n";
    let atom = |name: &str, element, stereo| ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge: 0,
        aromatic: false,
        leaving: false,
        stereo,
    };
    let bond = |other: &str| ComponentBond {
        atom_a: "CTR".into(),
        atom_b: other.into(),
        order: BondOrder::Single,
        aromatic: false,
        stereo: None,
    };
    let component = Component {
        id: "LIG".into(),
        name: "CHIRAL LIGAND".into(),
        kind: ComponentKind::NonPolymer,
        parent: None,
        one_letter_code: None,
        formula: None,
        atoms: vec![
            atom("CTR", Element::CARBON, Some(StereoConfiguration::R)),
            atom("A", Element::NITROGEN, None),
            atom("B", Element::OXYGEN, None),
            atom("C", Element::SULFUR, None),
        ]
        .into(),
        bonds: vec![bond("A"), bond("B"), bond("C")].into(),
        ideal_coordinates: Some(
            vec![
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
            .into(),
        ),
        model_coordinates: None,
    };
    let input = InputBuffer::from_bytes(source.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider = MemoryProvider::new(DictionaryVersion::new("stereo"), [component])
        .expect("component fixture is unique");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("chemistry applies");
    let Some(molframe_core::AtomAnnotation::Symbol(stereo)) = report
        .structure
        .annotations()
        .get(molframe_core::STEREO_CONFIGURATION_ANNOTATION)
    else {
        panic!("stereo annotation absent")
    };
    let Some((configuration, molframe_core::Presence::Present)) = stereo.get(0) else {
        panic!("centre configuration absent")
    };
    assert_eq!(report.structure.resolve(configuration), Some("S"));
}

#[test]
fn an_unknown_component_is_reported_once_and_yields_a_known_empty_graph() {
    let input = InputBuffer::from_bytes(STRUCTURE.as_bytes().to_vec());
    let (structure, _) = molframe_cif::read(&input, &ReadOptions::new()).expect("fixture reads");
    let provider =
        MemoryProvider::new(DictionaryVersion::new("empty"), []).expect("empty fixture is unique");
    let report = apply_component_chemistry(&structure, &provider, PolymerLinkPolicy::Disabled)
        .expect("provider works");
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
