use super::{ChiralityIssue, ChiralityOptions, chirality_outliers};
use molframe_chem::{
    Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider, StereoConfiguration,
};
use molframe_core::contract::DictionaryVersion;
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::{AnalysisPolicy, BondOrder, Element, Structure};
use std::sync::Arc;

const HEADER: &str = "data_c\n\
loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

const REFERENCE: [[f32; 3]; 4] = [
    [-0.966, 0.493, 1.500],
    [0.257, 0.418, 0.692],
    [-0.094, 0.017, -0.716],
    [1.204, -0.620, 1.296],
];

fn structure(rows: &str) -> Structure {
    let input = InputBuffer::from_bytes(format!("{HEADER}{rows}").into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn provider() -> MemoryProvider {
    let names = ["N1", "CT", "C1", "B1"];
    let elements = [
        Element::NITROGEN,
        Element::CARBON,
        Element::CARBON,
        Element::CARBON,
    ];
    let atoms = names
        .into_iter()
        .zip(elements)
        .enumerate()
        .map(|(index, (name, element))| ComponentAtom {
            name: name.into(),
            alternate_name: None,
            element,
            charge: 0,
            aromatic: false,
            leaving: false,
            stereo: (index == 1).then_some(StereoConfiguration::S),
        })
        .collect::<Vec<_>>();
    let bonds = ["N1", "C1", "B1"]
        .into_iter()
        .map(|name| ComponentBond {
            atom_a: "CT".into(),
            atom_b: name.into(),
            order: BondOrder::Single,
            aromatic: false,
            stereo: None,
        })
        .collect::<Vec<_>>();
    MemoryProvider::new(
        DictionaryVersion::new("ccd-test"),
        [Component {
            id: "CMP".into(),
            name: "test".into(),
            kind: ComponentKind::NonPolymer,
            parent: None,
            one_letter_code: None,
            formula: None,
            atoms: atoms.into(),
            bonds: bonds.into(),
            ideal_coordinates: Some(Arc::from(REFERENCE)),
            model_coordinates: None,
        }],
    )
    .expect("component fixture is unique")
}

fn options() -> ChiralityOptions {
    ChiralityOptions {
        minimum_abs_volume: 0.01,
    }
}

#[test]
fn ccd_reference_geometry_accepts_matching_handedness_without_standard_atom_names() {
    let rows = "1 N N1 CMP A 1 -0.966 0.493 1.500\n\
2 C CT CMP A 1 0.257 0.418 0.692\n\
3 C C1 CMP A 1 -0.094 0.017 -0.716\n\
4 C B1 CMP A 1 1.204 -0.620 1.296\n";
    let report = chirality_outliers(
        &structure(rows),
        &provider(),
        &AnalysisPolicy::default(),
        options(),
    )
    .expect("CCD validation succeeds");
    assert!(report.flags.is_empty());
    assert_eq!(report.dictionary_version.as_str(), "ccd-test");
    assert_eq!((report.intended, report.assessed), (1, 1));
}

#[test]
fn reflection_is_reported_against_the_same_ccd_neighbour_order() {
    let rows = "1 N N1 CMP A 1 0.966 0.493 1.500\n\
2 C CT CMP A 1 -0.257 0.418 0.692\n\
3 C C1 CMP A 1 0.094 0.017 -0.716\n\
4 C B1 CMP A 1 -1.204 -0.620 1.296\n";
    let report = chirality_outliers(
        &structure(rows),
        &provider(),
        &AnalysisPolicy::default(),
        options(),
    )
    .expect("CCD validation succeeds");
    assert_eq!(report.flags.len(), 1);
    assert_eq!(report.flags[0].issue, ChiralityIssue::Inverted);
    assert!(report.flags[0].observed_volume * report.flags[0].reference_volume < 0.0);
}

#[test]
fn missing_component_is_recorded_without_a_compatibility_fallback() {
    let rows = "1 C CT UNKNOWN A 1 0 0 0\n";
    let report = chirality_outliers(
        &structure(rows),
        &provider(),
        &AnalysisPolicy::default(),
        options(),
    )
    .expect("unknown chemistry is non-fatal");
    assert!(report.flags.is_empty());
    assert_eq!(report.findings.len(), 1);
    assert_eq!((report.intended, report.assessed), (1, 0));
}
