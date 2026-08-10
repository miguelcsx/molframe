use super::*;
use pdbiox_chem::{
    Component, ComponentAtom, ComponentKind, MemoryProvider, PolymerAtomRole, PolymerRoleProfile,
    PolymerRoleRule, apply_polymer_role_profile,
};
use pdbiox_core::Element;
use pdbiox_core::contract::DictionaryVersion;
use pdbiox_core::io::{InputBuffer, ReadOptions};
use std::sync::Arc;

const HEADER: &str = "data_n\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
const ATOMS: &str = "ATOM 1 O O4' C A 1 0 1 0\nATOM 2 C C1' C A 1 1 1 0.1\nATOM 3 C C2' C A 1 1.5 0 0\nATOM 4 C C3' C A 1 0.5 -1 0.1\nATOM 5 C C4' C A 1 -0.5 0 0\nATOM 6 N N1 C A 1 2.3 1 0.1\nATOM 7 C C2 C A 1 3.3 1 0.1\nATOM 8 N N3 C A 1 3.8 2 0.1\nATOM 9 C C4 C A 1 3.3 3 0.1\nATOM 10 C C5 C A 1 2.3 3 0.1\nATOM 11 C C6 C A 1 1.8 2 0.1\n";

#[test]
fn complete_nucleotide_reports_pucker_planarity_and_glycosidic_length() {
    let (structure, provider, profile) = context(ATOMS);
    let Ok(records) = nucleic_acid_geometry(&structure, &provider, &profile, policy()) else {
        panic!("valid explicit nucleic geometry");
    };
    assert_eq!(records.len(), 1);
    assert!(records[0].pucker.is_some());
    assert!(
        records[0]
            .base_plane_deviation
            .is_some_and(|value| value < 1.0e-6)
    );
    assert!(records[0].glycosidic_bond_length.is_some());
    assert!(!records[0].issues.iter().any(|issue| matches!(
        issue,
        NucleicGeometryIssue::MissingBaseAtoms { .. }
            | NucleicGeometryIssue::MissingSugarAtoms { .. }
    )));
}

#[test]
fn ccd_expected_and_observed_role_atoms_report_missing_and_warped_base() {
    let changed = ATOMS.replace("ATOM 9 C C4 C A 1 3.3 3 0.1\n", "").replace(
        "ATOM 10 C C5 C A 1 2.3 3 0.1",
        "ATOM 10 C C5 C A 1 2.3 3 1.5",
    );
    let (structure, provider, profile) = context(&changed);
    let mut strict = policy();
    strict.maximum_base_plane_deviation = 0.01;
    let Ok(records) = nucleic_acid_geometry(&structure, &provider, &profile, strict) else {
        panic!("valid incomplete nucleic geometry");
    };
    assert!(
        records[0]
            .issues
            .iter()
            .any(|issue| matches!(issue, NucleicGeometryIssue::MissingBaseAtoms { .. }))
    );
    assert!(
        records[0]
            .issues
            .iter()
            .any(|issue| matches!(issue, NucleicGeometryIssue::NonPlanarBase { .. }))
    );
}

#[test]
fn thresholds_have_no_hidden_default_or_reversal() {
    let (structure, provider, profile) = context(ATOMS);
    let mut invalid = policy();
    invalid.glycosidic_bond_range = [2.0, 1.0];
    assert!(matches!(
        nucleic_acid_geometry(&structure, &provider, &profile, invalid),
        Err(NucleicGeometryError::InvalidPolicy)
    ));
}

fn context(atoms: &str) -> (pdbiox_core::Structure, MemoryProvider, PolymerRoleProfile) {
    let input = InputBuffer::from_bytes(format!("{HEADER}{atoms}").into_bytes());
    let structure = match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    };
    let provider = MemoryProvider::new(DictionaryVersion::new("ccd-test"), [component()])
        .expect("component fixture is unique");
    let profile = profile();
    let structure = match apply_polymer_role_profile(&structure, &provider, &profile) {
        Ok(report) => report.structure,
        Err(error) => panic!("role application failed: {error}"),
    };
    (structure, provider, profile)
}

fn policy() -> NucleicGeometryPolicy {
    NucleicGeometryPolicy {
        maximum_base_plane_deviation: 0.15,
        glycosidic_bond_range: [1.2, 1.7],
        phosphodiester_bond_range: [1.3, 2.0],
        plane_fit: pdbiox_geom::EigenOptions::standard(),
    }
}

fn profile() -> PolymerRoleProfile {
    let assignments = [
        ("O4'", PolymerAtomRole::NUCLEIC_O4),
        ("C1'", PolymerAtomRole::NUCLEIC_C1),
        ("C2'", PolymerAtomRole::NUCLEIC_C2),
        ("C3'", PolymerAtomRole::NUCLEIC_C3),
        ("C4'", PolymerAtomRole::NUCLEIC_C4),
        ("N1", PolymerAtomRole::NUCLEIC_GLYCOSIDIC),
        ("C2", PolymerAtomRole::NUCLEIC_BASE_REFERENCE),
        ("N3", PolymerAtomRole::NUCLEIC_BASE),
        ("C4", PolymerAtomRole::NUCLEIC_BASE),
        ("C5", PolymerAtomRole::NUCLEIC_BASE),
        ("C6", PolymerAtomRole::NUCLEIC_BASE),
    ];
    PolymerRoleProfile {
        id: "nucleic-test-v1".into(),
        rules: assignments
            .into_iter()
            .map(|(atom_name, role)| PolymerRoleRule {
                component_id: Some("C".into()),
                component_kind: Some(ComponentKind::Nucleotide),
                atom_name: atom_name.into(),
                role,
            })
            .collect(),
    }
}

fn component() -> Component {
    let names = [
        "O4'", "C1'", "C2'", "C3'", "C4'", "N1", "C2", "N3", "C4", "C5", "C6",
    ];
    Component {
        id: "C".into(),
        name: "test cytidine".into(),
        kind: ComponentKind::Nucleotide,
        parent: None,
        one_letter_code: Some(b'C'),
        formula: None,
        atoms: names
            .into_iter()
            .map(|name| ComponentAtom {
                name: name.into(),
                alternate_name: None,
                element: if name.starts_with('N') {
                    Element::NITROGEN
                } else if name.starts_with('O') {
                    Element::OXYGEN
                } else {
                    Element::CARBON
                },
                charge: 0,
                aromatic: false,
                leaving: false,
                stereo: None,
            })
            .collect(),
        bonds: Arc::from([]),
        ideal_coordinates: None,
        model_coordinates: None,
    }
}
