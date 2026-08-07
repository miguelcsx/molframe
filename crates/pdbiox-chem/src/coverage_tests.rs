use super::component_coverage;
use crate::{Component, ComponentAtom, ComponentKind, MemoryProvider};
use pdbiox_core::contract::{AltlocPolicy, AnalysisPolicy, DictionaryVersion, HydrogenPolicy};
use pdbiox_core::{Element, InputBuffer, ReadOptions};
use std::sync::Arc;

const ENTRY: &str = "data_c\n\
loop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n\
_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n\
_atom_site.label_alt_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
1 N N GLY A 1 . 0 0 0\n\
2 C CA GLY A 1 A 1 0 0\n\
3 C CA GLY A 1 B 1 0 0\n";

#[test]
fn missing_and_ambiguous_counts_come_from_expected_ccd_atoms() {
    let structure = structure();
    let provider = provider();
    let keep_all = AnalysisPolicy {
        altloc: AltlocPolicy::KeepAll,
        hydrogens: HydrogenPolicy::IncludeInferred,
        ..AnalysisPolicy::default()
    };
    let report = match component_coverage(&structure, &provider, &keep_all) {
        Ok(report) => report,
        Err(finding) => panic!("coverage failed: {finding}"),
    };
    assert_eq!(report.coverage.intended, 4);
    assert_eq!(report.coverage.used, 1);
    assert_eq!(report.coverage.missing, 2);
    assert_eq!(report.coverage.ambiguous, 1);
    assert!((report.coverage.fraction() - 0.25).abs() < f32::EPSILON);
}

#[test]
fn policy_resolves_altlocs_and_defines_hydrogen_intent() {
    let structure = structure();
    let provider = provider();
    let policy = AnalysisPolicy {
        hydrogens: HydrogenPolicy::ExplicitOnly,
        ..AnalysisPolicy::default()
    };
    let report = match component_coverage(&structure, &provider, &policy) {
        Ok(report) => report,
        Err(finding) => panic!("coverage failed: {finding}"),
    };
    assert_eq!(report.coverage.intended, 3);
    assert_eq!(report.coverage.used, 2);
    assert_eq!(report.coverage.missing, 1);
    assert_eq!(report.coverage.ambiguous, 0);
}

fn structure() -> pdbiox_core::Structure {
    let input = InputBuffer::from_bytes(ENTRY.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn provider() -> MemoryProvider {
    MemoryProvider::new(
        DictionaryVersion::new("test"),
        [Component {
            id: "GLY".into(),
            name: "glycine".into(),
            kind: ComponentKind::AminoAcid,
            parent: None,
            formula: None,
            atoms: Arc::from([
                atom("N", Element::NITROGEN),
                atom("CA", Element::CARBON),
                atom("C", Element::CARBON),
                atom("H", Element::HYDROGEN),
            ]),
            bonds: Arc::from([]),
            ideal_coordinates: None,
            model_coordinates: None,
        }],
    )
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
