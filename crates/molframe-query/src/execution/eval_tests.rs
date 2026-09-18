use super::*;
use molframe_core::io::{InputBuffer, ReadOptions};

const SOURCE: &str = "data_q\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.auth_seq_id\n\
_atom_site.auth_asym_id\n\
ATOM 1 N N  GLY LONG 1 0 0 0 1.0 10 42 A\n\
ATOM 2 C CA GLY LONG 1 1 0 0 0.5 40 42 A\n\
HETATM 3 O O HOH W 2 3 0 0 1.0 20 7 Z\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn macro_structure() -> Structure {
    let structure = structure();
    let mut data = structure.data().clone();
    data.annotations.insert(
        molframe_core::COMPONENT_KIND_ANNOTATION,
        molframe_core::AtomAnnotation::Integer(
            molframe_core::AnnotationColumn::from_values(vec![
                molframe_chem::ComponentKind::AminoAcid.code(),
                molframe_chem::ComponentKind::AminoAcid.code(),
                molframe_chem::ComponentKind::Solvent.code(),
            ])
            .expect("small annotation column"),
        ),
    );
    data.annotations.insert(
        molframe_core::POLYMER_ATOM_ROLE_ANNOTATION,
        molframe_core::AtomAnnotation::Integer(
            molframe_core::AnnotationColumn::from_values(vec![
                molframe_chem::PolymerAtomRole::PROTEIN_NITROGEN.code(),
                molframe_chem::PolymerAtomRole::PROTEIN_ALPHA_CARBON.code(),
                molframe_chem::PolymerAtomRole::UNKNOWN.code(),
            ])
            .expect("small annotation column"),
        ),
    );
    Structure::new(data)
}

fn evaluate(source: &str) -> Evaluation {
    let query = match Query::compile(source) {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    match query.evaluate(
        &macro_structure(),
        &AnalysisPolicy::default().with_identifiers(molframe_core::contract::Namespace::Label),
        &Groups::new(),
        None,
    ) {
        Ok(evaluation) => evaluation,
        Err(findings) => panic!("evaluation failed: {findings:?}"),
    }
}

#[test]
fn boolean_membership_ranges_and_numeric_comparisons_execute() {
    assert_eq!(
        evaluate("name N CA and bfactor >= 30")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![1]
    );
    assert_eq!(
        evaluate("index 0:1 and not occupancy < 0.75")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![0]
    );
}

#[test]
fn hierarchy_expansion_and_macros_execute_over_one_store() {
    assert_eq!(
        evaluate("byres name CA")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(
        evaluate("same chain as name CA")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(
        evaluate("water").selection.iter().collect::<Vec<_>>(),
        vec![2]
    );
    assert_eq!(
        evaluate("protein and backbone")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
}

#[test]
fn a_chemistry_macro_refuses_absent_component_chemistry() {
    let query = Query::compile("protein").expect("query compiles");
    let result = query.evaluate(
        &structure(),
        &AnalysisPolicy::default(),
        &Groups::new(),
        None,
    );
    assert!(result.is_err());
}

#[test]
fn polymer_macro_does_not_infer_chemistry_from_atom_record_type() {
    let query = Query::compile("polymer").expect("query compiles");
    let evaluation = query
        .evaluate(
            &structure(),
            &AnalysisPolicy::default(),
            &Groups::new(),
            None,
        )
        .expect("topology query evaluates");

    assert!(evaluation.selection.is_empty());
}

#[test]
fn chemistry_numeric_columns_and_aromatic_edges_use_shared_rust_data() {
    assert_eq!(
        evaluate("mass > 13 and radius < 1.6")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![0, 2]
    );

    let mut data = structure().data().clone();
    let mut bonds = molframe_core::BondTableBuilder::new();
    bonds.push(molframe_core::BondRecord {
        atom_a: molframe_core::AtomIndex::new(0),
        atom_b: molframe_core::AtomIndex::new(1),
        order: molframe_core::BondOrder::Aromatic,
        provenance: molframe_core::BondProvenance::ChemicalComponentDictionary,
    });
    data.bonds = bonds.finish();
    let structure = Structure::new(data);
    let query = Query::compile("aromatic").expect("query compiles");
    let result = query
        .evaluate(&structure, &AnalysisPolicy::default(), &Groups::new(), None)
        .expect("annotations exist");
    assert_eq!(result.selection.iter().collect::<Vec<_>>(), vec![0, 1]);
}

#[test]
fn chemistry_macros_aromatic_charge_and_chirality_use_ccd_annotations() {
    let structure = chemistry_structure();
    for (source, expected) in [
        ("ligand", vec![0, 1, 2, 3]),
        ("aromatic", vec![2]),
        ("formalcharge < 0", vec![3]),
        ("chirality S", vec![0]),
        ("smarts '[C;D3]-[S-]'", vec![0, 3]),
    ] {
        let query = match Query::compile(source) {
            Ok(query) => query,
            Err(findings) => panic!("{source} compile failed: {findings:?}"),
        };
        let evaluation =
            match query.evaluate(&structure, &AnalysisPolicy::default(), &Groups::new(), None) {
                Ok(evaluation) => evaluation,
                Err(findings) => panic!("{source} evaluation failed: {findings:?}"),
            };
        assert_eq!(evaluation.selection.iter().collect::<Vec<_>>(), expected);
    }
}

fn chemistry_structure() -> Structure {
    use molframe_chem::{
        Component, ComponentAtom, ComponentBond, ComponentKind, MemoryProvider, StereoConfiguration,
    };
    use molframe_core::{BondOrder, Element};

    let source = "data_chem\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
HETATM 1 C CTR LIG A 1 0 0 0\n\
HETATM 2 N A LIG A 1 1 0 0\n\
HETATM 3 O B LIG A 1 0 1 0\n\
HETATM 4 S C LIG A 1 0 0 -1\n";
    let atom = |name: &str, element, charge, aromatic, stereo| ComponentAtom {
        name: name.into(),
        alternate_name: None,
        element,
        charge,
        aromatic,
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
            atom(
                "CTR",
                Element::CARBON,
                0,
                false,
                Some(StereoConfiguration::R),
            ),
            atom("A", Element::NITROGEN, 0, false, None),
            atom("B", Element::OXYGEN, 0, true, None),
            atom("C", Element::SULFUR, -1, false, None),
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
    let (structure, _) = match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok(result) => result,
        Err(findings) => panic!("chemistry fixture failed: {findings:?}"),
    };
    let provider = MemoryProvider::new(
        molframe_core::contract::DictionaryVersion::new("query-test"),
        [component],
    )
    .expect("component fixture is unique");
    match molframe_chem::apply_component_chemistry(
        &structure,
        &provider,
        molframe_chem::PolymerLinkPolicy::Disabled,
    ) {
        Ok(report) => report.structure,
        Err(finding) => panic!("chemistry annotation failed: {finding}"),
    }
}

#[test]
fn both_namespaces_and_named_groups_are_resolved_explicitly() {
    assert_eq!(
        evaluate("label_chain LONG")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(
        evaluate("auth_chain Z")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![2]
    );

    let query = match Query::compile("group pocket") {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    let groups = Groups::from([(
        Box::<str>::from("pocket"),
        AtomSelection::from_sorted(vec![1]),
    )]);
    let result = query.evaluate(&structure(), &AnalysisPolicy::default(), &groups, None);
    let result = match result {
        Ok(result) => result,
        Err(findings) => panic!("evaluation failed: {findings:?}"),
    };
    assert_eq!(result.selection.iter().collect::<Vec<_>>(), vec![1]);
}

#[test]
fn explicit_namespace_policy_rejects_unqualified_identity_selectors() {
    let query = match Query::compile("chain A") {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    let policy =
        AnalysisPolicy::default().with_identifiers(molframe_core::contract::Namespace::Explicit);
    let result = query.evaluate(&structure(), &policy, &Groups::new(), None);
    assert_eq!(
        result
            .err()
            .and_then(|findings| findings.first().map(Diagnostic::code)),
        Some(Code::E6001)
    );
}

#[test]
fn unavailable_semantics_fail_instead_of_returning_a_guess() {
    for source in ["bonded name CA", "aromatic", "smarts 'C=O'"] {
        let query = match Query::compile(source) {
            Ok(query) => query,
            Err(findings) => panic!("compile failed: {findings:?}"),
        };
        let result = query.evaluate(
            &structure(),
            &AnalysisPolicy::default(),
            &Groups::new(),
            None,
        );
        assert_eq!(
            result
                .err()
                .and_then(|findings| findings.first().map(Diagnostic::code)),
            Some(Code::E4003)
        );
    }
}

#[test]
fn segment_atom_and_predicted_annotations_execute_from_core_columns() {
    let source = structure();
    let mut data = source.data().clone();
    let system = match data.dictionary.intern("SYSTEM") {
        Ok(symbol) => symbol,
        Err(error) => panic!("dictionary failed: {error}"),
    };
    let solvent = match data.dictionary.intern("SOLVENT") {
        Ok(symbol) => symbol,
        Err(error) => panic!("dictionary failed: {error}"),
    };
    let _ = data.annotations.insert(
        molframe_core::SEGMENT_ID_ANNOTATION,
        molframe_core::AtomAnnotation::Symbol(
            molframe_core::AnnotationColumn::from_values(vec![system, system, solvent])
                .expect("small annotation column"),
        ),
    );
    let _ = data.annotations.insert(
        molframe_core::PLDDT_ANNOTATION,
        molframe_core::AtomAnnotation::Real(
            molframe_core::AnnotationColumn::from_values(vec![90.0, 40.0, 80.0])
                .expect("small annotation column"),
        ),
    );
    let structure = Structure::new(data);
    for (source, expected) in [
        ("segid SYSTEM", vec![0, 1]),
        ("atom SYSTEM 1 CA", vec![1]),
        ("same segment as name CA", vec![0, 1]),
        ("plddt > 70", vec![0, 2]),
    ] {
        let query = match Query::compile(source) {
            Ok(query) => query,
            Err(findings) => panic!("compile failed: {findings:?}"),
        };
        let policy =
            AnalysisPolicy::default().with_identifiers(molframe_core::contract::Namespace::Label);
        let result = query.evaluate(&structure, &policy, &Groups::new(), None);
        let result = match result {
            Ok(result) => result,
            Err(findings) => panic!("evaluation failed: {findings:?}"),
        };
        assert_eq!(result.selection.iter().collect::<Vec<_>>(), expected);
    }
}
