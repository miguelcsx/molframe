use super::*;
use pdbiox_core::io::{InputBuffer, ReadOptions};

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
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn evaluate(source: &str) -> Evaluation {
    let query = match Query::compile(source) {
        Ok(query) => query,
        Err(findings) => panic!("compile failed: {findings:?}"),
    };
    match query.evaluate(
        &structure(),
        &AnalysisPolicy::default(),
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
fn chemistry_numeric_columns_and_aromatic_edges_use_shared_rust_data() {
    assert_eq!(
        evaluate("mass > 13 and radius < 1.6")
            .selection
            .iter()
            .collect::<Vec<_>>(),
        vec![0, 2]
    );

    let mut data = structure().data().clone();
    let mut bonds = pdbiox_core::BondTableBuilder::new();
    bonds.push(pdbiox_core::BondRecord {
        atom_a: pdbiox_core::AtomIndex::new(0),
        atom_b: pdbiox_core::AtomIndex::new(1),
        order: pdbiox_core::BondOrder::Aromatic,
        provenance: pdbiox_core::BondProvenance::ChemicalComponentDictionary,
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
        AnalysisPolicy::default().with_identifiers(pdbiox_core::contract::Namespace::Explicit);
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
        pdbiox_core::SEGMENT_ID_ANNOTATION,
        pdbiox_core::AtomAnnotation::Symbol(pdbiox_core::AnnotationColumn::from_values(vec![
            system, system, solvent,
        ])),
    );
    let _ = data.annotations.insert(
        pdbiox_core::PLDDT_ANNOTATION,
        pdbiox_core::AtomAnnotation::Real(pdbiox_core::AnnotationColumn::from_values(vec![
            90.0, 40.0, 80.0,
        ])),
    );
    let structure = Structure::new(data);
    for (source, expected) in [
        ("segid SYSTEM", vec![0, 1]),
        ("atom SYSTEM 42 CA", vec![1]),
        ("same segment as name CA", vec![0, 1]),
        ("plddt > 70", vec![0, 2]),
    ] {
        let query = match Query::compile(source) {
            Ok(query) => query,
            Err(findings) => panic!("compile failed: {findings:?}"),
        };
        let result = query.evaluate(&structure, &AnalysisPolicy::default(), &Groups::new(), None);
        let result = match result {
            Ok(result) => result,
            Err(findings) => panic!("evaluation failed: {findings:?}"),
        };
        assert_eq!(result.selection.iter().collect::<Vec<_>>(), expected);
    }
}
