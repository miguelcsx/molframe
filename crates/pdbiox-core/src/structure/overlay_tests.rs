use super::*;
use crate::annotation::{AnnotationColumn, AtomAnnotation};
use crate::column::Presence;
use crate::structure::fixture;

#[test]
fn renaming_a_chain_publishes_a_new_snapshot_and_shares_atom_storage() {
    let original = crate::structure::fixture::sample();
    let original_chunks = original.data().chunks.clone();
    let original_positions = original.positions().as_ptr();
    let mut editor = original.edit();
    if let Err(finding) = editor.rename_chain(ChainIndex::new(0), "RENAMED") {
        panic!("rename failed: {finding}")
    }
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };

    assert_eq!(
        original
            .data()
            .chain_named("A")
            .and_then(crate::structure::ChainRef::label),
        Some("A")
    );
    assert_eq!(
        edited
            .data()
            .chain_named("RENAMED")
            .and_then(crate::structure::ChainRef::label),
        Some("RENAMED")
    );
    assert!(std::sync::Arc::ptr_eq(
        &original_chunks,
        &edited.data().chunks
    ));
    assert!(std::ptr::eq(
        original_positions,
        edited.positions().as_ptr()
    ));
}

#[test]
fn an_invalid_chain_edit_is_a_diagnostic_not_a_partial_mutation() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit();
    let result = editor.rename_chain(ChainIndex::new(99), "B");
    assert_eq!(
        result.err().map(|finding| finding.code()),
        Some(Code::E6006)
    );
}

#[test]
fn a_transform_is_transactional_and_advances_generation_once() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit();
    let selected = AtomSelection::from_sorted(vec![0, 23]);
    if let Err(finding) = editor.transform(&selected, |position| {
        [position[0] + 2.0, position[1] - 1.0, position[2] + 4.0]
    }) {
        panic!("transform failed: {finding}")
    }
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };

    assert_position(original.positions()[0], [0.0, 0.0, 0.0]);
    assert_position(edited.positions()[0], [2.0, -1.0, 4.0]);
    assert_position(edited.positions()[23], [3.0, 1.0, 7.0]);
    assert_eq!(Some(edited.generation()), original.generation().next());
}

#[test]
fn deleting_atoms_compacts_rows_residue_ranges_coordinates_and_bonds() {
    let fixture = crate::structure::fixture::sample();
    let original = with_annotation(&with_bonds(&fixture));
    let mut editor = original.edit();
    let deleted = AtomSelection::from_sorted(vec![1, 5, 23]);
    if let Err(finding) = editor.delete_atoms(&deleted) {
        panic!("delete failed: {finding}")
    }
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };

    assert_eq!(original.atom_count(), 24);
    assert_eq!(edited.atom_count(), 21);
    assert_position(edited.positions()[1], original.positions()[2]);
    assert_eq!(
        edited.data().topology.residues.atoms(ResidueIndex::new(0)),
        Some(0..3)
    );
    assert_eq!(
        edited.data().topology.residues.atoms(ResidueIndex::new(1)),
        Some(3..6)
    );
    assert_eq!(
        edited.data().topology.residues.atoms(ResidueIndex::new(5)),
        Some(18..21)
    );
    assert_eq!(edited.data().bonds.len(), 1);
    let Some(bond) = edited.data().bonds.get(crate::BondIndex::new(0)) else {
        panic!("remaining bond was lost")
    };
    assert_eq!((bond.atom_a.get(), bond.atom_b.get()), (1, 2));
    let Some(AtomAnnotation::Integer(column)) = edited.annotations().get("score") else {
        panic!("annotation was lost")
    };
    assert_eq!(column.len(), 21);
    assert_eq!(column.get(1), Some((2, Presence::Present)));
    assert_eq!(Some(edited.generation()), original.generation().next());
}

#[test]
fn annotation_edits_are_length_checked_and_copy_on_write() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit();
    let short = AtomAnnotation::Boolean(
        AnnotationColumn::from_values(vec![true]).expect("small annotation column"),
    );
    assert_eq!(
        editor
            .set_annotation("flag", short)
            .err()
            .map(|finding| finding.code()),
        Some(Code::E3011)
    );
    let values = (0..original.atom_count()).map(i64::from).collect();
    if let Err(finding) = editor.set_annotation(
        "score",
        AtomAnnotation::Integer(AnnotationColumn::from_values(values).expect("small column")),
    ) {
        panic!("annotation failed: {finding}")
    }
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };
    assert!(original.annotations().is_empty());
    assert!(edited.annotations().get("score").is_some());
    assert!(std::ptr::eq(
        original.positions().as_ptr(),
        edited.positions().as_ptr()
    ));
}

#[test]
fn a_failed_structural_operation_leaves_the_editors_snapshot_unchanged() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit();
    let invalid = AtomSelection::from_sorted(vec![original.atom_count()]);
    assert_eq!(
        editor
            .delete_atoms(&invalid)
            .err()
            .map(|finding| finding.code()),
        Some(Code::E6009)
    );
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };
    assert_eq!(edited.atom_count(), original.atom_count());
    assert!(std::ptr::eq(
        edited.positions().as_ptr(),
        original.positions().as_ptr()
    ));
}

#[test]
fn topology_edits_require_explicit_extension_removal() {
    let original = crate::structure::fixture::sample()
        .with_extension("example.topology.v1", String::from("derived"));
    let mut editor = original.edit();
    assert_eq!(
        editor
            .rename_chain(ChainIndex::new(0), "RENAMED")
            .err()
            .map(|finding| finding.code()),
        Some(Code::E3014)
    );

    editor.clear_extensions();
    if let Err(finding) = editor.rename_chain(ChainIndex::new(0), "RENAMED") {
        panic!("rename after extension removal failed: {finding}")
    }
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };
    assert!(edited.extensions().is_empty());
}

#[test]
fn coordinate_edits_preserve_domain_extensions() {
    let original = crate::structure::fixture::sample().with_extension("example.v1", 42_u32);
    let mut editor = original.edit();
    if let Err(finding) =
        editor.transform(&AtomSelection::from_sorted(vec![0]), |position| position)
    {
        panic!("transform failed: {finding}")
    }
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };
    assert_eq!(edited.extensions().get::<u32>("example.v1"), Some(&42));
}

fn with_bonds(structure: &Structure) -> Structure {
    let mut data = structure.data().clone();
    let mut bonds = BondTableBuilder::new();
    for (atom_a, atom_b) in [(0, 1), (2, 3)] {
        bonds.push(BondRecord {
            atom_a: AtomIndex::new(atom_a),
            atom_b: AtomIndex::new(atom_b),
            order: crate::BondOrder::Single,
            provenance: crate::BondProvenance::User,
        });
    }
    data.bonds = bonds.finish();
    Structure::new(data)
}

fn with_annotation(structure: &Structure) -> Structure {
    let mut editor = structure.edit();
    let values = (0..structure.atom_count()).map(i64::from).collect();
    if editor
        .set_annotation(
            "score",
            AtomAnnotation::Integer(AnnotationColumn::from_values(values).expect("small column")),
        )
        .is_err()
    {
        return structure.clone();
    }
    match editor.commit() {
        Ok(edited) => edited,
        Err(_) => structure.clone(),
    }
}

#[test]
fn materialization_compacts_hierarchy_bonds_annotations_and_extensions() {
    let source = with_annotation(&with_bonds(&fixture::sample()))
        .with_extension("selection.invalidated.v1", 7u32);
    let selection = AtomSelection::from_sorted(vec![0, 1]);
    let selected = match source.materialize(&selection) {
        Ok(structure) => structure,
        Err(findings) => panic!("materialization failed: {findings:?}"),
    };

    assert_eq!(selected.atom_count(), 2);
    assert_eq!(selected.residue_count(), 1);
    assert_eq!(selected.chain_count(), 1);
    assert_eq!(selected.entity_count(), 1);
    assert_eq!(selected.data().bonds.len(), 1);
    assert_eq!(
        selected.annotations().get("score").map(AtomAnnotation::len),
        Some(2)
    );
    assert!(selected.extensions().is_empty());
    assert_eq!(
        selected
            .atom(AtomIndex::new(0))
            .and_then(crate::structure::AtomRef::name),
        Some("N")
    );
    assert_eq!(
        selected
            .atom(AtomIndex::new(1))
            .and_then(crate::structure::AtomRef::name),
        Some("CA")
    );
}

#[test]
fn materialization_rejects_atoms_outside_the_structure() {
    let source = fixture::sample();
    let result = source.materialize(&AtomSelection::from_sorted(vec![source.atom_count()]));
    assert!(
        matches!(result, Err(findings) if findings.iter().any(|finding| finding.code() == Code::E6009))
    );
}

fn assert_position(actual: [f32; 3], expected: [f32; 3]) {
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!((actual - expected).abs() < f32::EPSILON);
    }
}
