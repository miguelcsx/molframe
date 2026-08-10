use super::INSTANCE_ID_ANNOTATION;
use crate::AssemblyExt;
use crate::view::tests::{ENTRY, attached};
use pdbiox_core::annotation::{AnnotationColumn, AtomAnnotation};
use pdbiox_core::bond::{BondRecord, BondTableBuilder};
use pdbiox_core::{AtomIndex, BondOrder, BondProvenance, Structure};

#[test]
fn materialization_copies_instances_and_preserves_their_identity() {
    let structure = with_bond_and_annotation(&attached(ENTRY));
    let view = match structure.assembly("1") {
        Ok(view) => view,
        Err(finding) => panic!("assembly failed: {finding}"),
    };
    let materialized = match view.materialize() {
        Ok(structure) => structure,
        Err(findings) => panic!("materialization failed: {findings:?}"),
    };

    assert_eq!(materialized.atom_count(), 2);
    assert_eq!(materialized.chain_count(), 2);
    assert_eq!(materialized.data().bonds.len(), 1);
    assert!(materialized.extensions().is_empty());
    let labels = materialized
        .data()
        .chains()
        .filter_map(pdbiox_core::structure::ChainRef::label)
        .collect::<Vec<_>>();
    assert_eq!(labels, ["A", "B"]);
    let Some(AtomAnnotation::Integer(instances)) =
        materialized.annotations().get(INSTANCE_ID_ANNOTATION)
    else {
        panic!("instance identity annotation absent")
    };
    assert_eq!(instances.values(), [0, 1]);
    let Some(AtomAnnotation::Integer(scores)) = materialized.annotations().get("score") else {
        panic!("source annotation absent")
    };
    assert_eq!(scores.values(), [7, 8]);
}

#[test]
fn repeated_source_chains_receive_collision_free_suffixes() {
    let repeated = ENTRY.replace("'(T)(R)'", "'(T,T)(R)'");
    let structure = attached(&repeated);
    let view = match structure.assembly("1") {
        Ok(view) => view,
        Err(finding) => panic!("assembly failed: {finding}"),
    };
    let materialized = match view.materialize() {
        Ok(structure) => structure,
        Err(findings) => panic!("materialization failed: {findings:?}"),
    };
    let labels = materialized
        .data()
        .chains()
        .filter_map(pdbiox_core::structure::ChainRef::label)
        .collect::<Vec<_>>();
    assert_eq!(labels, ["A", "B", "A-2", "B-2"]);
}

fn with_bond_and_annotation(structure: &Structure) -> Structure {
    let mut data = structure.data().clone();
    let mut bonds = BondTableBuilder::new();
    bonds.push(BondRecord {
        atom_a: AtomIndex::new(0),
        atom_b: AtomIndex::new(1),
        order: BondOrder::Single,
        provenance: BondProvenance::File,
    });
    data.bonds = bonds.finish();
    let Ok(scores) = AnnotationColumn::from_values(vec![7, 8]) else {
        panic!("two annotation values fit the column index domain");
    };
    let _ = data
        .annotations
        .insert("score", AtomAnnotation::Integer(scores));
    Structure::new(data)
}
