use super::*;
use crate::annotation::{AnnotationColumn, AtomAnnotation};
use crate::column::Presence;
use crate::structure::fixture;

fn annotated(label: &str) -> Structure {
    let source = fixture::sample();
    let mut data = source.data().clone();
    let symbol = match data.dictionary.intern(label) {
        Ok(symbol) => symbol,
        Err(error) => panic!("fixture dictionary failed: {error}"),
    };
    let _ = data.annotations.insert(
        "source",
        AtomAnnotation::Symbol(AnnotationColumn::from_values(vec![
            symbol;
            source.atom_count()
                as usize
        ])),
    );
    Structure::new(data)
}

#[test]
fn merge_remaps_hierarchy_coordinates_dictionary_and_annotations() {
    let left = annotated("left-system");
    let right = annotated("right-system");
    let merged = match Structure::merge(&[left, right]) {
        Ok(merged) => merged,
        Err(findings) => panic!("merge failed: {findings:?}"),
    };
    assert_eq!(merged.atom_count(), 48);
    assert_eq!(merged.residue_count(), 12);
    assert_eq!(merged.chain_count(), 4);
    assert_eq!(merged.entity_count(), 2);
    assert!(
        merged.positions()[24]
            .iter()
            .all(|coordinate| coordinate.abs() < f32::EPSILON)
    );

    let Some(AtomAnnotation::Symbol(column)) = merged.annotations().get("source") else {
        panic!("merged symbol annotation absent")
    };
    let Some((left, _)) = column.get(0) else {
        panic!("left value absent")
    };
    let Some((right, _)) = column.get(24) else {
        panic!("right value absent")
    };
    assert_eq!(merged.resolve(left), Some("left-system"));
    assert_eq!(merged.resolve(right), Some("right-system"));
    assert!(validate(merged.data()).is_empty());
}

#[test]
fn an_annotation_absent_from_one_input_is_inapplicable_there() {
    let left = annotated("left-system");
    let right = fixture::sample();
    let merged = match Structure::merge(&[left, right]) {
        Ok(merged) => merged,
        Err(findings) => panic!("merge failed: {findings:?}"),
    };
    let Some(AtomAnnotation::Symbol(column)) = merged.annotations().get("source") else {
        panic!("merged annotation absent")
    };
    assert_eq!(column.presence(23), Presence::Present);
    assert_eq!(column.presence(24), Presence::Inapplicable);
}

#[test]
fn unequal_frame_axes_are_refused_instead_of_padded() {
    let left = fixture::sample();
    let mut data = left.data().clone();
    let first = data.coords.block(ModelIndex::new(0)).cloned();
    let Some(first) = first else {
        panic!("fixture frame absent")
    };
    data.coords = CoordinateStore::Dense {
        frames: vec![first.clone(), first],
    };
    let result = Structure::merge(&[Structure::new(data), fixture::sample()]);
    let Err(findings) = result else {
        panic!("unequal axes accepted")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E3012));
}

#[test]
fn multi_source_merge_requires_explicit_extension_removal() {
    let extended = fixture::sample().with_extension("example.topology.v1", 1_u8);
    let result = Structure::merge(&[extended.clone(), fixture::sample()]);
    let Err(findings) = result else {
        panic!("merge silently discarded an extension")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E3014));

    let merged = Structure::merge(&[extended.without_extensions(), fixture::sample()]);
    assert!(merged.is_ok());
}

#[test]
fn single_source_merge_preserves_extensions() {
    let extended = fixture::sample().with_extension("example.v1", String::from("kept"));
    let merged = match Structure::merge(std::slice::from_ref(&extended)) {
        Ok(merged) => merged,
        Err(findings) => panic!("single-source merge failed: {findings:?}"),
    };
    assert_eq!(
        merged
            .extensions()
            .get::<String>("example.v1")
            .map(String::as_str),
        Some("kept")
    );
}
