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
        AtomAnnotation::Symbol(
            AnnotationColumn::from_values(vec![symbol; source.atom_count() as usize])
                .expect("small annotation column"),
        ),
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

#[test]
fn merge_offsets_anisotropy_rows_and_propagates_availability() {
    let left = with_anisotropy(&fixture::sample(), true);
    let right = with_anisotropy(&fixture::sample(), true);
    let merged = match Structure::merge(&[left, right]) {
        Ok(merged) => merged,
        Err(findings) => panic!("merge failed: {findings:?}"),
    };
    let anisotropy = &merged.data().anisotropy;
    assert!(anisotropy.is_available());
    assert_eq!(anisotropy.len(), 4);
    let Some(record) = anisotropy.get(crate::AnisotropyIndex::new(2)) else {
        panic!("second input's first ellipsoid absent")
    };
    // The second input's atoms moved by the first input's atom count.
    assert_eq!(record.atom, crate::AtomIndex::new(24));
    assert_tensor(record.u, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    let Some(later) = anisotropy.get(crate::AnisotropyIndex::new(3)) else {
        panic!("second input's last ellipsoid absent");
    };
    assert_eq!(later.atom, crate::AtomIndex::new(47));
}

#[test]
fn an_unavailable_anisotropy_source_makes_the_merge_unavailable() {
    let left = with_anisotropy(&fixture::sample(), true);
    let right = with_anisotropy(&fixture::sample(), false);
    let merged = match Structure::merge(&[left, right]) {
        Ok(merged) => merged,
        Err(findings) => panic!("merge failed: {findings:?}"),
    };
    let anisotropy = &merged.data().anisotropy;
    assert!(!anisotropy.is_available());
    // The rows from the available source survive; availability is about the
    // whole set being resolved, not about rows existing.
    assert_eq!(anisotropy.len(), 4);
}

fn with_anisotropy(structure: &Structure, available: bool) -> Structure {
    let mut data = structure.data().clone();
    let mut builder = crate::anisotropy::AnisotropyTableBuilder::new();
    builder.push(crate::AnisotropicDisplacement {
        atom: crate::AtomIndex::new(0),
        u: [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    });
    builder.push(crate::AnisotropicDisplacement {
        atom: crate::AtomIndex::new(structure.atom_count() - 1),
        u: [7.0; 6],
    });
    data.anisotropy = builder.finish_with_availability(available);
    Structure::new(data)
}

/// Asserts two six-component tensors agree componentwise.
fn assert_tensor(actual: [f32; 6], expected: [f32; 6]) {
    for (actual, expected) in actual.iter().zip(expected) {
        assert!((actual - expected).abs() <= f32::EPSILON);
    }
}
