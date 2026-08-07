use super::*;

#[test]
fn validity_is_retained_without_option_erasing_its_meaning() {
    let column = AnnotationColumn::from_entries([
        (7_i64, Presence::Present),
        (0, Presence::Unknown),
        (0, Presence::Inapplicable),
    ]);
    assert_eq!(column.get(0), Some((7, Presence::Present)));
    assert_eq!(column.get(1), Some((0, Presence::Unknown)));
    assert_eq!(column.get(2), Some((0, Presence::Inapplicable)));
    assert_eq!(column.get(3), None);
}

#[test]
fn annotation_names_iterate_deterministically() {
    let mut annotations = AtomAnnotations::default();
    let _ = annotations.insert(
        "zeta",
        AtomAnnotation::Boolean(AnnotationColumn::from_values(vec![true])),
    );
    let _ = annotations.insert(
        "alpha",
        AtomAnnotation::Real(AnnotationColumn::from_values(vec![1.0])),
    );
    assert_eq!(
        annotations.iter().map(|(name, _)| name).collect::<Vec<_>>(),
        vec!["alpha", "zeta"]
    );
}
