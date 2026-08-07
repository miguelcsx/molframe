use crate::diagnostic::Code;
use crate::index::ModelIndex;

#[test]
fn committing_coordinates_creates_a_new_generation_and_keeps_the_original() {
    let original = crate::structure::fixture::sample();
    let before = original.positions()[0];
    let mut editor = original.edit_coordinates();
    let Some(positions) = editor.positions_mut(ModelIndex::new(0)) else {
        panic!("first model missing")
    };
    positions[0] = [9.0, 8.0, 7.0];
    let edited = match editor.commit() {
        Ok(structure) => structure,
        Err(findings) => panic!("commit failed: {findings:?}"),
    };

    assert!(
        original.positions()[0]
            .iter()
            .zip(before)
            .all(|(actual, expected)| (actual - expected).abs() < f32::EPSILON)
    );
    assert!(
        edited.positions()[0]
            .iter()
            .zip([9.0, 8.0, 7.0])
            .all(|(actual, expected)| (actual - expected).abs() < f32::EPSILON)
    );
    assert_eq!(original.generation().get(), 0);
    assert_eq!(edited.generation().get(), 1);
}

#[test]
fn a_non_finite_coordinate_aborts_the_transaction() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit_coordinates();
    let Some(positions) = editor.positions_mut(ModelIndex::new(0)) else {
        panic!("first model missing")
    };
    positions[0][0] = f32::NAN;
    let Err(findings) = editor.commit() else {
        panic!("non-finite coordinate committed")
    };
    assert!(findings.iter().any(|finding| finding.code() == Code::E3008));
    assert!(original.positions()[0][0].is_finite());
}

#[test]
fn checked_coordinate_access_uses_the_diagnostic_registry() {
    let original = crate::structure::fixture::sample();
    let mut editor = original.edit_coordinates();
    let result = editor.try_positions_mut(ModelIndex::new(1));
    assert_eq!(
        result.err().map(|finding| finding.code()),
        Some(Code::E6003)
    );
}
