use super::*;
use molframe_core::ExecutionContext;

#[test]
fn identical_contact_areas_score_one() {
    let contacts = [
        ContactArea {
            first: 1,
            second: 2,
            area: 10.0,
        },
        ContactArea {
            first: 2,
            second: 3,
            area: 5.0,
        },
    ];
    let result =
        cad_score(&contacts, &contacts).unwrap_or_else(|error| panic!("CAD failed: {error}"));
    assert!((result.score - 1.0).abs() < f64::EPSILON);
    assert!(
        result
            .local
            .iter()
            .all(|local| (local.score - 1.0).abs() < f64::EPSILON)
    );
}

#[test]
fn missing_and_oversized_contacts_are_capped_per_reference_pair() {
    let reference = [
        ContactArea {
            first: 1,
            second: 2,
            area: 10.0,
        },
        ContactArea {
            first: 2,
            second: 3,
            area: 10.0,
        },
    ];
    let model = [ContactArea {
        first: 2,
        second: 1,
        area: 30.0,
    }];
    let result =
        cad_score(&reference, &model).unwrap_or_else(|error| panic!("CAD failed: {error}"));
    assert!((result.score - 0.0).abs() < f64::EPSILON);
    assert_eq!(result.contacts.len(), 2);
    assert!((result.local[1].score - 0.0).abs() < f64::EPSILON);
}

#[test]
fn duplicate_orientations_are_coalesced_deterministically() {
    let reference = [
        ContactArea {
            first: 4,
            second: 7,
            area: 2.0,
        },
        ContactArea {
            first: 7,
            second: 4,
            area: 3.0,
        },
    ];
    let model = [ContactArea {
        first: 4,
        second: 7,
        area: 4.0,
    }];
    let result =
        cad_score(&reference, &model).unwrap_or_else(|error| panic!("CAD failed: {error}"));
    assert!((result.score - 0.8).abs() < f64::EPSILON);
    assert!((result.contacts[0].reference_area - 5.0).abs() < f64::EPSILON);
}

#[test]
fn solvent_excluded_atom_patches_construct_residue_contact_areas() {
    let result = cad_contact_areas(
        &[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        &[1.5, 1.5],
        &[4, 7],
        0.0,
        50.0,
        &ExecutionContext::default(),
    )
    .unwrap_or_else(|error| panic!("CAD construction failed: {error}"));
    assert_eq!(result.len(), 1);
    assert_eq!((result[0].first, result[0].second), (4, 7));
    assert!(result[0].area > 0.0);
}
