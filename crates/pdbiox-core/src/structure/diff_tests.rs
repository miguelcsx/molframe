use super::*;
use crate::coords::CoordinateBlock;
use crate::structure::{CoordinateStore, fixture};

fn options(tolerance: f32) -> StructureDifferenceOptions {
    StructureDifferenceOptions {
        coordinate_tolerance: tolerance,
    }
}

#[test]
fn identical_snapshots_have_no_semantic_difference() {
    let structure = fixture::sample();
    let result = structure_difference(&structure, &structure, options(0.0));
    assert!(matches!(result, Ok(difference) if difference.is_empty()));
}

#[test]
fn metadata_and_coordinates_are_reported_independently() {
    let left = fixture::sample();
    let mut data = left.data().clone();
    data.entry.title = Some("changed".into());
    let mut positions = left.positions().to_vec();
    positions[0][0] += 0.25;
    data.coords = CoordinateStore::Single(positions.into_iter().collect::<CoordinateBlock>());
    let right = Structure::new(data);

    let result = structure_difference(&left, &right, options(0.1));
    let Ok(difference) = result else {
        panic!("valid comparison failed")
    };
    assert!(difference.metadata.title.is_some());
    assert_eq!(difference.changed_positions, 1);
    assert_eq!(difference.maximum_displacement, Some(0.25));
    assert_eq!(difference.changed_atoms, 0);
}

#[test]
fn coordinate_tolerance_must_be_explicitly_valid() {
    let structure = fixture::sample();
    assert_eq!(
        structure_difference(&structure, &structure, options(f32::NAN)),
        Err(DifferenceError::InvalidCoordinateTolerance)
    );
    assert_eq!(
        structure_difference(&structure, &structure, options(-0.1)),
        Err(DifferenceError::InvalidCoordinateTolerance)
    );
}
