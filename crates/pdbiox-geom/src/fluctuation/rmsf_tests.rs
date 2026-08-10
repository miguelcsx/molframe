use super::{FluctuationError, rmsf};

#[test]
fn an_empty_ensemble_has_no_mean_to_measure_against() {
    let frames: [&[[f32; 3]]; 0] = [];
    let Err(error) = rmsf(&frames) else {
        panic!("no frames should be rejected");
    };
    assert_eq!(error, FluctuationError::NoFrames);
}

#[test]
fn frames_that_disagree_on_atom_count_are_rejected() {
    let a: [[f32; 3]; 2] = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let b: [[f32; 3]; 1] = [[0.0, 0.0, 0.0]];
    let Err(error) = rmsf(&[&a, &b]) else {
        panic!("ragged frames should be rejected");
    };
    assert_eq!(error, FluctuationError::RaggedFrames);
}

#[test]
fn a_single_frame_cannot_depart_from_its_own_mean() {
    let only: [[f32; 3]; 2] = [[3.0, -1.0, 4.0], [0.0, 2.0, 0.0]];
    let Ok(values) = rmsf(&[&only]) else {
        panic!("a single frame is valid");
    };
    assert_eq!(values.len(), 2);
    assert!(values.iter().all(|value| value.abs() < 1e-12));
}

#[test]
fn two_frames_split_symmetrically_give_half_the_separation() {
    let a: [[f32; 3]; 1] = [[0.0, 0.0, 0.0]];
    let b: [[f32; 3]; 1] = [[2.0, 0.0, 0.0]];
    let Ok(values) = rmsf(&[&a, &b]) else {
        panic!("two frames are valid");
    };
    assert!((values[0] - 1.0).abs() < 1e-9);
}

#[test]
fn identical_frames_have_zero_fluctuation() {
    let frame: [[f32; 3]; 3] = [[1.0, 2.0, 3.0], [-4.0, 5.0, 6.0], [7.0, -8.0, 9.0]];
    let Ok(values) = rmsf(&[&frame, &frame, &frame, &frame]) else {
        panic!("repeated frames are valid");
    };
    assert!(values.iter().all(|value| value.abs() < 1e-9));
}

#[test]
fn fluctuation_is_independent_of_frame_order() {
    let a: [[f32; 3]; 1] = [[0.0, 0.0, 0.0]];
    let b: [[f32; 3]; 1] = [[3.0, 0.0, 0.0]];
    let c: [[f32; 3]; 1] = [[0.0, 6.0, 0.0]];
    let Ok(forward) = rmsf(&[&a, &b, &c]) else {
        panic!("valid");
    };
    let Ok(reversed) = rmsf(&[&c, &b, &a]) else {
        panic!("valid");
    };
    assert!((forward[0] - reversed[0]).abs() < 1e-9);
}
