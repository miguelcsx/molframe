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

#[test]
fn vector_blocks_and_tail_match_the_scalar_welford_order_bit_for_bit() {
    let frames = (0..4_i16)
        .map(|frame| {
            (0..7_i16)
                .map(|atom| {
                    let value = f32::from(frame * 7 + atom) * 0.125;
                    [value, value.mul_add(-0.5, 1.0), value * value]
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let views = frames.iter().map(Vec::as_slice).collect::<Vec<_>>();
    let actual = rmsf(&views).expect("SIMD RMSF");
    let mut mean = [[0.0_f64; 3]; 7];
    let mut square = [0.0_f64; 7];
    for (frame_index, frame) in frames.iter().enumerate() {
        let inverse = f64::from(u32::try_from(frame_index + 1).expect("small frame count")).recip();
        for (atom, point) in frame.iter().enumerate() {
            let point = point.map(f64::from);
            let before = core::array::from_fn::<_, 3, _>(|axis| point[axis] - mean[atom][axis]);
            for axis in 0..3 {
                mean[atom][axis] += before[axis] * inverse;
            }
            let after = core::array::from_fn::<_, 3, _>(|axis| point[axis] - mean[atom][axis]);
            square[atom] += before[0] * after[0] + before[1] * after[1] + before[2] * after[2];
        }
    }
    for (actual, square) in actual.iter().zip(square) {
        let expected = (square / 4.0).sqrt();
        assert_eq!(actual.to_bits(), expected.to_bits());
    }
}
