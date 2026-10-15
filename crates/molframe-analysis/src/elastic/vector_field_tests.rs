use super::{
    StreamlineDirection, StreamlineOptions, VectorFieldError, VectorFieldGrid,
    integrate_streamlines,
};

fn constant_field() -> VectorFieldGrid {
    let result = VectorFieldGrid::new([0.0; 3], [1.0; 3], [4, 4, 4], vec![[1.0, 0.0, 0.0]; 64]);
    let Ok(field) = result else {
        panic!("constant field should be valid")
    };
    field
}

#[test]
fn trilinear_sampling_preserves_a_constant_field() {
    assert_eq!(
        constant_field().sample([1.25, 1.5, 2.0]),
        Some([1.0, 0.0, 0.0])
    );
}

#[test]
fn bidirectional_integration_joins_once_at_the_seed() {
    let options = StreamlineOptions {
        step_size: 0.5,
        max_steps: 2,
        max_length: 10.0,
        min_speed: 0.0,
        direction: StreamlineDirection::Both,
    };
    let result = integrate_streamlines(&constant_field(), &[[1.5, 1.5, 1.5]], options);
    let Ok(lines) = result else {
        panic!("streamline should integrate")
    };
    assert_eq!(lines[0].len(), 5);
    assert!(
        lines[0][2]
            .iter()
            .all(|value| (*value - 1.5).abs() < 1.0e-6)
    );
    assert!(lines[0].windows(2).all(|pair| pair[0][0] < pair[1][0]));
}

#[test]
fn malformed_storage_is_rejected_before_sampling() {
    let result = VectorFieldGrid::new([0.0; 3], [1.0; 3], [2, 2, 2], vec![[0.0; 3]; 7]);
    assert_eq!(
        result,
        Err(VectorFieldError::LengthMismatch {
            expected: 8,
            actual: 7,
        })
    );
}

#[test]
fn integration_is_byte_deterministic() {
    let options = StreamlineOptions::default();
    let seeds = [[0.5, 0.5, 0.5], [1.0, 2.0, 1.0]];
    let left = integrate_streamlines(&constant_field(), &seeds, options);
    let right = integrate_streamlines(&constant_field(), &seeds, options);
    assert_eq!(left, right);
}
