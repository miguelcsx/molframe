use super::*;

fn identity(origin: [f64; 3]) -> BaseFrame {
    BaseFrame {
        origin,
        x: [1.0, 0.0, 0.0],
        y: [0.0, 1.0, 0.0],
        z: [0.0, 0.0, 1.0],
    }
}

#[test]
fn translation_and_twist_follow_the_explicit_first_frame() {
    let second = BaseFrame {
        origin: [1.0, 2.0, 3.0],
        x: [0.0, 1.0, 0.0],
        y: [-1.0, 0.0, 0.0],
        z: [0.0, 0.0, 1.0],
    };
    let result = helical_parameters(
        identity([0.0; 3]),
        second,
        HelicalOptions {
            frame_tolerance: 1e-12,
        },
    )
    .expect("orthonormal frames");
    assert!((result.x_displacement - 1.0).abs() < 1e-12);
    assert!((result.y_displacement - 2.0).abs() < 1e-12);
    assert!((result.z_displacement - 3.0).abs() < 1e-12);
    assert!((result.z_rotation_degrees - 90.0).abs() < 1e-10);
}

#[test]
fn non_orthonormal_frames_are_refused() {
    let mut invalid = identity([0.0; 3]);
    invalid.x = [2.0, 0.0, 0.0];
    assert_eq!(
        helical_parameters(
            invalid,
            identity([0.0; 3]),
            HelicalOptions {
                frame_tolerance: 1e-12,
            },
        ),
        Err(HelicalError::InvalidFrame)
    );
}
