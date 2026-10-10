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

fn options() -> HelicalOptions {
    HelicalOptions {
        frame_tolerance: 1e-9,
    }
}

fn turned(angle_degrees: f64, origin: [f64; 3]) -> BaseFrame {
    let (s, c) = angle_degrees.to_radians().sin_cos();
    BaseFrame {
        origin,
        x: [c, s, 0.0],
        y: [-s, c, 0.0],
        z: [0.0, 0.0, 1.0],
    }
}

#[test]
fn rise_and_twist_ladder_gives_exact_middle_frame_values() {
    let frames: Vec<BaseFrame> = (0..5)
        .map(|i| turned(36.0 * f64::from(i), [0.0, 0.0, 3.4 * f64::from(i)]))
        .collect();
    let steps = middle_frame_steps(&frames, options()).expect("valid frames");
    assert_eq!(steps.len(), 4);
    for step in steps {
        assert!((step.rise - 3.4).abs() < 1e-9, "{step:?}");
        assert!((step.twist_degrees - 36.0).abs() < 1e-9, "{step:?}");
        for zero in [step.shift, step.slide, step.tilt_degrees, step.roll_degrees] {
            assert!(zero.abs() < 1e-9, "{step:?}");
        }
    }
}

#[test]
fn middle_frame_slide_is_read_along_the_middle_y_axis() {
    // Twisting 90 degrees while moving along the average y direction: in the
    // first frame that move is split between x and y, in the middle frame it is
    // pure slide.
    let second = turned(90.0, [-1.0, 1.0, 0.0]);
    let step = middle_frame_step(identity([0.0; 3]), second, options()).expect("valid frames");
    let root = 2.0_f64.sqrt();
    assert!(step.shift.abs() < 1e-9, "{step:?}");
    assert!((step.slide - root).abs() < 1e-9, "{step:?}");
    assert!((step.twist_degrees - 90.0).abs() < 1e-9);
}

#[test]
fn pure_roll_and_tilt_are_bends_about_y_and_x() {
    let (s, c) = 10.0_f64.to_radians().sin_cos();
    let roll = BaseFrame {
        origin: [0.0; 3],
        x: [c, 0.0, -s],
        y: [0.0, 1.0, 0.0],
        z: [s, 0.0, c],
    };
    let step = middle_frame_step(identity([0.0; 3]), roll, options()).expect("valid frames");
    assert!((step.roll_degrees - 10.0).abs() < 1e-9, "{step:?}");
    assert!(step.tilt_degrees.abs() < 1e-9 && step.twist_degrees.abs() < 1e-9);
    let tilt = BaseFrame {
        origin: [0.0; 3],
        x: [1.0, 0.0, 0.0],
        y: [0.0, c, s],
        z: [0.0, -s, c],
    };
    let step = middle_frame_step(identity([0.0; 3]), tilt, options()).expect("valid frames");
    assert!((step.tilt_degrees - 10.0).abs() < 1e-9, "{step:?}");
    assert!(step.roll_degrees.abs() < 1e-9 && step.twist_degrees.abs() < 1e-9);
}
