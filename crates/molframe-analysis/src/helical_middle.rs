//! Middle-frame base-pair step parameters.
//!
//! The step parameters most nucleic-acid tools report (shift, slide, rise, tilt,
//! roll, twist) are not read in either base-pair frame. They are read in a
//! middle frame that sits halfway between the two, which makes them symmetric
//! under swapping the two frames' order and independent of which frame is
//! called first. The rotation between the frames is split into a twist about the
//! middle frame's z axis and a single bend about a hinge axis in its xy plane;
//! the bend's components along x and y are the tilt and the roll.

use super::{
    BaseFrame, HelicalError, HelicalOptions, dot, subtract, validate_frame, validate_options,
};

/// Six base-pair step parameters measured in the middle frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MiddleFrameStep {
    /// Translation along the middle frame's x axis.
    pub shift: f64,
    /// Translation along the middle frame's y axis.
    pub slide: f64,
    /// Translation along the middle frame's z axis.
    pub rise: f64,
    /// Bend about the middle frame's x axis, in degrees.
    pub tilt_degrees: f64,
    /// Bend about the middle frame's y axis, in degrees.
    pub roll_degrees: f64,
    /// Twist about the middle frame's z axis, in degrees.
    pub twist_degrees: f64,
}

/// Computes the middle-frame step parameters between two consecutive frames.
///
/// A ladder with a pure rise and a pure twist gives exactly those two values and
/// zeros elsewhere. A bend about the y axis is a pure roll and a bend about the
/// x axis a pure tilt.
///
/// # Errors
///
/// Returns an error for invalid tolerance or frame geometry.
pub fn middle_frame_step(
    first: BaseFrame,
    second: BaseFrame,
    options: HelicalOptions,
) -> Result<MiddleFrameStep, HelicalError> {
    validate_options(options)?;
    validate_frame(first, options.frame_tolerance)?;
    validate_frame(second, options.frame_tolerance)?;
    let relative = super::relative_rotation(first, second);
    let (twist, bend, hinge) = decompose(relative, options.frame_tolerance);
    // The middle frame is the first frame turned through half the twist and half
    // the bend, which is also what the second frame reaches by turning back.
    let alpha = twist / 2.0 + hinge;
    let turn = multiply(
        multiply(rotation_z(alpha), rotation_y(bend / 2.0)),
        rotation_z(-hinge),
    );
    let axes = [first.x, first.y, first.z];
    let middle: [[f64; 3]; 3] = std::array::from_fn(|column| {
        std::array::from_fn(|row| (0..3).map(|k| axes[k][row] * turn[k][column]).sum())
    });
    let displacement = subtract(second.origin, first.origin);
    Ok(MiddleFrameStep {
        shift: dot(displacement, middle[0]),
        slide: dot(displacement, middle[1]),
        rise: dot(displacement, middle[2]),
        tilt_degrees: (-bend * hinge.sin()).to_degrees(),
        roll_degrees: (bend * hinge.cos()).to_degrees(),
        twist_degrees: twist.to_degrees(),
    })
}

/// Computes middle-frame step parameters for every consecutive frame pair.
///
/// # Errors
///
/// Returns the first invalid frame or tolerance error.
pub fn middle_frame_steps(
    frames: &[BaseFrame],
    options: HelicalOptions,
) -> Result<Vec<MiddleFrameStep>, HelicalError> {
    frames
        .windows(2)
        .map(|pair| middle_frame_step(pair[0], pair[1], options))
        .collect()
}

/// Splits a rotation into `Rz(a) · Rhinge(bend) · Rz(a)` and returns
/// `(twist = 2a, bend, hinge angle from the y axis)`, the hinge axis being
/// `(-sin, cos, 0)` at that angle.
///
/// That product equals the Euler form `Rz(a + hinge) Ry(bend) Rz(a - hinge)`, so
/// the three numbers come from a z–y–z Euler decomposition.
fn decompose(matrix: [[f64; 3]; 3], tolerance: f64) -> (f64, f64, f64) {
    let bend = matrix[2][2].clamp(-1.0, 1.0).acos();
    if bend.sin().abs() <= tolerance {
        // No bend: only the total z rotation is defined.
        let twist = matrix[1][0].atan2(matrix[0][0]);
        return (twist, 0.0, 0.0);
    }
    let alpha = matrix[1][2].atan2(matrix[0][2]);
    let gamma = matrix[2][1].atan2(-matrix[2][0]);
    let mut half_twist = f64::midpoint(alpha, gamma);
    let mut hinge = (alpha - gamma) / 2.0;
    // Keep the twist inside (-180, 180]; a whole turn of the twist moves the
    // hinge by a half turn.
    if half_twist > std::f64::consts::FRAC_PI_2 {
        half_twist -= std::f64::consts::PI;
        hinge -= std::f64::consts::PI;
    } else if half_twist <= -std::f64::consts::FRAC_PI_2 {
        half_twist += std::f64::consts::PI;
        hinge += std::f64::consts::PI;
    }
    (2.0 * half_twist, bend, hinge)
}

fn rotation_z(angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]]
}

fn rotation_y(angle: f64) -> [[f64; 3]; 3] {
    let (s, c) = angle.sin_cos();
    [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]]
}

fn multiply(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|row| {
        std::array::from_fn(|column| (0..3).map(|k| left[row][k] * right[k][column]).sum())
    })
}
