//! Calibration-free nucleic parameters from explicit orthonormal base frames.

/// A caller-defined right-handed base frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaseFrame {
    /// Frame origin.
    pub origin: [f64; 3],
    /// Unit x axis.
    pub x: [f64; 3],
    /// Unit y axis.
    pub y: [f64; 3],
    /// Unit z axis.
    pub z: [f64; 3],
}

/// Explicit numerical tolerance for frame validation and zero rotations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelicalOptions {
    /// Maximum absolute orthonormality/right-handedness residual.
    pub frame_tolerance: f64,
}

/// Six rigid-body parameters in the first frame's coordinate system.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelicalParameters {
    /// Translation along the first frame's x axis.
    pub x_displacement: f64,
    /// Translation along the first frame's y axis.
    pub y_displacement: f64,
    /// Translation along the first frame's z axis.
    pub z_displacement: f64,
    /// Rotation-vector x component in degrees.
    pub x_rotation_degrees: f64,
    /// Rotation-vector y component in degrees.
    pub y_rotation_degrees: f64,
    /// Rotation-vector z component in degrees.
    pub z_rotation_degrees: f64,
}

/// Invalid frame geometry or tolerance.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum HelicalError {
    /// Frame tolerance must be positive and finite.
    #[error("frame tolerance must be positive and finite")]
    InvalidTolerance,
    /// A frame is non-finite, non-orthonormal or not right-handed.
    #[error("base frame must be finite, orthonormal and right-handed")]
    InvalidFrame,
}

/// Computes a calibration-free relative rigid transform between two base frames.
///
/// For a base pair the fields correspond to shear, stretch, stagger, buckle,
/// propeller and opening. For consecutive base-pair frames they correspond to
/// shift, slide, rise, tilt, roll and twist. The rotation values are components
/// of the exact SO(3) logarithm, avoiding empirical conformer tables.
///
/// # Errors
///
/// Returns an error for invalid tolerance or frame geometry.
pub fn helical_parameters(
    first: BaseFrame,
    second: BaseFrame,
    options: HelicalOptions,
) -> Result<HelicalParameters, HelicalError> {
    validate_options(options)?;
    validate_frame(first, options.frame_tolerance)?;
    validate_frame(second, options.frame_tolerance)?;
    let displacement = subtract(second.origin, first.origin);
    let translation = [
        dot(displacement, first.x),
        dot(displacement, first.y),
        dot(displacement, first.z),
    ];
    let rotation = relative_rotation(first, second);
    let vector = rotation_vector(rotation, options.frame_tolerance);
    Ok(HelicalParameters {
        x_displacement: translation[0],
        y_displacement: translation[1],
        z_displacement: translation[2],
        x_rotation_degrees: vector[0].to_degrees(),
        y_rotation_degrees: vector[1].to_degrees(),
        z_rotation_degrees: vector[2].to_degrees(),
    })
}

/// Computes parameters for every consecutive frame pair.
///
/// # Errors
///
/// Returns the first invalid frame/tolerance error.
pub fn helical_steps(
    frames: &[BaseFrame],
    options: HelicalOptions,
) -> Result<Vec<HelicalParameters>, HelicalError> {
    frames
        .windows(2)
        .map(|pair| helical_parameters(pair[0], pair[1], options))
        .collect()
}

fn validate_options(options: HelicalOptions) -> Result<(), HelicalError> {
    if options.frame_tolerance.is_finite() && options.frame_tolerance > 0.0 {
        Ok(())
    } else {
        Err(HelicalError::InvalidTolerance)
    }
}

fn validate_frame(frame: BaseFrame, tolerance: f64) -> Result<(), HelicalError> {
    let finite = [frame.origin, frame.x, frame.y, frame.z]
        .iter()
        .flatten()
        .all(|value| value.is_finite());
    let residuals = [
        norm_squared(frame.x) - 1.0,
        norm_squared(frame.y) - 1.0,
        norm_squared(frame.z) - 1.0,
        dot(frame.x, frame.y),
        dot(frame.x, frame.z),
        dot(frame.y, frame.z),
        dot(cross(frame.x, frame.y), frame.z) - 1.0,
    ];
    if finite && residuals.iter().all(|value| value.abs() <= tolerance) {
        Ok(())
    } else {
        Err(HelicalError::InvalidFrame)
    }
}

fn relative_rotation(first: BaseFrame, second: BaseFrame) -> [[f64; 3]; 3] {
    let left = [first.x, first.y, first.z];
    let right = [second.x, second.y, second.z];
    std::array::from_fn(|row| std::array::from_fn(|column| dot(left[row], right[column])))
}

fn rotation_vector(matrix: [[f64; 3]; 3], tolerance: f64) -> [f64; 3] {
    let antisymmetric = [
        matrix[2][1] - matrix[1][2],
        matrix[0][2] - matrix[2][0],
        matrix[1][0] - matrix[0][1],
    ];
    let sine = 0.5 * norm_squared(antisymmetric).sqrt();
    let cosine = ((matrix[0][0] + matrix[1][1] + matrix[2][2]) - 1.0) * 0.5;
    let angle = sine.atan2(cosine.clamp(-1.0, 1.0));
    if angle.abs() <= tolerance {
        return [0.0; 3];
    }
    let axis = if sine > tolerance {
        scale(antisymmetric, 0.5 / sine)
    } else {
        half_turn_axis(matrix)
    };
    scale(axis, angle)
}

fn half_turn_axis(matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let squared: [f64; 3] =
        std::array::from_fn(|axis| f64::midpoint(matrix[axis][axis], 1.0).max(0.0));
    let largest = if squared[1] > squared[0] {
        usize::from(squared[2] > squared[1]) + 1
    } else if squared[2] > squared[0] {
        2
    } else {
        0
    };
    let mut axis = [0.0; 3];
    axis[largest] = squared[largest].sqrt();
    if axis[largest] > 0.0 {
        for other in 0..3 {
            if other != largest {
                axis[other] =
                    (matrix[largest][other] + matrix[other][largest]) / (4.0 * axis[largest]);
            }
        }
    }
    axis
}

fn subtract(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|axis| left[axis] - right[axis])
}

fn dot(left: [f64; 3], right: [f64; 3]) -> f64 {
    (0..3).map(|axis| left[axis] * right[axis]).sum()
}

fn norm_squared(value: [f64; 3]) -> f64 {
    dot(value, value)
}

fn scale(value: [f64; 3], factor: f64) -> [f64; 3] {
    value.map(|component| component * factor)
}

fn cross(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]
}

#[cfg(test)]
#[path = "helical_tests.rs"]
mod tests;
