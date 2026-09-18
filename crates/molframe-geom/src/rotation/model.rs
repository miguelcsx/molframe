//! Validated rotations and intrinsic operations on `SO(3)`.
//!
//! Rotation matrices form a curved group rather than a vector space. Validation
//! prevents reflections and shears from entering APIs that promise rigid motion;
//! logarithms, exponentials and interpolation follow the shortest group geodesic.

use crate::numeric::exact_count;
use crate::{cross, dot, norm};
use std::f64::consts::PI;

/// Numerical controls for exponential, logarithmic and interpolation kernels.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotationOptions {
    /// Tolerance used to validate matrices constructed by the exponential map.
    pub matrix_tolerance: f64,
    /// Angle below which the stable small-angle formula is used.
    pub small_angle_tolerance: f64,
    /// Distance from π below which the stable half-turn logarithm is used.
    pub half_turn_tolerance: f64,
}

impl RotationOptions {
    /// Matrix tolerance in the named double-precision profile.
    pub const STANDARD_MATRIX_TOLERANCE: f64 = 1e-10;
    /// Small-angle boundary in the named double-precision profile.
    pub const STANDARD_SMALL_ANGLE_TOLERANCE: f64 = 1e-10;
    /// Half-turn boundary in the named double-precision profile.
    pub const STANDARD_HALF_TURN_TOLERANCE: f64 = 1e-8;

    /// Named double-precision profile used by convenience methods.
    #[must_use]
    pub const fn standard() -> Self {
        Self {
            matrix_tolerance: Self::STANDARD_MATRIX_TOLERANCE,
            small_angle_tolerance: Self::STANDARD_SMALL_ANGLE_TOLERANCE,
            half_turn_tolerance: Self::STANDARD_HALF_TURN_TOLERANCE,
        }
    }

    fn validate(self) -> Result<Self, RotationError> {
        if valid_tolerance(self.matrix_tolerance)
            && valid_tolerance(self.small_angle_tolerance)
            && valid_tolerance(self.half_turn_tolerance)
            && self.small_angle_tolerance * self.small_angle_tolerance <= self.matrix_tolerance
        {
            Ok(self)
        } else {
            Err(RotationError::InvalidOptions)
        }
    }
}

/// Explicit convergence and numerical controls for an intrinsic rotation mean.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotationMeanOptions {
    /// Geodesic update norm required for convergence.
    pub convergence_tolerance: f64,
    /// Maximum Karcher iterations.
    pub maximum_iterations: usize,
    /// Numerical profile for group exponential and logarithm operations.
    pub rotation: RotationOptions,
}

/// Why a matrix or intrinsic rotation computation was refused.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RotationError {
    /// A matrix or parameter contained NaN or infinity.
    NonFinite,
    /// The matrix was not orthonormal within the requested tolerance.
    NotOrthonormal,
    /// The matrix was a reflection rather than a rotation.
    Reflection,
    /// An iterative rotation mean did not converge.
    DidNotConverge,
    /// No rotations were supplied.
    Empty,
    /// A tolerance or iteration limit was invalid.
    InvalidOptions,
    /// The rotation count cannot be represented exactly by the numeric kernel.
    TooManyRotations,
}

/// A finite orthonormal three-dimensional rotation with determinant one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Rotation3 {
    matrix: [[f64; 3]; 3],
}

impl Rotation3 {
    /// Identity rotation.
    pub const IDENTITY: Self = Self {
        matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };

    /// Validates a rotation matrix.
    ///
    /// # Errors
    ///
    /// Refuses non-finite, non-orthonormal and reflecting matrices.
    pub fn from_matrix(matrix: [[f64; 3]; 3], tolerance: f64) -> Result<Self, RotationError> {
        if !valid_tolerance(tolerance) {
            return Err(RotationError::InvalidOptions);
        }
        if !matrix.iter().flatten().all(|v| v.is_finite()) {
            return Err(RotationError::NonFinite);
        }
        for row in &matrix {
            if (dot(*row, *row) - 1.0).abs() > tolerance {
                return Err(RotationError::NotOrthonormal);
            }
        }
        for left in 0..3 {
            for right in (left + 1)..3 {
                if dot(matrix[left], matrix[right]).abs() > tolerance {
                    return Err(RotationError::NotOrthonormal);
                }
            }
        }
        let determinant = determinant(matrix);
        if (determinant + 1.0).abs() <= tolerance {
            return Err(RotationError::Reflection);
        }
        if (determinant - 1.0).abs() > tolerance {
            return Err(RotationError::NotOrthonormal);
        }
        Ok(Self { matrix })
    }

    /// Constructs a rotation from an axis-angle tangent vector.
    ///
    /// The vector direction is the axis and its length is the angle in radians.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite vector.
    pub fn exp(tangent: [f64; 3]) -> Result<Self, RotationError> {
        Self::exp_unchecked(tangent, RotationOptions::standard())
    }

    /// Constructs a rotation with explicit numerical controls.
    ///
    /// # Errors
    ///
    /// Refuses non-finite vectors or invalid controls.
    pub fn exp_with_options(
        tangent: [f64; 3],
        options: RotationOptions,
    ) -> Result<Self, RotationError> {
        Self::exp_unchecked(tangent, options.validate()?)
    }

    fn exp_unchecked(tangent: [f64; 3], options: RotationOptions) -> Result<Self, RotationError> {
        if !tangent.iter().all(|value| value.is_finite()) {
            return Err(RotationError::NonFinite);
        }
        let angle = norm(tangent);
        if angle <= options.small_angle_tolerance {
            return Self::from_matrix(rodrigues_series(tangent), options.matrix_tolerance);
        }
        let axis = tangent.map(|value| value / angle);
        Self::from_matrix(rodrigues(axis, angle), options.matrix_tolerance)
    }

    /// The underlying row-major matrix.
    #[must_use]
    pub const fn matrix(self) -> [[f64; 3]; 3] {
        self.matrix
    }

    /// Applies the rotation to a double-precision vector.
    #[must_use]
    pub fn apply(self, vector: [f64; 3]) -> [f64; 3] {
        multiply_vector(self.matrix, vector)
    }

    /// The inverse rotation.
    #[must_use]
    pub fn inverse(self) -> Self {
        Self {
            matrix: transpose(self.matrix),
        }
    }

    /// Applies `self` and then `next`.
    #[must_use]
    pub fn then(self, next: Self) -> Self {
        Self {
            matrix: multiply(next.matrix, self.matrix),
        }
    }

    /// Axis-angle logarithm with magnitude in `[0, π]`.
    #[must_use]
    pub fn log(self) -> [f64; 3] {
        self.log_unchecked(RotationOptions::standard())
    }

    /// Axis-angle logarithm with explicit numerical controls.
    ///
    /// # Errors
    ///
    /// Refuses invalid controls.
    pub fn log_with_options(self, options: RotationOptions) -> Result<[f64; 3], RotationError> {
        Ok(self.log_unchecked(options.validate()?))
    }

    fn log_unchecked(self, options: RotationOptions) -> [f64; 3] {
        let cosine = ((trace(self.matrix) - 1.0) * 0.5).clamp(-1.0, 1.0);
        let angle = cosine.acos();
        if angle <= options.small_angle_tolerance {
            return [
                (self.matrix[2][1] - self.matrix[1][2]) * 0.5,
                (self.matrix[0][2] - self.matrix[2][0]) * 0.5,
                (self.matrix[1][0] - self.matrix[0][1]) * 0.5,
            ];
        }
        if (PI - angle).abs() <= options.half_turn_tolerance {
            return log_half_turn(self.matrix, angle, options.small_angle_tolerance);
        }
        let scale = angle / (2.0 * angle.sin());
        [
            scale * (self.matrix[2][1] - self.matrix[1][2]),
            scale * (self.matrix[0][2] - self.matrix[2][0]),
            scale * (self.matrix[1][0] - self.matrix[0][1]),
        ]
    }

    /// Bi-invariant geodesic distance in radians.
    #[must_use]
    pub fn distance(self, other: Self) -> f64 {
        norm(self.inverse().then(other).log())
    }

    /// Shortest constant-speed interpolation.
    ///
    /// # Errors
    ///
    /// Refuses non-finite interpolation parameters.
    pub fn interpolate(self, other: Self, amount: f64) -> Result<Self, RotationError> {
        self.interpolate_with_options(other, amount, RotationOptions::standard())
    }

    /// Shortest interpolation with explicit numerical controls.
    ///
    /// # Errors
    ///
    /// Refuses non-finite interpolation parameters or invalid controls.
    pub fn interpolate_with_options(
        self,
        other: Self,
        amount: f64,
        options: RotationOptions,
    ) -> Result<Self, RotationError> {
        let options = options.validate()?;
        if !amount.is_finite() {
            return Err(RotationError::NonFinite);
        }
        let tangent = self
            .inverse()
            .then(other)
            .log_unchecked(options)
            .map(|value| value * amount);
        let increment = Self::exp_unchecked(tangent, options)?;
        Ok(self.then(increment))
    }
}

/// Computes an intrinsic Karcher mean of rotations.
///
/// # Errors
///
/// Refuses empty input, invalid options or a mean that does not converge.
pub fn rotation_mean(
    rotations: &[Rotation3],
    tolerance: f64,
    maximum_iterations: usize,
) -> Result<Rotation3, RotationError> {
    rotation_mean_with_options(
        rotations,
        RotationMeanOptions {
            convergence_tolerance: tolerance,
            maximum_iterations,
            rotation: RotationOptions::standard(),
        },
    )
}

/// Computes an intrinsic Karcher mean with explicit numerical controls.
///
/// # Errors
///
/// Refuses empty input, invalid options or a mean that does not converge.
pub fn rotation_mean_with_options(
    rotations: &[Rotation3],
    options: RotationMeanOptions,
) -> Result<Rotation3, RotationError> {
    let Some(first) = rotations.first().copied() else {
        return Err(RotationError::Empty);
    };
    let rotation_options = options.rotation.validate()?;
    if !valid_tolerance(options.convergence_tolerance) || options.maximum_iterations == 0 {
        return Err(RotationError::InvalidOptions);
    }
    let count = exact_count(rotations.len()).ok_or(RotationError::TooManyRotations)?;
    let mut mean = first;
    for _ in 0..options.maximum_iterations {
        let mut average = [0.0; 3];
        for rotation in rotations {
            let tangent = mean
                .inverse()
                .then(*rotation)
                .log_unchecked(rotation_options);
            for axis in 0..3 {
                average[axis] += tangent[axis] / count;
            }
        }
        if norm(average) <= options.convergence_tolerance {
            return Ok(mean);
        }
        mean = mean.then(Rotation3::exp_unchecked(average, rotation_options)?);
    }
    Err(RotationError::DidNotConverge)
}

fn rodrigues(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
    let [x, y, z] = axis;
    let cosine = angle.cos();
    let sine = angle.sin();
    let one = 1.0 - cosine;
    [
        [
            cosine + x * x * one,
            x * y * one - z * sine,
            x * z * one + y * sine,
        ],
        [
            y * x * one + z * sine,
            cosine + y * y * one,
            y * z * one - x * sine,
        ],
        [
            z * x * one - y * sine,
            z * y * one + x * sine,
            cosine + z * z * one,
        ],
    ]
}

fn rodrigues_series(tangent: [f64; 3]) -> [[f64; 3]; 3] {
    let [x, y, z] = tangent;
    [[1.0, -z, y], [z, 1.0, -x], [-y, x, 1.0]]
}

fn log_half_turn(matrix: [[f64; 3]; 3], angle: f64, small_angle_tolerance: f64) -> [f64; 3] {
    let mut axis = [
        f64::midpoint(matrix[0][0], 1.0).max(0.0).sqrt(),
        f64::midpoint(matrix[1][1], 1.0).max(0.0).sqrt(),
        f64::midpoint(matrix[2][2], 1.0).max(0.0).sqrt(),
    ];
    if matrix[2][1] - matrix[1][2] < 0.0 {
        axis[0] = -axis[0];
    }
    if matrix[0][2] - matrix[2][0] < 0.0 {
        axis[1] = -axis[1];
    }
    if matrix[1][0] - matrix[0][1] < 0.0 {
        axis[2] = -axis[2];
    }
    if norm(axis) <= small_angle_tolerance {
        axis = [1.0, 0.0, 0.0];
    } else {
        let length = norm(axis);
        axis = axis.map(|value| value / length);
    }
    axis.map(|value| value * angle)
}

fn valid_tolerance(tolerance: f64) -> bool {
    tolerance.is_finite() && tolerance > 0.0 && tolerance < 1.0
}

fn multiply_vector(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| dot(row, vector))
}

fn multiply(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let columns = transpose(right);
    left.map(|row| columns.map(|column| dot(row, column)))
}

fn transpose(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [matrix[0][0], matrix[1][0], matrix[2][0]],
        [matrix[0][1], matrix[1][1], matrix[2][1]],
        [matrix[0][2], matrix[1][2], matrix[2][2]],
    ]
}

fn trace(matrix: [[f64; 3]; 3]) -> f64 {
    matrix[0][0] + matrix[1][1] + matrix[2][2]
}

fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    dot(matrix[0], cross(matrix[1], matrix[2]))
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
