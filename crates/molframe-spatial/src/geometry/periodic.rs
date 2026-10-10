//! Orthorhombic and triclinic minimum-image displacement.

use super::reduction::{self, Matrix};
use crate::SpatialError;
use crate::numeric::{f64_f32, i64_f64, rounded_i64};
use molframe_core::structure::UnitCell;

const ORTHOGONAL_ANGLE_TOLERANCE: f64 = 1.0e-10;

/// Safety cap on the per-axis offset search; reduced bases need at most 2.
const MAX_SEARCH_REACH: f64 = 16.0;

/// An invertible unit-cell basis and its inverse.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PeriodicBox {
    basis: [[f64; 3]; 3],
    inverse: [[f64; 3]; 3],
    orthogonal: bool,
    reduction: Reduction,
}

/// Integer change of basis to a lattice-reduced basis, `reduced = basis * T`.
///
/// Only the integer matrices are stored; the reduced basis is rebuilt per
/// search. That keeps the box small enough to pass by value.
#[derive(Clone, Copy, PartialEq, Debug)]
struct Reduction {
    transform: [[i32; 3]; 3],
    inverse_transform: [[i32; 3]; 3],
}

/// The reduced basis, its inverse and the row norms bounding the search.
struct ReducedFrame {
    basis: Matrix,
    inverse: Matrix,
    row_norms: [f64; 3],
}

impl Reduction {
    const IDENTITY: Self = Self {
        transform: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        inverse_transform: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    };

    fn new(basis: Matrix) -> Self {
        let reduced = reduction::reduce(basis);
        reduction::narrow(reduced.transform).map_or(
            Self::IDENTITY,
            |(transform, inverse_transform)| Self {
                transform,
                inverse_transform,
            },
        )
    }

    fn frame(&self, basis: Matrix, inverse: Matrix) -> ReducedFrame {
        let to_f64 = |matrix: [[i32; 3]; 3]| matrix.map(|row| row.map(f64::from));
        let reduced = multiply_matrices(basis, to_f64(self.transform));
        let reduced_inverse = multiply_matrices(to_f64(self.inverse_transform), inverse);
        ReducedFrame {
            basis: reduced,
            inverse: reduced_inverse,
            row_norms: reduced_inverse.map(|row| squared(row).sqrt()),
        }
    }

    /// Maps an integer shift in reduced coordinates to the original basis.
    fn original_shift(&self, shift: [i64; 3]) -> [i64; 3] {
        std::array::from_fn(|row| {
            (0..3).fold(0_i64, |total, column| {
                total.saturating_add(
                    i64::from(self.transform[row][column]).saturating_mul(shift[column]),
                )
            })
        })
    }
}

fn small_integer(value: f64) -> i32 {
    match i32::try_from(rounded_i64(value)) {
        Ok(converted) => converted,
        Err(_) => 0,
    }
}

fn multiply_matrices(left: Matrix, right: Matrix) -> Matrix {
    std::array::from_fn(|row| {
        std::array::from_fn(|column| (0..3).map(|k| left[row][k] * right[k][column]).sum())
    })
}

/// A shortest displacement on a periodic flat torus and the chosen image.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PeriodicImage {
    /// Cartesian displacement from the left point to the chosen right image.
    pub displacement: [f32; 3],
    /// Lattice translation subtracted from the un-wrapped right point.
    pub lattice_shift: [i64; 3],
}

impl PeriodicBox {
    /// Builds a periodic box from crystallographic lengths and angles.
    ///
    /// # Errors
    ///
    /// Returns [`SpatialError::InvalidCell`] when any parameter is non-finite,
    /// an edge is not positive, or the vectors are coplanar.
    pub fn from_cell(cell: UnitCell) -> Result<Self, SpatialError> {
        validate_cell(&cell)?;

        let orthogonal = cell
            .angles
            .iter()
            .all(|angle| (*angle - 90.0).abs() < ORTHOGONAL_ANGLE_TOLERANCE);

        let basis = basis_from_cell(&cell)?;
        let Some(inverse) = inverse3(basis) else {
            return Err(SpatialError::InvalidCell);
        };

        let reduction = Reduction::new(basis);

        Ok(Self {
            basis,
            inverse,
            orthogonal,
            reduction,
        })
    }

    /// The shortest periodic displacement from `left` to `right`.
    ///
    /// The nearest rounded fractional image is sufficient for an orthogonal
    /// cell. In a skewed cell a neighbouring lattice translation can be
    /// shorter, so a lattice-reduced basis is searched over a provably
    /// sufficient neighbourhood and the result is exact.
    #[must_use]
    pub fn displacement(&self, left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
        self.minimum_image(left, right).displacement
    }

    /// The shortest displacement and exact lattice image selected for it.
    ///
    /// Exact distance ties retain the deposited image. Triclinic candidates
    /// are visited in a fixed order, so the selected image is deterministic.
    #[must_use]
    pub fn minimum_image(&self, left: [f32; 3], right: [f32; 3]) -> PeriodicImage {
        let (best, lattice_shift) = self.image_f64(left.map(f64::from), right.map(f64::from));

        PeriodicImage {
            displacement: best.map(f64_f32),
            lattice_shift,
        }
    }

    /// Minimum-image displacement without rounding generated sample positions
    /// to stored-coordinate precision. Image selection matches `minimum_image`.
    #[must_use]
    pub fn displacement_f64(&self, left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
        self.image_f64(left, right).0
    }

    fn image_f64(&self, left: [f64; 3], right: [f64; 3]) -> ([f64; 3], [i64; 3]) {
        let cartesian = cartesian_delta(left, right);
        if self.orthogonal {
            let fractional = multiply(self.inverse, cartesian);
            orthogonal_image(self.basis, cartesian, fractional)
        } else {
            triclinic_image(&self.reduction, self.basis, self.inverse, cartesian)
        }
    }

    /// Squared minimum-image separation.
    #[must_use]
    pub fn distance_squared(&self, left: [f32; 3], right: [f32; 3]) -> f32 {
        let displacement = self.displacement(left, right);

        displacement[0] * displacement[0]
            + displacement[1] * displacement[1]
            + displacement[2] * displacement[2]
    }

    /// Converts a Cartesian position to fractional cell coordinates.
    #[must_use]
    pub fn fractional(&self, position: [f32; 3]) -> [f64; 3] {
        multiply(self.inverse, position.map(f64::from))
    }

    /// Converts fractional cell coordinates to Cartesian ångström.
    #[must_use]
    pub fn cartesian(&self, fractional: [f64; 3]) -> [f32; 3] {
        multiply(self.basis, fractional).map(f64_f32)
    }

    /// Wraps one position into the primary fractional interval `[0, 1)`.
    #[must_use]
    pub fn wrap(&self, position: [f32; 3]) -> [f32; 3] {
        let fractional = self.fractional(position).map(|value| value.rem_euclid(1.0));

        self.cartesian(fractional)
    }

    /// Interpolates along the shortest periodic geodesic and wraps the result.
    ///
    /// Returns `None` for a non-finite interpolation parameter. Values outside
    /// `[0, 1]` extrapolate along the same selected image.
    #[must_use]
    pub fn interpolate(&self, left: [f32; 3], right: [f32; 3], amount: f64) -> Option<[f32; 3]> {
        if !amount.is_finite() {
            return None;
        }
        let displacement = self.minimum_image(left, right).displacement;
        let point = std::array::from_fn(|axis| {
            f64_f32(f64::from(left[axis]) + amount * f64::from(displacement[axis]))
        });
        Some(self.wrap(point))
    }

    /// Maximum fractional-coordinate changes produced by a Cartesian radius.
    ///
    /// The bounds are the Euclidean norms of the reciprocal-basis rows times
    /// the radius. They conservatively bound cell and tree image searches.
    #[must_use]
    pub(crate) fn fractional_cutoff_bounds(&self, cutoff: f32) -> [f64; 3] {
        self.inverse
            .map(|row| squared(row).sqrt() * f64::from(cutoff))
    }

    /// Translates a Cartesian point by an exact lattice shift.
    pub(crate) fn translated(
        &self,
        point: [f32; 3],
        shift: [i64; 3],
    ) -> Result<[f32; 3], SpatialError> {
        let convert = |value| i64_f64(value).ok_or(SpatialError::NumericRangeExceeded);
        let translation = multiply(
            self.basis,
            [convert(shift[0])?, convert(shift[1])?, convert(shift[2])?],
        );
        Ok(std::array::from_fn(|axis| {
            f64_f32(f64::from(point[axis]) + translation[axis])
        }))
    }
}

/// Validates finite cell lengths/angles and positive edge lengths.
///
/// Runtime and auxiliary space are `O(1)`.
fn validate_cell(cell: &UnitCell) -> Result<(), SpatialError> {
    let finite = cell
        .lengths
        .iter()
        .chain(&cell.angles)
        .all(|value| value.is_finite());

    let positive = cell.lengths.iter().all(|length| *length > 0.0);

    if finite && positive {
        Ok(())
    } else {
        Err(SpatialError::InvalidCell)
    }
}

/// Builds the crystallographic Cartesian basis matrix.
///
/// # Errors
///
/// Returns [`SpatialError::InvalidCell`] for singular or numerically invalid
/// angle combinations.
fn basis_from_cell(cell: &UnitCell) -> Result<[[f64; 3]; 3], SpatialError> {
    let [a, b, c] = cell.lengths;
    let [alpha, beta, gamma] = cell.angles.map(f64::to_radians);

    let sin_gamma = gamma.sin();

    if !sin_gamma.is_finite() || sin_gamma.abs() <= f64::EPSILON {
        return Err(SpatialError::InvalidCell);
    }

    let cos_gamma = gamma.cos();
    let cos_beta = beta.cos();
    let cos_alpha = alpha.cos();

    let cx = c * cos_beta;
    let cy = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
    let cz_squared = c.mul_add(c, -cx.mul_add(cx, cy * cy));

    if !cz_squared.is_finite() || cz_squared <= f64::EPSILON {
        return Err(SpatialError::InvalidCell);
    }

    Ok([
        [a, b * cos_gamma, cx],
        [0.0, b * sin_gamma, cy],
        [0.0, 0.0, cz_squared.sqrt()],
    ])
}

/// Computes Cartesian `right - left` in `f64` arithmetic.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn cartesian_delta(left: [f64; 3], right: [f64; 3]) -> [f64; 3] {
    [right[0] - left[0], right[1] - left[1], right[2] - left[2]]
}

/// Computes the minimum image for an orthogonal cell.
///
/// The original Cartesian displacement is retained on exact distance ties,
/// preserving the displacement sign convention at half-cell boundaries.
fn orthogonal_image(
    basis: [[f64; 3]; 3],
    cartesian: [f64; 3],
    fractional: [f64; 3],
) -> ([f64; 3], [i64; 3]) {
    let rounded = fractional.map(f64::round);
    let wrapped = [
        fractional[0] - rounded[0],
        fractional[1] - rounded[1],
        fractional[2] - rounded[2],
    ];

    let candidate = multiply(basis, wrapped);

    if squared(candidate) < squared(cartesian) {
        (candidate, integer_shift(rounded))
    } else {
        (cartesian, [0; 3])
    }
}

/// Computes the minimum image for a skewed triclinic cell.
///
/// The search runs in a lattice-reduced basis. The rounded fractional position
/// there is a candidate at distance `r`; the true nearest image lies within
/// `r` of the target, hence within `2r` of that candidate, which bounds each
/// reduced-coordinate offset by `2r` times the norm of the matching row of the
/// reduced inverse. Every offset inside those bounds is evaluated, so the
/// result is exact for any valid cell. Runtime is bounded by the reduced
/// basis, not by the skew of the cell parameters.
fn triclinic_image(
    reduction: &Reduction,
    basis: Matrix,
    inverse: Matrix,
    cartesian: [f64; 3],
) -> ([f64; 3], [i64; 3]) {
    let frame = reduction.frame(basis, inverse);
    let fractional = multiply(frame.inverse, cartesian);
    let centre = fractional.map(f64::round);
    let base = [
        fractional[0] - centre[0],
        fractional[1] - centre[1],
        fractional[2] - centre[2],
    ];

    let rounded = multiply(frame.basis, base);
    let radius = squared(rounded).sqrt();
    if !radius.is_finite() {
        return (cartesian, [0; 3]);
    }

    let reach = frame.row_norms.map(|norm| {
        let bound = (2.0 * radius * norm).mul_add(1.0 + 1.0e-9, 1.0e-9).ceil();
        small_integer(bound.clamp(0.0, MAX_SEARCH_REACH))
    });

    let mut best = cartesian;
    let mut best_squared = squared(best);
    let mut best_shift = [0; 3];

    for first in -reach[0]..=reach[0] {
        for second in -reach[1]..=reach[1] {
            for third in -reach[2]..=reach[2] {
                let offset = [first, second, third].map(f64::from);
                let image = [
                    base[0] - offset[0],
                    base[1] - offset[1],
                    base[2] - offset[2],
                ];
                let candidate = multiply(frame.basis, image);
                let candidate_squared = squared(candidate);

                if candidate_squared < best_squared {
                    best = candidate;
                    best_squared = candidate_squared;
                    best_shift = reduction.original_shift([
                        rounded_i64(centre[0]).saturating_add(i64::from(first)),
                        rounded_i64(centre[1]).saturating_add(i64::from(second)),
                        rounded_i64(centre[2]).saturating_add(i64::from(third)),
                    ]);
                }
            }
        }
    }

    (best, best_shift)
}

fn integer_shift(shift: [f64; 3]) -> [i64; 3] {
    shift.map(rounded_i64)
}

/// Multiplies a three-by-three matrix by a three-component vector.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn multiply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    [
        matrix[0][0] * vector[0] + matrix[0][1] * vector[1] + matrix[0][2] * vector[2],
        matrix[1][0] * vector[0] + matrix[1][1] * vector[1] + matrix[1][2] * vector[2],
        matrix[2][0] * vector[0] + matrix[2][1] * vector[1] + matrix[2][2] * vector[2],
    ]
}

/// Returns a three-dimensional squared vector length.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn squared(vector: [f64; 3]) -> f64 {
    vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]
}

/// Computes all cofactors of a three-by-three matrix.
///
/// Runtime and auxiliary space are `O(1)`.
fn cofactors3(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1],
            matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2],
            matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0],
        ],
        [
            matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2],
            matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0],
            matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1],
        ],
        [
            matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1],
            matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2],
            matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0],
        ],
    ]
}

/// Computes a three-by-three inverse using the adjugate formula.
///
/// Returns `None` for non-finite or effectively singular matrices. Runtime and
/// auxiliary space are `O(1)`.
fn inverse3(matrix: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let cofactors = cofactors3(matrix);

    let determinant = matrix[0][0] * cofactors[0][0]
        + matrix[0][1] * cofactors[0][1]
        + matrix[0][2] * cofactors[0][2];

    if !determinant.is_finite() || determinant.abs() <= f64::EPSILON {
        return None;
    }

    let reciprocal = determinant.recip();

    Some([
        [
            cofactors[0][0] * reciprocal,
            cofactors[1][0] * reciprocal,
            cofactors[2][0] * reciprocal,
        ],
        [
            cofactors[0][1] * reciprocal,
            cofactors[1][1] * reciprocal,
            cofactors[2][1] * reciprocal,
        ],
        [
            cofactors[0][2] * reciprocal,
            cofactors[1][2] * reciprocal,
            cofactors[2][2] * reciprocal,
        ],
    ])
}

#[cfg(test)]
#[path = "periodic_tests.rs"]
mod tests;
