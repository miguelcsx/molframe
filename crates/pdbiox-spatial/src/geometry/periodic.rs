//! Orthorhombic and triclinic minimum-image displacement.

use crate::SpatialError;
use crate::numeric::{f64_f32, i64_f64, rounded_i64};
use pdbiox_core::structure::UnitCell;

const ORTHOGONAL_ANGLE_TOLERANCE: f64 = 1.0e-10;

const IMAGE_OFFSETS: [[f64; 3]; 27] = [
    [-1.0, -1.0, -1.0],
    [-1.0, -1.0, 0.0],
    [-1.0, -1.0, 1.0],
    [-1.0, 0.0, -1.0],
    [-1.0, 0.0, 0.0],
    [-1.0, 0.0, 1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, 1.0, 0.0],
    [-1.0, 1.0, 1.0],
    [0.0, -1.0, -1.0],
    [0.0, -1.0, 0.0],
    [0.0, -1.0, 1.0],
    [0.0, 0.0, -1.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 1.0, -1.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 1.0],
    [1.0, -1.0, -1.0],
    [1.0, -1.0, 0.0],
    [1.0, -1.0, 1.0],
    [1.0, 0.0, -1.0],
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, -1.0],
    [1.0, 1.0, 0.0],
    [1.0, 1.0, 1.0],
];

/// An invertible unit-cell basis and its inverse.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PeriodicBox {
    basis: [[f64; 3]; 3],
    inverse: [[f64; 3]; 3],
    orthogonal: bool,
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

        Ok(Self {
            basis,
            inverse,
            orthogonal,
        })
    }

    /// The shortest periodic displacement from `left` to `right`.
    ///
    /// The nearest rounded fractional image is sufficient for an orthogonal
    /// cell. In a skewed cell a neighbouring lattice translation can be
    /// shorter, so all twenty-seven candidates around it are checked in a
    /// fixed order.
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
        let fractional = multiply(self.inverse, cartesian);
        if self.orthogonal {
            orthogonal_image(self.basis, cartesian, fractional)
        } else {
            triclinic_image(self.basis, cartesian, fractional)
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
/// All 27 translations around the rounded fractional image are evaluated in a
/// deterministic order. Runtime and auxiliary space are `O(1)`.
fn triclinic_image(
    basis: [[f64; 3]; 3],
    cartesian: [f64; 3],
    fractional: [f64; 3],
) -> ([f64; 3], [i64; 3]) {
    let centre = fractional.map(f64::round);
    let base = [
        fractional[0] - centre[0],
        fractional[1] - centre[1],
        fractional[2] - centre[2],
    ];

    let mut best = cartesian;
    let mut best_squared = squared(best);
    let mut best_shift = [0; 3];

    for offset in IMAGE_OFFSETS {
        let image = [
            base[0] - offset[0],
            base[1] - offset[1],
            base[2] - offset[2],
        ];

        let candidate = multiply(basis, image);
        let candidate_squared = squared(candidate);

        if candidate_squared < best_squared {
            best = candidate;
            best_squared = candidate_squared;
            best_shift = integer_shift([
                centre[0] + offset[0],
                centre[1] + offset[1],
                centre[2] + offset[2],
            ]);
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
