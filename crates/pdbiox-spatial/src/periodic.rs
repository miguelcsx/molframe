//! Orthorhombic and triclinic minimum-image displacement.

use crate::SpatialError;
use pdbiox_core::structure::UnitCell;

/// An invertible unit-cell basis and its inverse.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct PeriodicBox {
    basis: [[f64; 3]; 3],
    inverse: [[f64; 3]; 3],
}

impl PeriodicBox {
    /// Builds a periodic box from crystallographic lengths and angles.
    ///
    /// # Errors
    ///
    /// Returns [`SpatialError::InvalidCell`] when any parameter is non-finite,
    /// an edge is not positive, or the vectors are coplanar.
    pub fn from_cell(cell: UnitCell) -> Result<Self, SpatialError> {
        if cell
            .lengths
            .iter()
            .chain(&cell.angles)
            .any(|value| !value.is_finite())
            || cell.lengths.iter().any(|length| *length <= 0.0)
        {
            return Err(SpatialError::InvalidCell);
        }
        let [a, b, c] = cell.lengths;
        let [alpha, beta, gamma] = cell.angles.map(f64::to_radians);
        let sin_gamma = gamma.sin();
        if sin_gamma.abs() <= f64::EPSILON {
            return Err(SpatialError::InvalidCell);
        }
        let cx = c * beta.cos();
        let cy = c * (alpha.cos() - beta.cos() * gamma.cos()) / sin_gamma;
        let cz_squared = c.mul_add(c, -cx.mul_add(cx, cy * cy));
        if cz_squared <= f64::EPSILON {
            return Err(SpatialError::InvalidCell);
        }
        let basis = [
            [a, b * gamma.cos(), cx],
            [0.0, b * sin_gamma, cy],
            [0.0, 0.0, cz_squared.sqrt()],
        ];
        let Some(inverse) = inverse3(basis) else {
            return Err(SpatialError::InvalidCell);
        };
        Ok(Self { basis, inverse })
    }

    /// The shortest periodic displacement from `left` to `right`.
    ///
    /// The nearest rounded fractional image is sufficient for an orthogonal
    /// cell. In a skewed cell a neighbouring lattice translation can be
    /// shorter, so all twenty-seven candidates around it are checked in a
    /// fixed order.
    #[must_use]
    pub fn displacement(&self, left: [f32; 3], right: [f32; 3]) -> [f32; 3] {
        let cartesian = [
            f64::from(right[0] - left[0]),
            f64::from(right[1] - left[1]),
            f64::from(right[2] - left[2]),
        ];
        let fractional = multiply(self.inverse, cartesian);
        let centre = fractional.map(f64::round);
        let mut best = cartesian;
        let mut best_squared = squared(best);
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let image = [
                        fractional[0] - centre[0] - f64::from(i),
                        fractional[1] - centre[1] - f64::from(j),
                        fractional[2] - centre[2] - f64::from(k),
                    ];
                    let candidate = multiply(self.basis, image);
                    let candidate_squared = squared(candidate);
                    if candidate_squared < best_squared {
                        best = candidate;
                        best_squared = candidate_squared;
                    }
                }
            }
        }
        best.map(|component| component as f32)
    }

    /// Squared minimum-image separation.
    #[must_use]
    pub fn distance_squared(&self, left: [f32; 3], right: [f32; 3]) -> f32 {
        self.displacement(left, right)
            .iter()
            .map(|component| component * component)
            .sum()
    }
}

fn multiply(matrix: [[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    matrix.map(|row| row.iter().zip(vector).map(|(a, b)| a * b).sum())
}

fn squared(vector: [f64; 3]) -> f64 {
    vector.iter().map(|component| component * component).sum()
}

fn inverse3(matrix: [[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let determinant = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1])
        - matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0])
        + matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
    if !determinant.is_finite() || determinant.abs() <= f64::EPSILON {
        return None;
    }
    let reciprocal = determinant.recip();
    Some([
        [
            (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * reciprocal,
            (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * reciprocal,
            (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * reciprocal,
        ],
        [
            (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * reciprocal,
            (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * reciprocal,
            (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * reciprocal,
        ],
        [
            (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * reciprocal,
            (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * reciprocal,
            (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * reciprocal,
        ],
    ])
}

#[cfg(test)]
#[path = "periodic_tests.rs"]
mod tests;
