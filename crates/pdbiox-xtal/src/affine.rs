//! General Cartesian affine transforms.

use crate::numeric::f64_to_f32;

/// A three-by-three linear map followed by a translation.
///
/// Unlike a rigid transform, an NCS operator is not required by the `PDBx`
/// dictionary to be a proper rotation. This type therefore preserves the
/// recorded matrix exactly instead of projecting it onto a rotation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AffineTransform {
    /// Linear matrix, as rows.
    pub matrix: [[f64; 3]; 3],
    /// Translation applied after the matrix.
    pub translation: [f64; 3],
}

impl AffineTransform {
    /// Transform that changes nothing.
    pub const IDENTITY: Self = Self {
        matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translation: [0.0; 3],
    };

    /// Creates a transform from its recorded parts.
    #[must_use]
    pub const fn new(matrix: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            matrix,
            translation,
        }
    }

    /// Whether every recorded coefficient is finite.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        self.matrix.iter().flatten().all(|value| value.is_finite())
            && self.translation.iter().all(|value| value.is_finite())
    }

    /// Applies the transform to one Cartesian position.
    #[must_use]
    pub fn apply(&self, position: [f32; 3]) -> [f32; 3] {
        let point = position.map(f64::from);
        let mut output = [0.0_f32; 3];
        for (row, value) in output.iter_mut().enumerate() {
            *value = f64_to_f32(
                self.matrix[row][0] * point[0]
                    + self.matrix[row][1] * point[1]
                    + self.matrix[row][2] * point[2]
                    + self.translation[row],
            );
        }
        output
    }

    /// Transform that applies `self` and then `next`.
    #[must_use]
    pub fn then(&self, next: &Self) -> Self {
        let mut matrix = [[0.0; 3]; 3];
        for (row, values) in matrix.iter_mut().enumerate() {
            for (column, value) in values.iter_mut().enumerate() {
                *value = (0..3)
                    .map(|axis| next.matrix[row][axis] * self.matrix[axis][column])
                    .sum();
            }
        }
        let translated = next.apply_f64(self.translation);
        Self::new(matrix, translated)
    }

    fn apply_f64(&self, point: [f64; 3]) -> [f64; 3] {
        let mut output = [0.0; 3];
        for (row, value) in output.iter_mut().enumerate() {
            *value = self.matrix[row][0] * point[0]
                + self.matrix[row][1] * point[1]
                + self.matrix[row][2] * point[2]
                + self.translation[row];
        }
        output
    }
}
