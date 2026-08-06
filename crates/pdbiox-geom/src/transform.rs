//! Rigid transforms: a rotation and a translation.
//!
//! Rigid means distance-preserving, which is the property that makes these safe
//! to apply to a structure: no transform here can stretch a bond or invert a
//! centre. Composition and inversion are exact rather than approximate, so a
//! transform applied and then undone returns the original positions.

/// A rotation followed by a translation.
///
/// # Examples
///
/// ```
/// use pdbiox_geom::Rigid;
///
/// let shift = Rigid::translation([1.0, 2.0, 3.0]);
/// assert_eq!(shift.apply([0.0, 0.0, 0.0]), [1.0, 2.0, 3.0]);
/// assert_eq!(shift.inverse().apply([1.0, 2.0, 3.0]), [0.0, 0.0, 0.0]);
/// ```
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Rigid {
    /// The rotation, as rows.
    pub rotation: [[f64; 3]; 3],
    /// The translation applied after rotating.
    pub translation: [f64; 3],
}

impl Rigid {
    /// The transform that changes nothing.
    pub const IDENTITY: Self = Self {
        rotation: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        translation: [0.0; 3],
    };

    /// Creates a transform from its parts.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub const fn new(rotation: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    /// A pure translation.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub const fn translation(offset: [f64; 3]) -> Self {
        Self {
            rotation: Self::IDENTITY.rotation,
            translation: offset,
        }
    }

    /// Applies the transform to one position.
    ///
    /// Runs in `O(1)` time with fixed-size arithmetic and no allocation.
    #[must_use]
    #[inline]
    pub fn apply(&self, position: [f32; 3]) -> [f32; 3] {
        let point = [
            f64::from(position[0]),
            f64::from(position[1]),
            f64::from(position[2]),
        ];
        let transformed = rotate_and_translate(&self.rotation, point, self.translation);

        [
            transformed[0] as f32,
            transformed[1] as f32,
            transformed[2] as f32,
        ]
    }

    /// Applies the transform to every position in place.
    ///
    /// One pass, no allocation.
    ///
    /// Runs in `O(n)` time and `O(1)` auxiliary space.
    pub fn apply_all(&self, positions: &mut [[f32; 3]]) {
        for position in positions {
            *position = self.apply(*position);
        }
    }

    /// The transform that undoes this one.
    ///
    /// A rotation's inverse is its transpose, which is exact — so undoing a
    /// transform introduces no error beyond the rounding of the arithmetic.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    pub fn inverse(&self) -> Self {
        let rotation = transpose(&self.rotation);
        let rotated_translation = rotate(&rotation, self.translation);
        let translation = [
            -rotated_translation[0],
            -rotated_translation[1],
            -rotated_translation[2],
        ];

        Self {
            rotation,
            translation,
        }
    }

    /// The transform that applies `self` and then `next`.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    pub fn then(&self, next: &Self) -> Self {
        let rotation = multiply(&next.rotation, &self.rotation);
        let translation = rotate_and_translate(&next.rotation, self.translation, next.translation);

        Self {
            rotation,
            translation,
        }
    }

    /// The determinant of the rotation.
    ///
    /// One for a rotation, minus one for a reflection. Nothing in this crate
    /// produces a reflection, and checking is how that stays true.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn determinant(&self) -> f64 {
        let r = &self.rotation;

        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }
}

/// The transpose of a three-by-three, which is a rotation's inverse.
///
/// Runs in `O(1)` time and returns a stack-allocated matrix.
#[inline]
fn transpose(matrix: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [matrix[0][0], matrix[1][0], matrix[2][0]],
        [matrix[0][1], matrix[1][1], matrix[2][1]],
        [matrix[0][2], matrix[1][2], matrix[2][2]],
    ]
}

/// Multiplies two three-by-three matrices as `left * right`.
///
/// Runs in `O(1)` time with a fixed 27 scalar products and allocates no memory.
fn multiply(left: &[[f64; 3]; 3], right: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0; 3]; 3];

    for (row_index, out_row) in out.iter_mut().enumerate() {
        let left_row = &left[row_index];

        for (column, slot) in out_row.iter_mut().enumerate() {
            *slot = left_row[0] * at(&right[0], column)
                + left_row[1] * at(&right[1], column)
                + left_row[2] * at(&right[2], column);
        }
    }

    out
}

/// Applies a three-by-three rotation to a double-precision vector.
///
/// Runs in `O(1)` time and allocates no heap memory.
#[inline]
fn rotate(rotation: &[[f64; 3]; 3], vector: [f64; 3]) -> [f64; 3] {
    [
        rotation[0][0] * vector[0] + rotation[0][1] * vector[1] + rotation[0][2] * vector[2],
        rotation[1][0] * vector[0] + rotation[1][1] * vector[1] + rotation[1][2] * vector[2],
        rotation[2][0] * vector[0] + rotation[2][1] * vector[1] + rotation[2][2] * vector[2],
    ]
}

/// Applies a rotation and then adds a translation.
///
/// Runs in `O(1)` time and allocates no heap memory.
#[inline]
fn rotate_and_translate(
    rotation: &[[f64; 3]; 3],
    vector: [f64; 3],
    translation: [f64; 3],
) -> [f64; 3] {
    let rotated = rotate(rotation, vector);

    [
        rotated[0] + translation[0],
        rotated[1] + translation[1],
        rotated[2] + translation[2],
    ]
}

/// One entry of a row, reading past the end as zero.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn at(row: &[f64; 3], column: usize) -> f64 {
    match row.get(column) {
        Some(value) => *value,
        None => 0.0,
    }
}

impl Default for Rigid {
    /// Returns [`Rigid::IDENTITY`].
    ///
    /// Runs in `O(1)` time and allocates no memory.
    fn default() -> Self {
        Self::IDENTITY
    }
}

#[cfg(test)]
#[path = "transform_tests.rs"]
mod tests;
