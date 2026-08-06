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
    #[must_use]
    pub const fn new(rotation: [[f64; 3]; 3], translation: [f64; 3]) -> Self {
        Self {
            rotation,
            translation,
        }
    }

    /// A pure translation.
    #[must_use]
    pub const fn translation(offset: [f64; 3]) -> Self {
        Self {
            rotation: Self::IDENTITY.rotation,
            translation: offset,
        }
    }

    /// Applies the transform to one position.
    #[must_use]
    pub fn apply(&self, position: [f32; 3]) -> [f32; 3] {
        let mut out = [0.0f32; 3];
        for (slot, (row, shift)) in out
            .iter_mut()
            .zip(self.rotation.iter().zip(self.translation))
        {
            let mut value = shift;
            for (coefficient, component) in row.iter().zip(position) {
                value += coefficient * f64::from(component);
            }
            *slot = value as f32;
        }
        out
    }

    /// Applies the transform to every position in place.
    ///
    /// One pass, no allocation.
    pub fn apply_all(&self, positions: &mut [[f32; 3]]) {
        for position in positions {
            *position = self.apply(*position);
        }
    }

    /// The transform that undoes this one.
    ///
    /// A rotation's inverse is its transpose, which is exact — so undoing a
    /// transform introduces no error beyond the rounding of the arithmetic.
    #[must_use]
    pub fn inverse(&self) -> Self {
        let rotation = transpose(&self.rotation);
        let mut translation = [0.0f64; 3];
        for (slot, row) in translation.iter_mut().zip(&rotation) {
            let mut value = 0.0;
            for (coefficient, component) in row.iter().zip(self.translation) {
                value += coefficient * component;
            }
            *slot = -value;
        }
        Self {
            rotation,
            translation,
        }
    }

    /// The transform that applies `self` and then `next`.
    #[must_use]
    pub fn then(&self, next: &Self) -> Self {
        let mut rotation = [[0.0f64; 3]; 3];
        for (out_row, next_row) in rotation.iter_mut().zip(&next.rotation) {
            for (column, slot) in out_row.iter_mut().enumerate() {
                let mut value = 0.0;
                for (coefficient, own_row) in next_row.iter().zip(&self.rotation) {
                    value += coefficient * at(own_row, column);
                }
                *slot = value;
            }
        }
        let mut translation = next.translation;
        for (slot, row) in translation.iter_mut().zip(&next.rotation) {
            for (coefficient, component) in row.iter().zip(self.translation) {
                *slot += coefficient * component;
            }
        }
        Self {
            rotation,
            translation,
        }
    }

    /// The determinant of the rotation.
    ///
    /// One for a rotation, minus one for a reflection. Nothing in this crate
    /// produces a reflection, and checking is how that stays true.
    #[must_use]
    pub fn determinant(&self) -> f64 {
        let r = &self.rotation;
        r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
    }
}

/// The transpose of a three-by-three, which is a rotation's inverse.
fn transpose(matrix: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0f64; 3]; 3];
    for (row, source) in matrix.iter().enumerate() {
        for (column, value) in source.iter().enumerate() {
            if let Some(slot) = out.get_mut(column).and_then(|out_row| out_row.get_mut(row)) {
                *slot = *value;
            }
        }
    }
    out
}

/// One entry of a row, reading past the end as zero.
fn at(row: &[f64; 3], column: usize) -> f64 {
    match row.get(column) {
        Some(value) => *value,
        None => 0.0,
    }
}

impl Default for Rigid {
    fn default() -> Self {
        Self::IDENTITY
    }
}

#[cfg(test)]
#[path = "transform_tests.rs"]
mod tests;
