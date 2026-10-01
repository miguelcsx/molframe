//! Reciprocal-space classification using the existing exact symmetry operators.

use crate::{SymmetryOperation, SymmetrySet};
use molframe_core::Diagnostic;

/// Space-group constraints for one reciprocal-lattice reflection.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReflectionSymmetry {
    /// Some rotational operation maps the reflection to its Friedel mate.
    pub centric: bool,
    /// A stabilizing operation introduces a non-integral translation phase.
    pub systematically_absent: bool,
    /// Number of complete group operations fixing the indices, including centring.
    pub epsilon_factor: usize,
}

impl SymmetryOperation {
    /// Applies the transposed direct-space rotation to Miller indices.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic if the transformed indices cannot fit `i32`.
    pub fn apply_to_hkl(&self, hkl: [i32; 3]) -> Result<[i32; 3], Diagnostic> {
        let mut result = [0; 3];
        for (axis, output) in result.iter_mut().enumerate() {
            let value: i128 = (0..3)
                .map(|row| i128::from(self.rotation[row][axis]) * i128::from(hkl[row]))
                .sum();
            *output = i32::try_from(value)
                .map_err(|_| crate::symmetry::symmetry_error("Miller index overflow"))?;
        }
        Ok(result)
    }

    /// Whether the translation changes this reflection's phase modulo a full turn.
    ///
    /// Integer arithmetic preserves screw/glide extinction conditions even for
    /// large negative indices and independently reduced translation fractions.
    #[must_use]
    pub fn has_reflection_phase_shift(&self, hkl: [i32; 3]) -> bool {
        let mut numerator = 0_i128;
        let mut denominator = 1_i128;
        for (index, translation) in hkl.into_iter().zip(self.translation) {
            let next_denominator = i128::from(translation.denominator());
            let next_numerator = (i128::from(index) * i128::from(translation.numerator()))
                .rem_euclid(next_denominator);
            let divisor = gcd(denominator, next_denominator);
            let scale = next_denominator / divisor;
            numerator = numerator * scale + next_numerator * (denominator / divisor);
            denominator *= scale;
            numerator = numerator.rem_euclid(denominator);
        }
        numerator != 0
    }
}

impl SymmetrySet {
    /// Classifies a reflection using every explicit operation, including centring.
    ///
    /// The origin is centric, never absent and has epsilon equal to group order.
    /// The operations must represent the complete space group, not only its
    /// primitive representatives.
    ///
    /// # Errors
    ///
    /// Returns a diagnostic for absent operations or unrepresentable indices.
    pub fn reflection_symmetry(&self, hkl: [i32; 3]) -> Result<ReflectionSymmetry, Diagnostic> {
        if self.operations().is_empty() {
            return Err(crate::symmetry::symmetry_error(
                "missing symmetry operations",
            ));
        }
        let mut result = ReflectionSymmetry {
            centric: false,
            systematically_absent: false,
            epsilon_factor: 0,
        };
        for operation in self.operations() {
            let transformed = operation.apply_to_hkl(hkl)?;
            result.centric |= transformed
                .into_iter()
                .zip(hkl)
                .all(|(a, b)| i64::from(a) == -i64::from(b));
            if transformed == hkl {
                result.epsilon_factor += 1;
                result.systematically_absent |= operation.has_reflection_phase_shift(hkl);
            }
        }
        Ok(result)
    }
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a
}

#[cfg(test)]
#[path = "reflection_symmetry_tests.rs"]
mod tests;
