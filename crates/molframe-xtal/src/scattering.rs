//! X-ray scattering factors from the Cromer–Mann Gaussian expansion.

use crate::scattering_table::IT92;
use molframe_core::Element;

/// Neutral-atom X-ray form factor as four Gaussians and a constant.
///
/// `f(s) = c + Σ aᵢ exp(−bᵢ s²)` with `s = sin θ / λ` in reciprocal ångström.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GaussianFormFactor {
    a: [f64; 4],
    b: [f64; 4],
    c: f64,
}

impl GaussianFormFactor {
    /// The tabulated neutral-atom form factor of `element`, absent beyond
    /// californium and for the placeholder element.
    #[must_use]
    pub fn xray(element: Element) -> Option<Self> {
        let row = IT92.get(usize::from(element.atomic_number()))?;
        if row.iter().all(|value| *value == 0.0) {
            return None;
        }
        Some(Self {
            a: [row[0], row[1], row[2], row[3]],
            b: [row[4], row[5], row[6], row[7]],
            c: row[8],
        })
    }

    /// Scattering amplitude at `(sin θ / λ)²`, in electrons.
    #[must_use]
    pub fn value(&self, stol2: f64) -> f64 {
        self.a
            .iter()
            .zip(&self.b)
            .fold(self.c, |sum, (a, b)| sum + a * (-b * stol2).exp())
    }
}

#[cfg(test)]
#[path = "scattering_tests.rs"]
mod tests;
