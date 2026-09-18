//! Sugar-ring pucker from the five endocyclic torsions.
//!
//! A furanose ring does not lie flat; it puckers, and the pucker is captured by
//! two numbers derived from the ring's five endocyclic torsions ν0–ν4. The phase
//! angle P says which atoms rise out of the mean plane — near 18° is the C3'-endo
//! form common in A-form nucleic acids, near 162° the C2'-endo of B-form — and
//! the amplitude says how far. Both come straight from the torsions, so this is a
//! pure function with no structure knowledge.
//!
//! Runs in `O(1)` time and allocates nothing.

/// The pucker of a five-membered sugar ring.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pucker {
    /// The pseudorotation phase angle, in degrees over `[0, 360)`.
    pub phase_degrees: f64,
    /// The pucker amplitude, in the same angular unit as the input torsions.
    pub amplitude: f64,
}

/// Computes the pseudorotation phase and amplitude from the torsions ν0–ν4.
///
/// The torsions are given in degrees, `nu[k]` being νk; the phase is returned in
/// degrees and the amplitude in degrees. Equal torsions give a phase of zero.
///
/// # Examples
///
/// ```
/// use molframe_analysis::sugar_pucker;
///
/// // Five equal torsions have no preferred direction: the phase is zero and the
/// // amplitude is the central torsion.
/// let pucker = sugar_pucker([5.0, 5.0, 5.0, 5.0, 5.0]);
/// assert!(pucker.phase_degrees.abs() < 1e-9);
/// assert!((pucker.amplitude - 5.0).abs() < 1e-9);
/// ```
#[must_use]
pub fn sugar_pucker(nu: [f64; 5]) -> Pucker {
    let numerator = (nu[4] + nu[1]) - (nu[3] + nu[0]);
    let spread = 36.0_f64.to_radians().sin() + 72.0_f64.to_radians().sin();
    let denominator = 2.0 * nu[2] * spread;
    let phase = numerator.atan2(denominator);
    let mut phase_degrees = phase.to_degrees();
    if phase_degrees < 0.0 {
        phase_degrees += 360.0;
    }
    Pucker {
        phase_degrees,
        amplitude: nu[2] / phase.cos(),
    }
}

#[cfg(test)]
#[path = "pucker_tests.rs"]
mod tests;
