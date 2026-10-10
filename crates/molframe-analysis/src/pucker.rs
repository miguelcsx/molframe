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
/// degrees and the amplitude in degrees. The phase follows the Altona-Sundaralingam
/// tangent form. The amplitude is the least-squares fit of
/// `nu_j = amplitude * cos(P + 144 * (j - 2))` over all five torsions, i.e. the
/// length of the projection onto the two ring-puckering basis vectors, so it
/// stays finite and accurate at every phase, including 90 and 270 degrees where
/// ν2 vanishes. A ring whose torsions are all equal has no puckering component:
/// the amplitude is zero and the phase is reported as zero.
///
/// # Examples
///
/// ```
/// use molframe_analysis::sugar_pucker;
///
/// // Five equal torsions have no puckering component.
/// let pucker = sugar_pucker([5.0, 5.0, 5.0, 5.0, 5.0]);
/// assert!(pucker.phase_degrees.abs() < 1e-9);
/// assert!(pucker.amplitude.abs() < 1e-9);
/// ```
#[must_use]
pub fn sugar_pucker(nu: [f64; 5]) -> Pucker {
    let numerator = (nu[4] + nu[1]) - (nu[3] + nu[0]);
    let spread = 36.0_f64.to_radians().sin() + 72.0_f64.to_radians().sin();
    let denominator = 2.0 * nu[2] * spread;
    let phase = if numerator == 0.0 && denominator == 0.0 {
        0.0
    } else {
        numerator.atan2(denominator)
    };
    let mut phase_degrees = phase.to_degrees();
    if phase_degrees < 0.0 {
        phase_degrees += 360.0;
    }
    // Projections of the torsions onto cos/sin of the 144-degree stepped angles.
    let (mut cos_part, mut sin_part) = (0.0, 0.0);
    for (step, value) in [-2.0_f64, -1.0, 0.0, 1.0, 2.0].into_iter().zip(nu.iter()) {
        let angle = (144.0 * step).to_radians();
        cos_part += value * angle.cos();
        sin_part -= value * angle.sin();
    }
    let amplitude = 0.4 * cos_part.hypot(sin_part);
    Pucker {
        phase_degrees,
        amplitude,
    }
}

#[cfg(test)]
#[path = "pucker_tests.rs"]
mod tests;
