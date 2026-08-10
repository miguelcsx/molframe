//! Bond-length deviation from summed single-bond covalent radii.
//!
//! Several checks share the same question: is the separation of two bonded
//! atoms close to the sum of their covalent radii? The connectivity check over a
//! whole bond table, the ligand-scoped check, and the nucleic backbone check all
//! ask it, so the arithmetic lives here once. The reference is the single-bond
//! covalent radii, which makes it a coarse sanity check on gross errors rather
//! than a per-bond-type restraint comparison.

use pdbiox_core::element::Element;

use crate::numeric::f64_to_f32;

/// A measured bond length set against the sum of the two covalent radii.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct CovalentMeasurement {
    /// The measured separation, in ångström.
    pub observed: f32,
    /// The sum of the two covalent radii, in ångström.
    pub expected: f32,
    /// The signed difference `observed - expected`.
    pub delta: f32,
}

/// The covalent-radius deviation of a bond, when both radii are known.
///
/// Returns `None` when either element has no tabulated covalent radius, so a
/// caller skips the pair rather than comparing against a guessed length.
///
/// Runs in `O(1)` time and allocates no memory.
pub(crate) fn covalent_deviation(
    position_a: [f32; 3],
    element_a: Element,
    position_b: [f32; 3],
    element_b: Element,
) -> Option<CovalentMeasurement> {
    let radius_a = pdbiox_chem::element_properties(element_a)?.covalent_radius?;
    let radius_b = pdbiox_chem::element_properties(element_b)?.covalent_radius?;
    let observed = f64_to_f32(pdbiox_geom::distance(position_a, position_b));
    let expected = radius_a + radius_b;
    Some(CovalentMeasurement {
        observed,
        expected,
        delta: observed - expected,
    })
}
