//! Bond-length deviation from summed single-bond covalent radii.
//!
//! Several checks share the same question: is the separation of two bonded
//! atoms close to the sum of their covalent radii? The connectivity check over a
//! whole bond table, the ligand-scoped check, and the nucleic backbone check all
//! ask it, so the arithmetic lives here once. The reference is the single-bond
//! covalent radii, which makes it a coarse sanity check on gross errors rather
//! than a per-bond-type restraint comparison.

use molframe_chem::CovalentRadii;
use molframe_core::element::Element;
use molframe_core::index::AtomIndex;
use molframe_core::structure::Structure;

use crate::covalent_length::BondDeviation;
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

impl CovalentMeasurement {
    /// The finding for the bond between `atom_a` and `atom_b` this measures.
    pub(crate) const fn into_deviation(
        self,
        atom_a: AtomIndex,
        atom_b: AtomIndex,
    ) -> BondDeviation {
        BondDeviation {
            atom_a,
            atom_b,
            observed: self.observed,
            expected: self.expected,
            deviation: self.delta,
        }
    }
}

/// The covalent-radius deviation of a bond, when both radii are known.
///
/// Returns `None` when either element has no tabulated covalent radius, so a
/// caller skips the pair rather than comparing against a guessed length.
///
/// `radii` is the covalent-radius table in force.
///
/// Runs in `O(1)` time and allocates no memory.
pub(crate) fn covalent_deviation(
    radii: &CovalentRadii,
    position_a: [f32; 3],
    element_a: Element,
    position_b: [f32; 3],
    element_b: Element,
) -> Option<CovalentMeasurement> {
    let radius_a = radii.radius(element_a)?;
    let radius_b = radii.radius(element_b)?;
    let observed = f64_to_f32(molframe_geom::distance(position_a, position_b));
    let expected = radius_a + radius_b;
    Some(CovalentMeasurement {
        observed,
        expected,
        delta: observed - expected,
    })
}

/// The covalent-radius deviation of the bond between two atoms of `structure`.
///
/// Returns `None` when either atom lacks an element, a tabulated radius or a
/// coordinate, so callers skip the bond rather than guess.
pub(crate) fn bond_deviation(
    structure: &Structure,
    positions: &[[f32; 3]],
    radii: &CovalentRadii,
    atom_a: AtomIndex,
    atom_b: AtomIndex,
) -> Option<CovalentMeasurement> {
    let element_a = structure.data().atom(atom_a)?.element()?;
    let element_b = structure.data().atom(atom_b)?.element()?;
    let a = *positions.get(atom_a.as_usize())?;
    let b = *positions.get(atom_b.as_usize())?;
    covalent_deviation(radii, a, element_a, b, element_b)
}
