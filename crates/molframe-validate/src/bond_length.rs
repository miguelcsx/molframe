//! Bond lengths that stray from what the two elements expect.
//!
//! A quick, dictionary-free sanity check on connectivity: the length of a bond
//! should be close to the sum of the two atoms' covalent radii. A large
//! departure means a bond that is stretched, compressed, or simply wrong. The
//! reference is the single-bond covalent radii, so the tolerance has to be
//! generous — this catches gross errors, not the fine deviations a restraint
//! library measures against ideal geometry per bond type.
//!
//! Needs a bond table and elements with known covalent radii; bonds missing
//! either are skipped rather than guessed at.

use molframe_core::index::AtomIndex;
use molframe_core::structure::Structure;

/// A bond whose length departs from the sum of covalent radii.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BondDeviation {
    /// One endpoint.
    pub atom_a: AtomIndex,
    /// The other endpoint.
    pub atom_b: AtomIndex,
    /// The measured bond length, in ångström.
    pub observed: f32,
    /// The sum of the two covalent radii.
    pub expected: f32,
    /// The signed difference `observed - expected`.
    pub deviation: f32,
}

/// Flags bonds whose length differs from the expected by more than `tolerance`.
///
/// Findings keep the bond table's order (by endpoint indices). Runs in
/// `O(bonds)` time.
#[must_use]
pub fn bond_length_deviations(structure: &Structure, tolerance: f32) -> Vec<BondDeviation> {
    if !structure.data().bonds.is_available() {
        return Vec::new();
    }
    let positions = structure.positions();
    let mut deviations = Vec::new();
    for bond in structure.data().bonds.iter() {
        let (Some(element_a), Some(element_b)) = (
            element_of(structure, bond.atom_a),
            element_of(structure, bond.atom_b),
        ) else {
            continue;
        };
        let (Some(&a), Some(&b)) = (
            positions.get(bond.atom_a.as_usize()),
            positions.get(bond.atom_b.as_usize()),
        ) else {
            continue;
        };
        let Some(deviation) = crate::covalent::covalent_deviation(a, element_a, b, element_b)
        else {
            continue;
        };
        if deviation.delta.abs() > tolerance {
            deviations.push(BondDeviation {
                atom_a: bond.atom_a,
                atom_b: bond.atom_b,
                observed: deviation.observed,
                expected: deviation.expected,
                deviation: deviation.delta,
            });
        }
    }
    deviations
}

/// The element of an atom by index, if it is recorded.
fn element_of(structure: &Structure, atom: AtomIndex) -> Option<molframe_core::element::Element> {
    structure.data().atom(atom)?.element()
}

#[cfg(test)]
#[path = "bond_length_tests.rs"]
mod tests;
