//! Over-coordination: atoms with more bonds than their element allows.
//!
//! A carbon with five bonds or an oxygen with three is chemically impossible and
//! signals a modelling error — a spurious connection, a misassigned element. This
//! counts each atom's bonds and flags any that exceed the most its element can
//! form. It deliberately does not flag *under*-coordination: a heavy-atom model
//! legitimately omits hydrogens, so too few bonds is the normal case, not a
//! fault. Elements whose coordination is too variable to bound are left unchecked.
//!
//! The check needs a bond table; without one it reports nothing rather than
//! guessing connectivity.

use pdbiox_core::element::Element;
use pdbiox_core::index::AtomIndex;
use pdbiox_core::structure::Structure;

/// An atom bonded more times than its element permits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ValenceError {
    /// The over-coordinated atom.
    pub atom: AtomIndex,
    /// How many bonds it has.
    pub bonds: usize,
    /// The most its element can form.
    pub maximum: usize,
}

/// Flags atoms with more bonds than their element's maximum coordination.
///
/// Returns an empty list when the structure has no bond table. Findings are
/// ordered by atom index.
///
/// Runs in `O(atoms)` time.
#[must_use]
pub fn overvalent_atoms(structure: &Structure) -> Vec<ValenceError> {
    if !structure.data().bonds.is_available() {
        return Vec::new();
    }
    let adjacency = structure.data().bonds.adjacency(structure.atom_count());
    let mut errors = Vec::new();
    for atom in structure.data().atoms() {
        let Some(element) = atom.element() else {
            continue;
        };
        let Some(maximum) = maximum_bonds(element) else {
            continue;
        };
        let bonds = adjacency.neighbours(atom.index()).len();
        if bonds > maximum {
            errors.push(ValenceError {
                atom: atom.index(),
                bonds,
                maximum,
            });
        }
    }
    errors
}

/// The most bonds an element forms, or `None` when it is too variable to bound.
fn maximum_bonds(element: Element) -> Option<usize> {
    match element {
        Element::HYDROGEN => Some(1),
        Element::OXYGEN => Some(2),
        Element::NITROGEN | Element::CARBON => Some(4),
        Element::PHOSPHORUS => Some(5),
        Element::SULFUR => Some(6),
        _ => None,
    }
}

#[cfg(test)]
#[path = "valence_tests.rs"]
mod tests;
