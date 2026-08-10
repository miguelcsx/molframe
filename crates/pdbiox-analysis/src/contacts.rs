//! Atom pairs closer than a cutoff.
//!
//! A contact is the simplest interaction: two atoms within a distance. What the
//! contact *means* — a bond, a clash, a hydrogen bond — is a question for a
//! layer that knows chemistry; this only reports who is near whom, sorted, so a
//! caller can filter it however the question demands.
//!
//! The pairs come from the shared spatial search, so the cost tracks the local
//! density rather than the square of the atom count.

use pdbiox_core::index::AtomIndex;
use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_spatial::{SpatialBackend, SpatialError, pairs_within};

use crate::numeric::f64_to_f32;

/// Two atoms found within the cutoff, with the distance between them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Contact {
    /// The lower-indexed atom of the pair.
    pub first: AtomIndex,
    /// The higher-indexed atom of the pair.
    pub second: AtomIndex,
    /// The distance between the two atoms.
    pub distance: f32,
}

/// Finds every pair of atoms no further apart than `cutoff`.
///
/// Pairs are unordered and unique, sorted by `(first, second)`, and never
/// include an atom with itself. Only atoms with finite coordinates take part.
///
/// Runs in `O(atoms · local density)` time via the shared spatial search.
///
/// # Errors
///
/// Returns [`SpatialError`] for a non-finite or negative cutoff, propagated from
/// the neighbour search.
pub fn atom_contacts(
    structure: &Structure,
    cutoff: f32,
    backend: SpatialBackend,
) -> Result<Vec<Contact>, SpatialError> {
    let positions = structure.positions();
    let all = AtomSelection::from_sorted((0..structure.atom_count()).collect());
    let pairs = pairs_within(positions, &all, &all, cutoff, backend, None)?;

    let mut contacts = Vec::with_capacity(pairs.len());
    for pair in pairs {
        let first = pair.first as usize;
        let second = pair.second as usize;
        let (Some(&a), Some(&b)) = (positions.get(first), positions.get(second)) else {
            continue;
        };
        contacts.push(Contact {
            first: AtomIndex::new(pair.first),
            second: AtomIndex::new(pair.second),
            distance: f64_to_f32(pdbiox_geom::distance(a, b)),
        });
    }
    Ok(contacts)
}

#[cfg(test)]
#[path = "contacts_tests.rs"]
mod tests;
