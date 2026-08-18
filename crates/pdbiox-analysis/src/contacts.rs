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
use pdbiox_spatial::{SpatialBackend, SpatialError, StructureSpatial, pairs_within};

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
    let all = AtomSelection::from_sorted((0..structure.atom_count()).collect());
    atom_contacts_between(structure, &all, &all, cutoff, backend)
}

/// Finds contacts between two arbitrary atom selections.
///
/// The selections are already typed and sorted by the query layer, so this
/// operation does not parse or iterate over Python objects. The spatial
/// backend only visits candidate neighbours and the returned contacts are
/// deterministic by atom index.
///
/// # Errors
///
/// Returns [`SpatialError`] when the cutoff or spatial workload is invalid.
pub fn atom_contacts_between(
    structure: &Structure,
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
) -> Result<Vec<Contact>, SpatialError> {
    let pairs = pairs_within(structure.positions(), left, right, cutoff, backend, None)?;
    Ok(contacts_from_pairs(structure, pairs))
}

/// Finds contacts between two selections using a reusable structure-bound
/// spatial resolver.
///
/// This is the plan-facing form: compatible workloads reuse the resolver's
/// bounded index cache while the contact projection remains the same native
/// kernel as [`atom_contacts_between`].
///
/// # Errors
///
/// Returns a spatial diagnostic when the resolver rejects the workload or a
/// pair cannot be projected to a finite contact.
pub fn atom_contacts_between_with_spatial(
    structure: &Structure,
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    spatial: &StructureSpatial<'_>,
) -> Result<Vec<Contact>, pdbiox_core::diagnostic::Diagnostic> {
    let pairs = spatial.pairs_with_backend(left, right, cutoff, backend)?;
    Ok(contacts_from_pairs(structure, pairs))
}

fn contacts_from_pairs(
    structure: &Structure,
    pairs: Vec<pdbiox_spatial::NeighborPair>,
) -> Vec<Contact> {
    let positions = structure.positions();

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
    contacts
}

#[cfg(test)]
#[path = "contacts_tests.rs"]
mod tests;
