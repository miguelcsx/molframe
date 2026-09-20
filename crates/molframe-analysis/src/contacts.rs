//! Atom pairs closer than a cutoff.
//!
//! A contact is the simplest interaction: two atoms within a distance. What the
//! contact *means* — a bond, a clash, a hydrogen bond — is a question for a
//! layer that knows chemistry; this only reports who is near whom, sorted, so a
//! caller can filter it however the question demands.
//!
//! The pairs come from the shared spatial search, so the cost tracks the local
//! density rather than the square of the atom count.

use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use molframe_core::{ExecutionContext, index::AtomIndex};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, StructureSpatial,
    for_each_pairs_within_unsorted,
};

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

/// Native structure-of-arrays storage for an unbounded contact result.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ContactTable {
    first: Vec<AtomIndex>,
    second: Vec<AtomIndex>,
    distance: Vec<f32>,
}

impl ContactTable {
    /// Starts an empty table with capacity for `capacity` rows.
    #[must_use]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            first: Vec::with_capacity(capacity),
            second: Vec::with_capacity(capacity),
            distance: Vec::with_capacity(capacity),
        }
    }

    /// Number of aligned rows.
    #[must_use]
    pub fn len(&self) -> usize {
        debug_assert_eq!(self.first.len(), self.second.len());
        debug_assert_eq!(self.first.len(), self.distance.len());
        self.first.len()
    }

    /// Whether the table has no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.first.is_empty()
    }

    /// Lower atom-index column.
    #[must_use]
    pub fn first(&self) -> &[AtomIndex] {
        &self.first
    }

    /// Higher atom-index column.
    #[must_use]
    pub fn second(&self) -> &[AtomIndex] {
        &self.second
    }

    /// Distance column in ångström.
    #[must_use]
    pub fn distances(&self) -> &[f32] {
        &self.distance
    }

    /// One row reconstructed without allocating.
    #[must_use]
    pub fn row(&self, index: usize) -> Option<Contact> {
        Some(Contact {
            first: *self.first.get(index)?,
            second: *self.second.get(index)?,
            distance: *self.distance.get(index)?,
        })
    }

    /// Rows in deterministic table order.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Contact> + '_ {
        self.first
            .iter()
            .copied()
            .zip(self.second.iter().copied())
            .zip(self.distance.iter().copied())
            .map(|((first, second), distance)| Contact {
                first,
                second,
                distance,
            })
    }

    /// Appends one row while preserving column alignment.
    pub fn push(&mut self, contact: Contact) {
        self.first.push(contact.first);
        self.second.push(contact.second);
        self.distance.push(contact.distance);
    }

    /// Moves a complete canonical reducer block into this table.
    pub fn append(&mut self, other: &mut Self) {
        self.first.append(&mut other.first);
        self.second.append(&mut other.second);
        self.distance.append(&mut other.distance);
    }
}

impl FromIterator<Contact> for ContactTable {
    fn from_iter<T: IntoIterator<Item = Contact>>(iter: T) -> Self {
        let iterator = iter.into_iter();
        let (lower, _) = iterator.size_hint();
        let mut table = Self::with_capacity(lower);
        for contact in iterator {
            table.push(contact);
        }
        table
    }
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
    context: &ExecutionContext,
) -> Result<ContactTable, SpatialError> {
    let all = AtomSelection::All(structure.atom_count());
    collect_contacts(structure, &all, &all, cutoff, backend, context)
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
    context: &ExecutionContext,
) -> Result<ContactTable, SpatialError> {
    collect_contacts(structure, left, right, cutoff, backend, context)
}

/// Visits every unique contact without retaining an output vector.
///
/// Emission order is deterministic for a fixed backend. Callers that require
/// globally sorted output should use [`atom_contacts`], which is the explicit
/// materialising operation.
///
/// # Errors
///
/// Returns [`SpatialError`] when the cutoff, selections or spatial plan are
/// invalid.
pub fn visit_atom_contacts(
    structure: &Structure,
    cutoff: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
    emit: impl FnMut(Contact),
) -> Result<(), SpatialError> {
    let all = AtomSelection::All(structure.atom_count());
    visit_atom_contacts_between(structure, &all, &all, cutoff, backend, context, emit)
}

/// Visits contacts between two selections without retaining them.
///
/// # Errors
///
/// Returns [`SpatialError`] when the cutoff, selections or spatial plan are
/// invalid.
pub fn visit_atom_contacts_between(
    structure: &Structure,
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
    mut emit: impl FnMut(Contact),
) -> Result<(), SpatialError> {
    let positions = structure.positions();
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions,
            left,
            right,
            cutoff,
            options: SpatialSearchOptions::with_backend(backend),
            periodic: None,
            context,
        },
        |pair| {
            if let Some(contact) = contact_from_pair(positions, pair) {
                emit(contact);
            }
        },
    )
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
) -> Result<ContactTable, molframe_core::diagnostic::Diagnostic> {
    let pairs = spatial.pairs_with_backend(left, right, cutoff, backend)?;
    let positions = structure.positions();
    Ok(pairs
        .into_iter()
        .filter_map(|pair| contact_from_pair(positions, pair))
        .collect())
}

fn collect_contacts(
    structure: &Structure,
    left: &AtomSelection,
    right: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<ContactTable, SpatialError> {
    let positions = structure.positions();
    let mut contacts = Vec::new();
    visit_atom_contacts_between(
        structure,
        left,
        right,
        cutoff,
        backend,
        context,
        |contact| contacts.push(contact),
    )?;
    contacts.sort_unstable_by_key(|contact| (contact.first, contact.second));
    debug_assert!(contacts.iter().all(|contact| {
        positions.get(contact.first.get() as usize).is_some()
            && positions.get(contact.second.get() as usize).is_some()
    }));
    Ok(contacts.into_iter().collect())
}

fn contact_from_pair(
    positions: &[[f32; 3]],
    pair: molframe_spatial::NeighborPair,
) -> Option<Contact> {
    let first = pair.first as usize;
    let second = pair.second as usize;
    let (&a, &b) = (positions.get(first)?, positions.get(second)?);
    Some(Contact {
        first: AtomIndex::new(pair.first),
        second: AtomIndex::new(pair.second),
        distance: f64_to_f32(molframe_geom::distance(a, b)),
    })
}

#[cfg(test)]
#[path = "contacts_tests.rs"]
mod tests;
