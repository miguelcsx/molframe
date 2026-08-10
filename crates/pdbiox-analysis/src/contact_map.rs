//! Which residues touch which, reduced from atom contacts.
//!
//! Two residues are in contact when any of their atoms are within the cutoff.
//! The map is built from the atom-level contacts rather than by scanning every
//! residue pair, so a large structure costs the same neighbour search the atom
//! contacts already pay for, plus a reduction over the pairs it returns.
//!
//! Same-residue pairs are never reported, and a topological separation filter
//! drops residues that sit too close in the chain ordering to count as a
//! long-range contact.

use pdbiox_core::index::ResidueIndex;
use pdbiox_core::structure::Structure;
use pdbiox_spatial::{SpatialBackend, SpatialError};
use std::collections::BTreeMap;

use crate::atom_pairs::atom_contacts;

/// A residue pair in contact and the closest approach between them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResidueContact {
    /// The lower-indexed residue.
    pub first: ResidueIndex,
    /// The higher-indexed residue.
    pub second: ResidueIndex,
    /// The smallest atom-to-atom distance found between the two residues.
    pub min_distance: f32,
}

/// The residue contacts of a structure, sorted by residue pair.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactMap {
    residue_count: usize,
    contacts: Vec<ResidueContact>,
}

impl ContactMap {
    /// The residue pairs in contact, ordered by `(first, second)`.
    #[must_use]
    pub fn contacts(&self) -> &[ResidueContact] {
        &self.contacts
    }

    /// How many residues the structure holds, contact or not.
    #[must_use]
    pub fn residue_count(&self) -> usize {
        self.residue_count
    }
}

/// Builds the residue contact map at a distance cutoff.
///
/// `min_separation` drops any pair whose residue indices differ by less than it,
/// so passing `0` keeps every inter-residue contact and passing `1` keeps them
/// all except adjacent residues. The separation is topological — a difference of
/// residue positions in the structure — so cross-chain pairs, which have no
/// shared numbering, are always kept.
///
/// Runs in the cost of one atom-contact search plus a reduction over its pairs.
///
/// # Errors
///
/// Returns [`SpatialError`] propagated from the neighbour search.
pub fn residue_contact_map(
    structure: &Structure,
    cutoff: f32,
    min_separation: u32,
    backend: SpatialBackend,
) -> Result<ContactMap, SpatialError> {
    let residue_of = atom_to_residue(structure);
    let atom_pairs = atom_contacts(structure, cutoff, backend)?;

    let mut closest: BTreeMap<(u32, u32), f32> = BTreeMap::new();
    for contact in atom_pairs {
        let (Some(&Some(first)), Some(&Some(second))) = (
            residue_of.get(contact.first.as_usize()),
            residue_of.get(contact.second.as_usize()),
        ) else {
            continue;
        };
        if first == second {
            continue;
        }
        let (low, high) = if first <= second {
            (first, second)
        } else {
            (second, first)
        };
        if high - low < min_separation {
            continue;
        }
        closest
            .entry((low, high))
            .and_modify(|distance| {
                if contact.distance < *distance {
                    *distance = contact.distance;
                }
            })
            .or_insert(contact.distance);
    }

    let contacts = closest
        .into_iter()
        .map(|((low, high), min_distance)| ResidueContact {
            first: ResidueIndex::new(low),
            second: ResidueIndex::new(high),
            min_distance,
        })
        .collect();

    Ok(ContactMap {
        residue_count: structure.residue_count(),
        contacts,
    })
}

/// Maps each atom index to the residue that owns it.
///
/// Atoms outside any residue — which a malformed input can produce — map to
/// `None` and take no part in a contact.
fn atom_to_residue(structure: &Structure) -> Vec<Option<u32>> {
    let mut residue_of = vec![None; structure.atom_count() as usize];
    for residue in structure.data().residues() {
        let ordinal = residue.index().get();
        for atom in residue.atoms() {
            if let Some(slot) = residue_of.get_mut(atom.index().as_usize()) {
                *slot = Some(ordinal);
            }
        }
    }
    residue_of
}

#[cfg(test)]
#[path = "contact_map_tests.rs"]
mod tests;
