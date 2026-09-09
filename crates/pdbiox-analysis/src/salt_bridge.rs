//! Salt bridges between oppositely charged side-chain groups.
//!
//! A salt bridge is a close approach between an atom of a negatively charged
//! group and an atom of a positively charged group. Membership comes from the
//! CCD-derived formal-charge annotation; residue and atom names are never used
//! as a substitute for chemistry annotations.
//!
//! Only anion–cation pairs are searched, so two carboxylates or two amines are
//! never reported. Results are sorted by `(anion, cation)`.

use pdbiox_core::selection::AtomSelection;
use pdbiox_core::structure::Structure;
use pdbiox_core::{ExecutionContext, index::AtomIndex};
use pdbiox_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, reduce_pairs_within_unsorted,
};

use crate::numeric::f64_to_f32;

/// A charged pair within salt-bridge range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaltBridge {
    /// The negatively charged atom.
    pub anion: AtomIndex,
    /// The positively charged atom.
    pub cation: AtomIndex,
    /// The distance between them.
    pub distance: f32,
}

/// Finds salt bridges no further apart than `max_distance` (about 4 Å typically).
///
/// Runs in `O(charged atoms · local density)` time.
///
/// # Errors
///
/// Returns [`SpatialError`] for a non-finite or negative cutoff.
pub fn salt_bridges(
    structure: &Structure,
    max_distance: f32,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<Vec<SaltBridge>, SpatialError> {
    let positions = structure.positions();
    let mut anions = Vec::new();
    let mut cations = Vec::new();
    for atom in structure.data().atoms() {
        match crate::chemistry::formal_charge(structure, atom.index().get()) {
            Some(charge) if charge < 0 => anions.push(atom.index().get()),
            Some(charge) if charge > 0 => cations.push(atom.index().get()),
            _ => {}
        }
    }
    anions.sort_unstable();
    cations.sort_unstable();

    let anion_set = AtomSelection::from_sorted(anions);
    let cation_set = AtomSelection::from_sorted(cations);
    // The reduction keeps only the bridges, so pairs are visited as they are
    // produced rather than collected into a vector sized by the candidate count.
    let query = PairQuery {
        positions,
        left: &anion_set,
        right: &cation_set,
        cutoff: max_distance,
        options: SpatialSearchOptions::with_backend(backend),
        periodic: None,
        context,
    };
    let parts =
        reduce_pairs_within_unsorted(&query, Vec::new, |result: &mut Vec<SaltBridge>, pair| {
            // A pair joins one atom from each set; whichever is the anion is
            // `first` only when it happened to have the lower index, so classify
            // explicitly.
            let (anion, cation) = if anion_set.contains(pair.first) {
                (pair.first, pair.second)
            } else {
                (pair.second, pair.first)
            };
            let (Some(&a), Some(&b)) = (
                positions.get(anion as usize),
                positions.get(cation as usize),
            ) else {
                return;
            };
            result.push(SaltBridge {
                anion: AtomIndex::new(anion),
                cation: AtomIndex::new(cation),
                distance: f64_to_f32(pdbiox_geom::distance(a, b)),
            });
        })?;

    let mut result: Vec<SaltBridge> = Vec::new();
    for part in parts {
        result.extend(part);
    }

    result.sort_by_key(|bridge| (bridge.anion.get(), bridge.cation.get()));
    Ok(result)
}

#[cfg(test)]
#[path = "salt_bridge_tests.rs"]
mod tests;
