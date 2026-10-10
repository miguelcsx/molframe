//! Salt bridges between oppositely charged side-chain groups.
//!
//! A salt bridge is a close approach between an atom of a negatively charged
//! group and an atom of a positively charged group. Membership comes from the
//! CCD-derived formal-charge annotation; residue and atom names are never used
//! as a substitute for chemistry annotations.
//!
//! Only anion–cation pairs are searched, so two carboxylates or two amines are
//! never reported. Results are sorted by `(anion, cation)`.
//!
//! By default a pair is skipped when its atoms are covalently bonded (a nitro
//! group's N+ and O- are one group, not two) or lie in the same residue (a
//! zwitterionic amino acid's own termini), since neither is an interaction
//! between separate groups. [`SaltBridgeOptions`] can include either kind.
//! Atoms in incompatible alternate locations are never paired.

use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use molframe_core::{ExecutionContext, index::AtomIndex};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, reduce_pairs_within_unsorted,
};

use molframe_spatial::PeriodicBox;

use crate::hbond::altlocs_compatible;
use crate::numeric::f64_to_f32;

/// Explicit policy for salt-bridge detection.
#[allow(
    clippy::struct_excessive_bools,
    reason = "each flag is an independent, documented switch"
)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SaltBridgeOptions {
    /// Largest anion–cation distance in ångström (about 4 Å typically).
    pub max_distance: f32,
    /// Spatial implementation.
    pub backend: SpatialBackend,
    /// Apply the structure unit cell and minimum-image convention; reported
    /// distances are then the periodic ones.
    pub periodic: bool,
    /// Accept a placeholder unit cell for a periodic request. Placeholder
    /// cells (a unit cube) are refused otherwise.
    pub allow_placeholder_cell: bool,
    /// Keep pairs whose atoms are covalently bonded. Default `false`.
    pub include_bonded_pairs: bool,
    /// Keep pairs whose atoms are in the same residue. Default `false`.
    pub include_same_residue_pairs: bool,
}

impl SaltBridgeOptions {
    /// Non-periodic options with bonded and same-residue pairs excluded.
    #[must_use]
    pub const fn new(max_distance: f32, backend: SpatialBackend) -> Self {
        Self {
            max_distance,
            backend,
            periodic: false,
            allow_placeholder_cell: false,
            include_bonded_pairs: false,
            include_same_residue_pairs: false,
        }
    }
}

/// Salt-bridge detection failure.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum SaltBridgeError {
    /// Periodic geometry was requested without a unit cell.
    #[error("periodic salt bridges require a unit cell")]
    MissingCell,
    /// Periodic geometry was requested with a placeholder unit cell.
    #[error("periodic salt bridges refuse a placeholder unit cell")]
    PlaceholderCell,
    /// Spatial indexing rejected the request.
    #[error(transparent)]
    Spatial(#[from] SpatialError),
}

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

define_soa_table! {
    /// Native columnar storage for salt bridges.
    pub struct SaltBridgeTable for SaltBridge {
        /// Anion atom indices.
        anion: AtomIndex,
        /// Cation atom indices.
        cation: AtomIndex,
        /// Interatomic distances.
        distance: f32,
    }
}

/// Finds non-periodic salt bridges no further apart than `max_distance`
/// (about 4 Å typically), excluding bonded and same-residue pairs.
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
) -> Result<SaltBridgeTable, SpatialError> {
    salt_bridges_with_options(
        structure,
        SaltBridgeOptions::new(max_distance, backend),
        context,
    )
    .map_err(|error| match error {
        SaltBridgeError::Spatial(error) => error,
        SaltBridgeError::MissingCell | SaltBridgeError::PlaceholderCell => {
            SpatialError::InvalidCell
        }
    })
}

/// Finds salt bridges under an explicit [`SaltBridgeOptions`] policy.
///
/// # Errors
///
/// Returns [`SaltBridgeError`] for a missing or placeholder cell on a periodic
/// request, or a spatial-search failure.
pub fn salt_bridges_with_options(
    structure: &Structure,
    options: SaltBridgeOptions,
    context: &ExecutionContext,
) -> Result<SaltBridgeTable, SaltBridgeError> {
    let max_distance = options.max_distance;
    let backend = options.backend;
    let periodic_box = if options.periodic {
        let cell = structure.data().cell.ok_or(SaltBridgeError::MissingCell)?;
        if cell.is_placeholder() && !options.allow_placeholder_cell {
            return Err(SaltBridgeError::PlaceholderCell);
        }
        Some(PeriodicBox::from_cell(cell)?)
    } else {
        None
    };
    let bonds = &structure.data().bonds;
    let adjacency = (!options.include_bonded_pairs && bonds.is_available())
        .then(|| bonds.adjacency(structure.atom_count()));
    let residue_of = |atom: u32| {
        structure
            .data()
            .atom(AtomIndex::new(atom))
            .and_then(molframe_core::structure::AtomRef::residue)
            .map(molframe_core::structure::ResidueRef::index)
    };
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
        periodic: periodic_box.as_ref(),
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
            let bonded = adjacency.is_some_and(|adjacency| {
                adjacency
                    .neighbours(AtomIndex::new(anion))
                    .binary_search(&AtomIndex::new(cation))
                    .is_ok()
            });
            if bonded
                || !altlocs_compatible(structure, anion, cation)
                || (!options.include_same_residue_pairs
                    && residue_of(anion).is_some_and(|r| Some(r) == residue_of(cation)))
            {
                return;
            }
            let (Some(&a), Some(&b)) = (
                positions.get(anion as usize),
                positions.get(cation as usize),
            ) else {
                return;
            };
            result.push(SaltBridge {
                anion: AtomIndex::new(anion),
                cation: AtomIndex::new(cation),
                distance: match periodic_box.as_ref() {
                    Some(periodic) => {
                        let d = periodic.displacement(a, b);
                        f64_to_f32(molframe_geom::distance([0.0; 3], d))
                    }
                    None => f64_to_f32(molframe_geom::distance(a, b)),
                },
            });
        })?;

    let mut result: Vec<SaltBridge> = Vec::new();
    for part in parts {
        result.extend(part);
    }

    result.sort_by_key(|bridge| (bridge.anion.get(), bridge.cation.get()));
    Ok(result.into_iter().collect())
}

#[cfg(test)]
#[path = "salt_bridge_tests.rs"]
mod tests;
