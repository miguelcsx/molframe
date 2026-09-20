//! Steric clashes: atoms closer than their van der Waals radii allow.
//!
//! Two non-bonded atoms clash when the sum of their van der Waals radii exceeds
//! the distance between them by more than a tolerance. The tolerance keeps the
//! ordinary give of a contact from being reported as an error while still
//! catching the overlaps that signal a modelling mistake.
//!
//! Directly bonded atoms are never a clash — they are supposed to overlap — so
//! they are excluded whenever bond information is available. The candidate pairs
//! come from the shared spatial search bounded by the widest radius in play.

use molframe_chem::{RadiusSet, vdw_radius};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::Structure;
use molframe_core::{ExecutionContext, index::AtomIndex};
use molframe_spatial::{
    PairQuery, SpatialBackend, SpatialError, SpatialSearchOptions, reduce_pairs_within_unsorted,
};

use crate::numeric::f64_to_f32;

/// Two atoms overlapping more than the tolerance allows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Clash {
    /// The lower-indexed atom.
    pub first: AtomIndex,
    /// The higher-indexed atom.
    pub second: AtomIndex,
    /// How far the van der Waals spheres interpenetrate, in ångström.
    pub overlap: f32,
}

/// Native structure-of-arrays storage for an unbounded clash result.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ClashTable {
    first: Vec<AtomIndex>,
    second: Vec<AtomIndex>,
    overlap: Vec<f32>,
}

impl ClashTable {
    /// Number of aligned rows.
    #[must_use]
    pub fn len(&self) -> usize {
        debug_assert_eq!(self.first.len(), self.second.len());
        debug_assert_eq!(self.first.len(), self.overlap.len());
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

    /// Van der Waals overlap column.
    #[must_use]
    pub fn overlaps(&self) -> &[f32] {
        &self.overlap
    }

    /// Rows in deterministic table order.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Clash> + '_ {
        self.first
            .iter()
            .copied()
            .zip(self.second.iter().copied())
            .zip(self.overlap.iter().copied())
            .map(|((first, second), overlap)| Clash {
                first,
                second,
                overlap,
            })
    }

    /// Appends one row while preserving column alignment.
    pub fn push(&mut self, clash: Clash) {
        self.first.push(clash.first);
        self.second.push(clash.second);
        self.overlap.push(clash.overlap);
    }

    /// Moves a complete reducer block into this table.
    pub fn append(&mut self, other: &mut Self) {
        self.first.append(&mut other.first);
        self.second.append(&mut other.second);
        self.overlap.append(&mut other.overlap);
    }
}

impl FromIterator<Clash> for ClashTable {
    fn from_iter<T: IntoIterator<Item = Clash>>(rows: T) -> Self {
        let iterator = rows.into_iter();
        let (lower, _) = iterator.size_hint();
        let mut table = Self {
            first: Vec::with_capacity(lower),
            second: Vec::with_capacity(lower),
            overlap: Vec::with_capacity(lower),
        };
        for row in iterator {
            table.push(row);
        }
        table
    }
}

/// Finds steric clashes under a van der Waals radius set.
///
/// `tolerance` is the overlap tolerated before a pair counts as a clash; a
/// common choice is about 0.4 Å. Atoms whose element has no radius in the set
/// take no part. Results are sorted by `(first, second)`.
///
/// Runs in `O(atoms · local density)` time.
///
/// # Errors
///
/// Returns [`SpatialError`] from the neighbour search.
pub fn clashes(
    structure: &Structure,
    tolerance: f32,
    radius_set: RadiusSet,
    backend: SpatialBackend,
    context: &ExecutionContext,
) -> Result<ClashTable, SpatialError> {
    let count = structure.atom_count() as usize;
    let mut radii = vec![f32::NAN; count];
    let mut widest = 0.0f32;
    for atom in structure.data().atoms() {
        let Some(element) = atom.element() else {
            continue;
        };
        if let Some(radius) = vdw_radius(element, radius_set) {
            if let Some(slot) = radii.get_mut(atom.index().as_usize()) {
                *slot = radius;
            }
            widest = widest.max(radius);
        }
    }

    let positions = structure.positions();
    let all = AtomSelection::All(structure.atom_count());
    let cutoff = 2.0 * widest - tolerance;
    if cutoff <= 0.0 {
        return Ok(ClashTable::default());
    }
    let bonds = structure.data().bonds.adjacency(structure.atom_count());

    // Clashes are rare relative to candidate pairs, so pairs are reduced as they
    // are produced rather than collected: retained bytes track the clash count,
    // not the quadratic candidate count.
    let query = PairQuery {
        positions,
        left: &all,
        right: &all,
        cutoff,
        options: SpatialSearchOptions::with_backend(backend),
        periodic: None,
        context,
    };
    let parts = reduce_pairs_within_unsorted(&query, Vec::new, |result: &mut Vec<Clash>, pair| {
        let first = AtomIndex::new(pair.first);
        let second = AtomIndex::new(pair.second);
        let (Some(&radius_a), Some(&radius_b)) = (
            radii.get(pair.first as usize),
            radii.get(pair.second as usize),
        ) else {
            return;
        };
        if radius_a.is_nan() || radius_b.is_nan() {
            return;
        }
        let (Some(&a), Some(&b)) = (
            positions.get(pair.first as usize),
            positions.get(pair.second as usize),
        ) else {
            return;
        };
        let overlap = (radius_a + radius_b) - f64_to_f32(molframe_geom::distance(a, b));
        if overlap <= tolerance {
            return;
        }
        if bonds.neighbours(first).binary_search(&second).is_ok() {
            return;
        }
        result.push(Clash {
            first,
            second,
            overlap,
        });
    })?;

    let mut result: Vec<Clash> = Vec::new();
    for part in parts {
        result.extend(part);
    }

    // The streaming query does not order its pairs, so the clashes are ordered
    // here. Sorting the survivors is far cheaper than materialising and sorting
    // every candidate.
    result.sort_by_key(|clash| (clash.first.get(), clash.second.get()));
    Ok(result.into_iter().collect())
}

#[cfg(test)]
#[path = "clash_tests.rs"]
mod tests;
