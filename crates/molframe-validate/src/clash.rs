//! Steric clashes: atoms closer than their van der Waals radii allow.
//!
//! Two non-bonded atoms clash when the sum of their van der Waals radii exceeds
//! the distance between them by more than a tolerance. The tolerance keeps the
//! ordinary give of a contact from being reported as an error while still
//! catching the overlaps that signal a modelling mistake.
//!
//! Directly bonded atoms (1-2) and atoms sharing a bonded neighbour (1-3) are
//! never a clash — their separation is fixed by bond lengths and angles, and
//! their van der Waals spheres overlap by construction — so they are excluded
//! whenever bond information is available. Atoms three bonds apart (1-4) are
//! reported by default; [`ClashOptions::exclude_one_four`] excludes them too.
//! Atoms in different alternate locations never coexist, so such pairs are
//! skipped. The candidate pairs
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

/// Options for [`clashes_with`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClashOptions {
    /// Overlap tolerated before a pair counts as a clash, in ångström.
    pub tolerance: f32,
    /// Van der Waals radius set.
    pub radius_set: RadiusSet,
    /// Neighbour-search backend.
    pub backend: SpatialBackend,
    /// Also exclude atoms exactly three bonds apart. Off by default: 1-4 pairs
    /// can legitimately sit close, but a real overlap there is still an error
    /// for most models.
    pub exclude_one_four: bool,
}

impl ClashOptions {
    /// Default policy (1-2 and 1-3 excluded, 1-4 reported) with the given
    /// tolerance, radius set and backend.
    #[must_use]
    pub const fn new(tolerance: f32, radius_set: RadiusSet, backend: SpatialBackend) -> Self {
        Self {
            tolerance,
            radius_set,
            backend,
            exclude_one_four: false,
        }
    }
}

/// Finds steric clashes under a van der Waals radius set.
///
/// `tolerance` is the overlap tolerated before a pair counts as a clash; a
/// common choice is about 0.4 Å. Atoms whose element has no radius in the set
/// take no part. Bonded (1-2) and angle-related (1-3) pairs and pairs from
/// incompatible alternate locations are skipped; 1-4 pairs are kept (see
/// [`clashes_with`]). Results are sorted by `(first, second)`.
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
    clashes_with(
        structure,
        ClashOptions::new(tolerance, radius_set, backend),
        context,
    )
}

/// Finds steric clashes with explicit exclusion policy.
///
/// # Errors
///
/// Returns [`SpatialError`] from the neighbour search.
pub fn clashes_with(
    structure: &Structure,
    options: ClashOptions,
    context: &ExecutionContext,
) -> Result<ClashTable, SpatialError> {
    let ClashOptions {
        tolerance,
        radius_set,
        backend,
        exclude_one_four,
    } = options;
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
        if bonded_within(bonds, first, second, exclude_one_four) {
            return;
        }
        if let (Some(atom_a), Some(atom_b)) =
            (structure.data().atom(first), structure.data().atom(second))
            && !crate::backbone::alt_compatible(atom_a, atom_b)
        {
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

/// Whether two atoms are 1-2 or 1-3 related (and 1-4 when requested).
fn bonded_within(
    bonds: &molframe_core::bond::BondAdjacency,
    a: AtomIndex,
    b: AtomIndex,
    include_one_four: bool,
) -> bool {
    let near_a = bonds.neighbours(a);
    if near_a.binary_search(&b).is_ok() {
        return true;
    }
    let near_b = bonds.neighbours(b);
    if near_a
        .iter()
        .any(|middle| near_b.binary_search(middle).is_ok())
    {
        return true;
    }
    include_one_four
        && near_a.iter().any(|&x| {
            bonds
                .neighbours(x)
                .iter()
                .any(|&y| y != a && near_b.binary_search(&y).is_ok())
        })
}

#[cfg(test)]
#[path = "clash_tests.rs"]
mod tests;
