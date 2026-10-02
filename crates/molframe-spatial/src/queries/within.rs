//! `within` selection reduction: which query atoms lie near a target set.

use super::visit::for_each_pairs_within_unsorted;
use crate::{
    NeighborPair, PairQuery, PeriodicBox, SpatialBackend, SpatialError, SpatialSearchOptions,
};
use molframe_core::{ExecutionContext, selection::AtomSelection};

/// Selects query atoms within `cutoff` of any target atom.
///
/// Target atoms are included when they also belong to `query`, matching the
/// conventional meaning of `within` rather than `around`.
///
/// # Errors
///
/// Returns the same errors as [`crate::pairs_within`].
pub fn within(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    target: &AtomSelection,
    cutoff: f32,
    backend: SpatialBackend,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<AtomSelection, SpatialError> {
    within_with_options(
        positions,
        query,
        target,
        cutoff,
        SpatialSearchOptions::with_backend(backend),
        periodic,
        context,
    )
}

/// Selects query atoms under a complete, inspectable planning profile.
///
/// # Errors
///
/// Returns the same errors as [`crate::pairs_within_with_options`].
pub fn within_with_options(
    positions: &[[f32; 3]],
    query: &AtomSelection,
    target: &AtomSelection,
    cutoff: f32,
    options: SpatialSearchOptions,
    periodic: Option<&PeriodicBox>,
    context: &ExecutionContext,
) -> Result<AtomSelection, SpatialError> {
    let mut matched = Matched::new(positions.len());
    let mut reduction_error = None;
    for_each_pairs_within_unsorted(
        &PairQuery {
            positions,
            left: query,
            right: target,
            cutoff,
            options,
            periodic,
            context,
        },
        |pair| {
            if reduction_error.is_none()
                && let Err(error) = mark_pair_matches(&mut matched, query, target, pair)
            {
                reduction_error = Some(error);
            }
        },
    )?;
    if let Some(error) = reduction_error {
        return Err(error);
    }

    collect_matches(query, target, &matched)
}

/// One bit per atom: which query atoms some pair has already matched.
///
/// Packed words rather than `bool`s, so a 100,000-atom structure needs 12.5 kB
/// instead of 100 kB; plain words rather than a shared bit vector, so marking a
/// bit costs no copy-on-write check.
struct Matched {
    words: Vec<u64>,
    len: usize,
}

impl Matched {
    fn new(len: usize) -> Self {
        Self {
            words: vec![0; len.div_ceil(64)],
            len,
        }
    }

    /// Marks one atom.
    ///
    /// # Errors
    ///
    /// Returns an out-of-bounds error for an invalid pair endpoint.
    fn mark(&mut self, atom: u32) -> Result<(), SpatialError> {
        let index = self.index(atom)?;
        self.words[index / 64] |= 1 << (index % 64);
        Ok(())
    }

    fn get(&self, atom: u32) -> Result<bool, SpatialError> {
        let index = self.index(atom)?;
        Ok(self.words[index / 64] & (1 << (index % 64)) != 0)
    }

    fn index(&self, atom: u32) -> Result<usize, SpatialError> {
        let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
        if index >= self.len {
            return Err(SpatialError::AtomOutOfBounds(atom));
        }
        Ok(index)
    }
}

/// Marks query endpoints participating in query-target neighbour pairs.
///
/// The dense bitmap removes the `O(M log M)` sort/dedup stage previously
/// required by `within`.
fn mark_pair_matches(
    matched: &mut Matched,
    query: &AtomSelection,
    target: &AtomSelection,
    pair: NeighborPair,
) -> Result<(), SpatialError> {
    if query.contains(pair.first) && target.contains(pair.second) {
        matched.mark(pair.first)?;
    }
    if query.contains(pair.second) && target.contains(pair.first) {
        matched.mark(pair.second)?;
    }
    Ok(())
}

/// Collects marked atoms in query-selection order.
///
/// Query-target overlap is included without requiring self-pairs.
fn collect_matches(
    query: &AtomSelection,
    target: &AtomSelection,
    matched: &Matched,
) -> Result<AtomSelection, SpatialError> {
    let mut selected = Vec::new();

    for atom in query {
        if matched.get(atom)? || target.contains(atom) {
            selected.push(atom);
        }
    }

    Ok(AtomSelection::from_sorted(selected))
}
