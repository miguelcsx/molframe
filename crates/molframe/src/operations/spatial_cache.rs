//! Reusable coordinate-array spatial indexes for one native plan execution.
//!
//! A cached entry is keyed on the target atom list it was built for, and that
//! list is as long as the selection. Comparing it element by element against
//! every cached entry would make a cache *hit* cost more than rebuilding the
//! index — sixty-four full-length comparisons to avoid one build. Entries
//! therefore carry a length and a fingerprint of their targets, and the
//! element-wise comparison runs only when both agree, which is once.
//!
//! The retained target lists are still proportional to the selection size,
//! so a context holding many large distinct selections holds them all. That
//! bound belongs with the selection representation rather than here: the cache
//! is handed `&[u32]`, and it cannot store a selection more compactly than the
//! caller materialised it.

use super::plan::value::ExecutionPlanError;
use super::spatial::SpatialRequest;
use molframe_core::ExecutionContext;
use molframe_core::selection::AtomSelection;
use molframe_spatial::{
    CellList, KdTree, NeighborPair, PeriodicBox, SpatialBackend, SpatialSearchOptions,
};

const MAX_CACHED_INDICES: usize = 64;

/// One coordinate-array workload and its bounded index cache.
pub(super) struct SpatialContext<'a> {
    slot: usize,
    positions: &'a [[f32; 3]],
    options: SpatialSearchOptions,
    periodic: Option<PeriodicBox>,
    cache: Vec<CachedIndex<'a>>,
}

enum CachedIndex<'a> {
    Cell {
        cutoff: f32,
        key: TargetKey,
        targets: Box<[u32]>,
        index: CellList<'a>,
    },
    Kd {
        key: TargetKey,
        targets: Box<[u32]>,
        index: KdTree<'a>,
    },
}

/// A cheap summary of a target list, compared before the list itself.
///
/// Two different lists can share a key, which is why a match here is a licence
/// to compare rather than a conclusion. What it buys is that the comparison
/// happens once instead of once per cached entry.
#[derive(Clone, Copy, PartialEq, Eq)]
struct TargetKey {
    len: usize,
    fingerprint: u64,
}

impl TargetKey {
    /// The offset basis and prime of FNV-1a, which is fixed so that a key is
    /// the same in every run.
    const BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    fn of(targets: &[u32]) -> Self {
        let mut fingerprint = Self::BASIS;
        for &target in targets {
            fingerprint ^= u64::from(target);
            fingerprint = fingerprint.wrapping_mul(Self::PRIME);
        }
        Self {
            len: targets.len(),
            fingerprint,
        }
    }
}

impl<'a> SpatialContext<'a> {
    pub(super) fn new(
        slot: usize,
        positions: &'a [[f32; 3]],
        options: SpatialSearchOptions,
        periodic: Option<PeriodicBox>,
    ) -> Result<Self, ExecutionPlanError> {
        options
            .validate()
            .map_err(|error| ExecutionPlanError::SpatialKernel(error.to_string().into()))?;
        Ok(Self {
            slot,
            positions,
            options,
            periodic,
            cache: Vec::new(),
        })
    }

    pub(super) fn matches(&self, request: &SpatialRequest) -> bool {
        match request {
            SpatialRequest::NeighborPairs {
                positions,
                options,
                periodic,
                ..
            }
            | SpatialRequest::AtomsWithin {
                positions,
                options,
                periodic,
                ..
            } => *positions == self.slot && *options == self.options && *periodic == self.periodic,
        }
    }

    pub(super) fn cached_index_count(&self) -> usize {
        self.cache.len()
    }

    pub(super) fn pairs(
        &mut self,
        left: &AtomSelection,
        right: &AtomSelection,
        cutoff: f32,
        context: &ExecutionContext,
    ) -> Result<Vec<NeighborPair>, molframe_spatial::SpatialError> {
        let plan = self.options.plan(
            left_count(left)?,
            left_count(right)?,
            self.periodic.is_some(),
            cutoff,
        )?;
        if self.periodic.is_some() {
            return molframe_spatial::pairs_within_with_options(
                self.positions,
                left,
                right,
                cutoff,
                self.options,
                self.periodic.as_ref(),
                context,
            );
        }

        match plan.backend {
            SpatialBackend::CellList => {
                let left_indices = indices(left);
                let right_indices = indices(right);
                self.cell_pairs(&left_indices, &right_indices, cutoff)
            }
            SpatialBackend::KdTree => {
                let left_indices = indices(left);
                let right_indices = indices(right);
                self.kd_pairs(&left_indices, &right_indices, cutoff)
            }
            backend => molframe_spatial::pairs_within_with_options(
                self.positions,
                left,
                right,
                cutoff,
                SpatialSearchOptions {
                    backend,
                    ..self.options
                },
                None,
                context,
            ),
        }
    }

    pub(super) fn within(
        &mut self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
        context: &ExecutionContext,
    ) -> Result<AtomSelection, molframe_spatial::SpatialError> {
        let pairs = self.pairs(query, target, cutoff, context)?;
        let mut matched = vec![false; self.positions.len()];
        for pair in pairs {
            mark_pair(&mut matched, query, target, pair)?;
        }

        let mut selected = Vec::new();
        for atom in query {
            let index = usize::try_from(atom)
                .map_err(|_| molframe_spatial::SpatialError::NumericRangeExceeded)?;
            if matched
                .get(index)
                .copied()
                .ok_or(molframe_spatial::SpatialError::AtomOutOfBounds(atom))?
                || target.contains(atom)
            {
                selected.push(atom);
            }
        }
        Ok(AtomSelection::from_sorted(selected))
    }

    fn cell_pairs(
        &mut self,
        left: &[u32],
        right: &[u32],
        cutoff: f32,
    ) -> Result<Vec<NeighborPair>, molframe_spatial::SpatialError> {
        let key = TargetKey::of(right);
        if let Some(index) = self.cache.iter().find_map(|entry| match entry {
            CachedIndex::Cell {
                cutoff: built_cutoff,
                key: built_key,
                targets,
                index,
            } if *built_key == key && *built_cutoff >= cutoff && targets.as_ref() == right => {
                Some(index)
            }
            CachedIndex::Cell { .. } | CachedIndex::Kd { .. } => None,
        }) {
            return index.pairs(left, cutoff);
        }
        let index = CellList::build_with_options(
            self.positions,
            right,
            cutoff,
            None,
            self.options.cell_grid,
        )?;
        let pairs = index.pairs(left, cutoff)?;
        self.retain(CachedIndex::Cell {
            cutoff,
            key,
            targets: right.to_vec().into_boxed_slice(),
            index,
        });
        Ok(pairs)
    }

    fn kd_pairs(
        &mut self,
        left: &[u32],
        right: &[u32],
        cutoff: f32,
    ) -> Result<Vec<NeighborPair>, molframe_spatial::SpatialError> {
        let key = TargetKey::of(right);
        if let Some(index) = self.cache.iter().find_map(|entry| match entry {
            CachedIndex::Kd {
                key: built_key,
                targets,
                index,
            } if *built_key == key && targets.as_ref() == right => Some(index),
            CachedIndex::Cell { .. } | CachedIndex::Kd { .. } => None,
        }) {
            return index.pairs(left, cutoff);
        }
        let index =
            KdTree::build_with_options(self.positions, right, None, self.options.kd_periodic)?;
        let pairs = index.pairs(left, cutoff)?;
        self.retain(CachedIndex::Kd {
            key,
            targets: right.to_vec().into_boxed_slice(),
            index,
        });
        Ok(pairs)
    }

    fn retain(&mut self, index: CachedIndex<'a>) {
        if self.cache.len() < MAX_CACHED_INDICES {
            self.cache.push(index);
        }
    }
}

fn left_count(selection: &AtomSelection) -> Result<usize, molframe_spatial::SpatialError> {
    usize::try_from(selection.len())
        .map_err(|_| molframe_spatial::SpatialError::NumericRangeExceeded)
}

fn indices(selection: &AtomSelection) -> Vec<u32> {
    selection.iter().collect()
}

fn mark_pair(
    matched: &mut [bool],
    query: &AtomSelection,
    target: &AtomSelection,
    pair: NeighborPair,
) -> Result<(), molframe_spatial::SpatialError> {
    if query.contains(pair.first) && target.contains(pair.second) {
        mark(matched, pair.first)?;
    }
    if query.contains(pair.second) && target.contains(pair.first) {
        mark(matched, pair.second)?;
    }
    Ok(())
}

fn mark(matched: &mut [bool], atom: u32) -> Result<(), molframe_spatial::SpatialError> {
    let index =
        usize::try_from(atom).map_err(|_| molframe_spatial::SpatialError::NumericRangeExceeded)?;
    let slot = matched
        .get_mut(index)
        .ok_or(molframe_spatial::SpatialError::AtomOutOfBounds(atom))?;
    *slot = true;
    Ok(())
}
