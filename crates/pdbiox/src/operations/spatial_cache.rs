//! Reusable coordinate-array spatial indexes for one native plan execution.

use super::requests::ExecutionPlanError;
use super::spatial::SpatialRequest;
use pdbiox_core::selection::AtomSelection;
use pdbiox_spatial::{
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
        targets: Box<[u32]>,
        index: CellList<'a>,
    },
    Kd {
        targets: Box<[u32]>,
        index: KdTree<'a>,
    },
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
    ) -> Result<Vec<NeighborPair>, pdbiox_spatial::SpatialError> {
        let plan = self.options.plan(
            left_count(left)?,
            left_count(right)?,
            self.periodic.is_some(),
            cutoff,
        )?;
        if self.periodic.is_some() {
            return pdbiox_spatial::pairs_within_with_options(
                self.positions,
                left,
                right,
                cutoff,
                self.options,
                self.periodic.as_ref(),
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
            backend => pdbiox_spatial::pairs_within_with_options(
                self.positions,
                left,
                right,
                cutoff,
                SpatialSearchOptions {
                    backend,
                    ..self.options
                },
                None,
            ),
        }
    }

    pub(super) fn within(
        &mut self,
        query: &AtomSelection,
        target: &AtomSelection,
        cutoff: f32,
    ) -> Result<AtomSelection, pdbiox_spatial::SpatialError> {
        let pairs = self.pairs(query, target, cutoff)?;
        let mut matched = vec![false; self.positions.len()];
        for pair in pairs {
            mark_pair(&mut matched, query, target, pair)?;
        }

        let mut selected = Vec::new();
        for atom in query {
            let index = usize::try_from(atom)
                .map_err(|_| pdbiox_spatial::SpatialError::NumericRangeExceeded)?;
            if matched
                .get(index)
                .copied()
                .ok_or(pdbiox_spatial::SpatialError::AtomOutOfBounds(atom))?
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
    ) -> Result<Vec<NeighborPair>, pdbiox_spatial::SpatialError> {
        if let Some(index) = self.cache.iter().find_map(|entry| match entry {
            CachedIndex::Cell {
                cutoff: built_cutoff,
                targets,
                index,
            } if targets.as_ref() == right && *built_cutoff >= cutoff => Some(index),
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
    ) -> Result<Vec<NeighborPair>, pdbiox_spatial::SpatialError> {
        if let Some(index) = self.cache.iter().find_map(|entry| match entry {
            CachedIndex::Kd { targets, index } if targets.as_ref() == right => Some(index),
            CachedIndex::Cell { .. } | CachedIndex::Kd { .. } => None,
        }) {
            return index.pairs(left, cutoff);
        }
        let index =
            KdTree::build_with_options(self.positions, right, None, self.options.kd_periodic)?;
        let pairs = index.pairs(left, cutoff)?;
        self.retain(CachedIndex::Kd {
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

fn left_count(selection: &AtomSelection) -> Result<usize, pdbiox_spatial::SpatialError> {
    usize::try_from(selection.len()).map_err(|_| pdbiox_spatial::SpatialError::NumericRangeExceeded)
}

fn indices(selection: &AtomSelection) -> Vec<u32> {
    selection.iter().collect()
}

fn mark_pair(
    matched: &mut [bool],
    query: &AtomSelection,
    target: &AtomSelection,
    pair: NeighborPair,
) -> Result<(), pdbiox_spatial::SpatialError> {
    if query.contains(pair.first) && target.contains(pair.second) {
        mark(matched, pair.first)?;
    }
    if query.contains(pair.second) && target.contains(pair.first) {
        mark(matched, pair.second)?;
    }
    Ok(())
}

fn mark(matched: &mut [bool], atom: u32) -> Result<(), pdbiox_spatial::SpatialError> {
    let index =
        usize::try_from(atom).map_err(|_| pdbiox_spatial::SpatialError::NumericRangeExceeded)?;
    let slot = matched
        .get_mut(index)
        .ok_or(pdbiox_spatial::SpatialError::AtomOutOfBounds(atom))?;
    *slot = true;
    Ok(())
}
