//! Frame-aware reuse of the shared Verlet neighbour-list implementation.

use crate::{Timestep, TrajectoryError};
use pdbiox_core::coords::CoordinateGeneration;
use pdbiox_core::structure::UnitCell;
use pdbiox_spatial::{NeighborList, NeighborPair, PeriodicBox};

/// Observable cost of a frame neighbour-list cache.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NeighborStatistics {
    /// Frames queried.
    pub frames: usize,
    /// Candidate-list rebuilds.
    pub rebuilds: usize,
}

/// Verlet candidates reused across trajectory frames until skin exhaustion.
#[derive(Debug)]
pub struct FrameNeighborList {
    left: Box<[u32]>,
    right: Box<[u32]>,
    cutoff: f32,
    skin: f32,
    list: Option<NeighborList>,
    generation: CoordinateGeneration,
    statistics: NeighborStatistics,
}

impl FrameNeighborList {
    /// Creates an empty cache over explicit atom sets.
    #[must_use]
    pub fn new(
        left: impl Into<Box<[u32]>>,
        right: impl Into<Box<[u32]>>,
        cutoff: f32,
        skin: f32,
    ) -> Self {
        Self {
            left: left.into(),
            right: right.into(),
            cutoff,
            skin,
            list: None,
            generation: CoordinateGeneration::INITIAL,
            statistics: NeighborStatistics::default(),
        }
    }

    /// Filters pairs at the inner cutoff, rebuilding only when required.
    ///
    /// # Errors
    ///
    /// Returns a shared spatial error for invalid cutoffs, selections or cells.
    pub fn pairs(&mut self, timestep: &Timestep) -> Result<Vec<NeighborPair>, TrajectoryError> {
        self.pairs_positions(&timestep.positions, timestep.cell)
    }

    /// Filters pairs from borrowed coordinates without constructing a frame.
    ///
    /// This is the native array-oriented entry point. It lets bindings retain a
    /// caller-owned coordinate buffer for the duration of the spatial query
    /// instead of cloning it into a [`Timestep`].
    ///
    /// # Errors
    ///
    /// Returns a shared spatial error for invalid cutoffs, selections or cells.
    pub fn pairs_positions(
        &mut self,
        positions: &[[f32; 3]],
        cell: Option<UnitCell>,
    ) -> Result<Vec<NeighborPair>, TrajectoryError> {
        let periodic = cell.map(PeriodicBox::from_cell).transpose()?;
        self.statistics.frames += 1;
        let reusable = self
            .list
            .as_ref()
            .is_some_and(|list| list.can_reuse(positions, periodic.as_ref()));
        if !reusable {
            self.generation = self
                .generation
                .next()
                .ok_or(TrajectoryError::CoordinateGenerationExhausted)?;
            self.list = Some(NeighborList::build(
                positions,
                &self.left,
                &self.right,
                self.cutoff,
                self.skin,
                periodic.as_ref(),
                self.generation,
            )?);
            self.statistics.rebuilds += 1;
        }
        match &self.list {
            Some(list) => Ok(list.pairs(positions, self.cutoff, periodic.as_ref())?),
            None => Ok(Vec::new()),
        }
    }

    /// Frames and rebuild count so skin choices remain observable.
    #[must_use]
    pub const fn statistics(&self) -> NeighborStatistics {
        self.statistics
    }
}

#[cfg(test)]
#[path = "neighbors_tests.rs"]
mod tests;
