//! Reusable Verlet neighbour candidates.

use crate::brute::{canonicalise, distance_squared};
use crate::{CellList, NeighborPair, PeriodicBox, SpatialError};
use pdbiox_core::coords::CoordinateGeneration;

/// Candidate pairs built with a displacement skin for reuse across frames.
#[derive(Debug)]
pub struct NeighborList {
    reference: Vec<[f32; 3]>,
    candidates: Vec<(u32, u32)>,
    skin: f32,
    generation: CoordinateGeneration,
}

impl NeighborList {
    /// Builds candidates out to `cutoff + skin`.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid distances or atom indices.
    pub fn build(
        positions: &[[f32; 3]],
        left: &[u32],
        right: &[u32],
        cutoff: f32,
        skin: f32,
        periodic: Option<&PeriodicBox>,
        generation: CoordinateGeneration,
    ) -> Result<Self, SpatialError> {
        if !skin.is_finite() || skin < 0.0 {
            return Err(SpatialError::InvalidCutoff);
        }
        let outer = cutoff + skin;
        if !outer.is_finite() || cutoff < 0.0 {
            return Err(SpatialError::InvalidCutoff);
        }
        let pairs = CellList::build(positions, right, outer, periodic)?.pairs(left, outer)?;
        Ok(Self {
            reference: positions.to_vec(),
            candidates: pairs
                .into_iter()
                .map(|pair| (pair.first, pair.second))
                .collect(),
            skin,
            generation,
        })
    }

    /// The coordinate generation used to build this candidate set.
    #[must_use]
    pub const fn generation(&self) -> CoordinateGeneration {
        self.generation
    }

    /// Whether a structure-local cache entry belongs to `generation`.
    #[must_use]
    pub const fn is_current(&self, generation: CoordinateGeneration) -> bool {
        self.generation.get() == generation.get()
    }

    /// Whether no atom has moved more than half the skin.
    #[must_use]
    pub fn can_reuse(&self, positions: &[[f32; 3]], periodic: Option<&PeriodicBox>) -> bool {
        if positions.len() != self.reference.len() {
            return false;
        }
        let limit = self.skin * self.skin * 0.25;
        self.reference
            .iter()
            .zip(positions)
            .all(|(before, after)| distance_squared(*before, *after, periodic) <= limit)
    }

    /// Filters the candidates at the requested inner cutoff.
    ///
    /// # Errors
    ///
    /// Returns an error when positions moved far enough that an unseen pair
    /// could have entered the cutoff, or when the cutoff is invalid.
    pub fn pairs(
        &self,
        positions: &[[f32; 3]],
        cutoff: f32,
        periodic: Option<&PeriodicBox>,
    ) -> Result<Vec<NeighborPair>, SpatialError> {
        if !cutoff.is_finite() || cutoff < 0.0 {
            return Err(SpatialError::InvalidCutoff);
        }
        if !self.can_reuse(positions, periodic) {
            return Err(SpatialError::StaleNeighborList);
        }
        let cutoff_squared = cutoff * cutoff;
        let mut found = Vec::new();
        for (left, right) in &self.candidates {
            let (Some(left_position), Some(right_position)) = (
                positions.get(*left as usize),
                positions.get(*right as usize),
            ) else {
                return Err(SpatialError::AtomOutOfBounds((*left).max(*right)));
            };
            let squared = distance_squared(*left_position, *right_position, periodic);
            if squared <= cutoff_squared {
                found.push(NeighborPair::new(*left, *right, squared));
            }
        }
        canonicalise(&mut found);
        Ok(found)
    }
}

#[cfg(test)]
#[path = "neighbor_tests.rs"]
mod tests;
