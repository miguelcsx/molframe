//! Reusable Verlet neighbour candidates.

use crate::brute::distance_squared;
use crate::{CellList, NeighborListOptions, NeighborPair, PeriodicBox, SpatialError};
use molframe_core::coords::CoordinateGeneration;

const SAFE_DISPLACEMENT_FRACTION: f32 = 0.5;

/// Candidate pairs built with a displacement skin for reuse across frames.
#[derive(Debug)]
pub struct NeighborList {
    reference: Vec<[f32; 3]>,
    candidates: Vec<(u32, u32)>,
    cutoff: f32,
    skin: f32,
    periodic: Option<PeriodicBox>,
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
        Self::build_with_options(
            positions,
            left,
            right,
            cutoff,
            periodic,
            generation,
            NeighborListOptions::with_skin(skin),
        )
    }

    /// Builds candidates under an explicit cell-grid memory policy.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid options, distances or atom indices.
    pub fn build_with_options(
        positions: &[[f32; 3]],
        left: &[u32],
        right: &[u32],
        cutoff: f32,
        periodic: Option<&PeriodicBox>,
        generation: CoordinateGeneration,
        options: NeighborListOptions,
    ) -> Result<Self, SpatialError> {
        validate_distance(cutoff)?;
        options.validate()?;

        let outer = cutoff + options.skin;

        if !outer.is_finite() {
            return Err(SpatialError::InvalidCutoff);
        }

        let index = CellList::build_with_options(
            positions,
            right,
            outer,
            periodic.copied(),
            options.cell_grid,
        )?;
        let pairs = if left == right {
            crate::backends::cell::pairs_same_selection(&index, left, outer)?
        } else {
            index.pairs(left, outer)?
        };

        Ok(Self {
            reference: positions.to_vec(),
            candidates: pairs
                .into_iter()
                .map(|pair| (pair.first, pair.second))
                .collect(),
            cutoff,
            skin: options.skin,
            periodic: periodic.copied(),
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
    ///
    /// Reuse additionally requires the same coordinate count and the same
    /// periodic box used to construct the candidate set.
    #[must_use]
    pub fn can_reuse(&self, positions: &[[f32; 3]], periodic: Option<&PeriodicBox>) -> bool {
        if positions.len() != self.reference.len()
            || !same_periodic_box(self.periodic.as_ref(), periodic)
        {
            return false;
        }

        let limit = self.skin * SAFE_DISPLACEMENT_FRACTION;
        let limit_squared = limit * limit;

        self.reference.iter().zip(positions).all(|(before, after)| {
            displacement_within_limit(*before, *after, limit_squared, periodic)
        })
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
        validate_query_cutoff(cutoff, self.cutoff)?;

        if !self.can_reuse(positions, periodic) {
            return Err(SpatialError::StaleNeighborList);
        }

        let cutoff_squared = cutoff * cutoff;
        let mut found = Vec::new();

        for &(left, right) in &self.candidates {
            append_candidate_if_within(
                positions,
                left,
                right,
                cutoff_squared,
                periodic,
                &mut found,
            )?;
        }

        // `candidates` originates from canonicalised CellList output.
        // Filtering preserves both its ordering and uniqueness.
        Ok(found)
    }
}

/// Adds one cached candidate when its current distance remains within cutoff.
///
/// # Errors
///
/// Returns an out-of-bounds error if a cached candidate cannot be resolved in
/// the supplied coordinate frame.
fn append_candidate_if_within(
    positions: &[[f32; 3]],
    left: u32,
    right: u32,
    cutoff_squared: f32,
    periodic: Option<&PeriodicBox>,
    found: &mut Vec<NeighborPair>,
) -> Result<(), SpatialError> {
    let left_position = candidate_position(positions, left)?;
    let right_position = candidate_position(positions, right)?;

    let squared = distance_squared(left_position, right_position, periodic);

    if squared <= cutoff_squared {
        found.push(NeighborPair::new(left, right, squared));
    }

    Ok(())
}

/// Returns a cached candidate position or an out-of-bounds error.
///
/// Runtime and auxiliary space are `O(1)`.
fn candidate_position(positions: &[[f32; 3]], atom: u32) -> Result<[f32; 3], SpatialError> {
    let index = usize::try_from(atom).map_err(|_| SpatialError::NumericRangeExceeded)?;
    positions
        .get(index)
        .copied()
        .ok_or(SpatialError::AtomOutOfBounds(atom))
}

/// Tests whether one atom moved within the permitted squared displacement.
///
/// Runtime and auxiliary space are `O(1)`.
#[inline]
fn displacement_within_limit(
    before: [f32; 3],
    after: [f32; 3],
    limit_squared: f32,
    periodic: Option<&PeriodicBox>,
) -> bool {
    distance_squared(before, after, periodic) <= limit_squared
}

/// Compares the periodic cell used to build and query the neighbour list.
///
/// Value equality is used rather than pointer identity.
#[inline]
fn same_periodic_box(built: Option<&PeriodicBox>, current: Option<&PeriodicBox>) -> bool {
    built == current
}

/// Validates a non-negative finite distance.
///
/// Runtime and auxiliary space are `O(1)`.
fn validate_distance(value: f32) -> Result<(), SpatialError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

/// Validates that a requested cutoff is covered by the candidate build cutoff.
///
/// Restricting queries to the original inner cutoff is required for the
/// half-skin Verlet reuse guarantee.
fn validate_query_cutoff(cutoff: f32, build_cutoff: f32) -> Result<(), SpatialError> {
    if cutoff.is_finite() && cutoff >= 0.0 && cutoff <= build_cutoff {
        Ok(())
    } else {
        Err(SpatialError::InvalidCutoff)
    }
}

#[cfg(test)]
#[path = "neighbor_tests.rs"]
mod tests;
