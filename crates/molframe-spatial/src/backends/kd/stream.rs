//! Charged tree construction and serial radius visitation.

use super::KdTree;
use crate::{KdPeriodicOptions, PeriodicBox, SpatialError};
use molframe_core::ExecutionContext;

impl<'a> KdTree<'a> {
    /// Builds a tree whose peak and retained bytes are charged to `context`.
    ///
    /// Construction partitions on the context's workers; the peak reservation
    /// covers the entry scratch and the finished arrays, and only the
    /// retained arrays stay charged for the life of the tree.
    pub(crate) fn build_in(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        periodic: Option<PeriodicBox>,
        options: KdPeriodicOptions,
        context: &ExecutionContext,
    ) -> Result<Self, SpatialError> {
        if context.cancellation().is_cancelled() {
            return Err(SpatialError::Cancelled);
        }
        super::validate_indices(targets, positions.len())?;
        options.validate()?;
        let bytes = Self::build_bytes(targets.len()).ok_or(SpatialError::NumericRangeExceeded)?;
        let mut reservation = context.try_reserve(bytes)?;
        let mut tree = Self::build_parallel(positions, targets, periodic, options, Some(context))?;
        reservation.shrink_to(tree.retained_bytes());
        tree.reservation = Some(reservation);
        Ok(tree)
    }

    /// Visits every target within `cutoff` of each query atom, serially.
    ///
    /// Periodic targets are emitted once per query with their exact
    /// minimum-image distance. Nothing is allocated per query.
    pub(crate) fn for_each_candidate(
        &self,
        query: &[u32],
        cutoff: f32,
        context: &ExecutionContext,
        mut emit: impl FnMut(u32, u32, f32),
    ) -> Result<(), SpatialError> {
        super::validate_cutoff(cutoff)?;
        super::validate_indices(query, self.positions.len())?;
        let mut collector = self.radius_collector(cutoff)?;
        for &atom in query {
            if context.cancellation().is_cancelled() {
                return Err(SpatialError::Cancelled);
            }
            collector.collect(self, atom, |target, squared| emit(atom, target, squared))?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
