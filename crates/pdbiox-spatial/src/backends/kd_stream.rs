//! Charged tree construction and radius visitation with O(N) retained storage.

use super::{
    Entry, KdTree, Node, finite, for_each_shift, periodic_image_limits, validate_image_budget,
};
use crate::{KdPeriodicOptions, PeriodicBox, SpatialError};
use pdbiox_core::ExecutionContext;

impl<'a> KdTree<'a> {
    pub(crate) fn build_in(
        positions: &'a [[f32; 3]],
        targets: &[u32],
        periodic: Option<&'a PeriodicBox>,
        options: KdPeriodicOptions,
        context: &ExecutionContext,
    ) -> Result<Self, SpatialError> {
        if context.cancellation().is_cancelled() {
            return Err(SpatialError::Cancelled);
        }
        super::validate_indices(targets, positions.len())?;
        options.validate()?;
        let bytes = targets
            .len()
            .checked_mul(size_of::<Entry>() + size_of::<Node>())
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let mut reservation = context.try_reserve(bytes)?;
        let mut tree = Self::build_with_options(positions, targets, periodic, options)?;
        reservation.shrink_to(tree.nodes.capacity() * size_of::<Node>());
        tree.reservation = Some(reservation);
        Ok(tree)
    }

    // A generation per tree node suppresses repeated periodic images without a
    // pair vector or a full bitmap clear per query atom. No allocation per atom.
    pub(crate) fn for_each_candidate(
        &self,
        query: &[u32],
        cutoff: f32,
        context: &ExecutionContext,
        mut emit: impl FnMut(u32, u32, f32),
    ) -> Result<(), SpatialError> {
        super::validate_cutoff(cutoff)?;
        super::validate_indices(query, self.positions.len())?;
        let Some(root) = self.root else {
            return Ok(());
        };
        let bytes = if self.periodic.is_some() {
            self.nodes
                .len()
                .checked_mul(size_of::<usize>())
                .ok_or(SpatialError::NumericRangeExceeded)?
        } else {
            0
        };
        let _reservation = context.try_reserve(bytes)?;
        let mut seen = vec![0_usize; bytes / size_of::<usize>()];
        let squared_cutoff = cutoff * cutoff;
        let limits = self
            .periodic
            .map(|periodic| periodic_image_limits(periodic, cutoff))
            .transpose()?;
        if let Some(limits) = limits {
            validate_image_budget(limits, self.periodic_options)?;
        }
        for (ordinal, &atom) in query.iter().enumerate() {
            if context.cancellation().is_cancelled() {
                return Err(SpatialError::Cancelled);
            }
            let position = self.positions[atom as usize];
            if !finite(position) {
                continue;
            }
            let Some(periodic) = self.periodic else {
                self.radius_search(root, position, squared_cutoff, &mut |_, target, squared| {
                    if atom != target {
                        emit(atom, target, squared);
                    }
                });
                continue;
            };
            let epoch = ordinal
                .checked_add(1)
                .ok_or(SpatialError::NumericRangeExceeded)?;
            let wrapped = periodic.wrap(position);
            let Some(limits) = limits else {
                continue;
            };
            for_each_shift(limits, |shift| {
                let image = periodic.translated(wrapped, shift)?;
                self.radius_search(root, image, squared_cutoff, &mut |node, target, _| {
                    if atom == target || seen[node] == epoch {
                        return;
                    }
                    seen[node] = epoch;
                    let squared =
                        periodic.distance_squared(position, self.positions[target as usize]);
                    if squared <= squared_cutoff {
                        emit(atom, target, squared);
                    }
                });
                Ok(())
            })?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "kd_stream_tests.rs"]
mod tests;
