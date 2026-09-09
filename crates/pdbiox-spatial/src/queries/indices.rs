//! One charged index expansion shared by spatial dispatch and reductions.

use super::search::checked_indices;
use crate::SpatialError;
use pdbiox_core::{AtomSelection, ExecutionContext, MemoryReservation};
use std::borrow::Cow;

pub(super) struct IndexWorkspace<'a> {
    left: Cow<'a, [u32]>,
    right: Option<Cow<'a, [u32]>>,
    _reservation: MemoryReservation,
}

impl<'a> IndexWorkspace<'a> {
    pub(super) fn new(
        left: &'a AtomSelection,
        right: &'a AtomSelection,
        positions: usize,
        context: &ExecutionContext,
    ) -> Result<Self, SpatialError> {
        if context.cancellation().is_cancelled() {
            return Err(SpatialError::Cancelled);
        }
        let same = std::ptr::eq(left, right) || left == right;
        let bytes = index_bytes(left)?
            .checked_add(if same { 0 } else { index_bytes(right)? })
            .ok_or(SpatialError::NumericRangeExceeded)?;
        let reservation = context.try_reserve(bytes)?;
        let left = checked_indices(left, positions)?;
        let right = if same {
            None
        } else {
            Some(checked_indices(right, positions)?)
        };
        Ok(Self {
            left,
            right,
            _reservation: reservation,
        })
    }

    pub(super) fn left(&self) -> &[u32] {
        &self.left
    }

    pub(super) fn right(&self) -> &[u32] {
        self.right.as_deref().map_or(&self.left, |right| right)
    }
}

fn index_bytes(selection: &AtomSelection) -> Result<usize, SpatialError> {
    if matches!(selection, AtomSelection::Sparse(_) | AtomSelection::Empty) {
        return Ok(0);
    }
    usize::try_from(selection.len())
        .ok()
        .and_then(|length| length.checked_mul(size_of::<u32>()))
        .ok_or(SpatialError::NumericRangeExceeded)
}

#[cfg(test)]
#[path = "indices_tests.rs"]
mod tests;
