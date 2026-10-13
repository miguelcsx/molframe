//! Row-block streaming of the distance matrix under an execution budget.

use super::{MatrixError, distance_rows_f32_into, distance_rows_into};
use molframe_core::ExecutionContext;
use std::ops::Range;

/// The consecutive row ranges that cover an `size`-row matrix, `rows` at a time.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RowBlocks {
    size: usize,
    rows: usize,
    next: usize,
}

impl RowBlocks {
    /// Blocks of `rows` rows covering `size` rows; `None` when `rows` is zero.
    #[must_use]
    pub fn new(size: usize, rows: usize) -> Option<Self> {
        (rows > 0).then_some(Self {
            size,
            rows,
            next: 0,
        })
    }
}

impl Iterator for RowBlocks {
    type Item = Range<usize>;

    fn next(&mut self) -> Option<Range<usize>> {
        if self.next >= self.size {
            return None;
        }
        let start = self.next;
        let stop = start.saturating_add(self.rows).min(self.size);
        self.next = stop;
        Some(start..stop)
    }
}

fn reserve_block(
    context: &ExecutionContext,
    rows: usize,
    size: usize,
    element: usize,
) -> Result<molframe_core::MemoryReservation, MatrixError> {
    let bytes = rows.saturating_mul(size).saturating_mul(element);
    context
        .try_reserve(bytes)
        .map_err(|_| MatrixError::MemoryBudget {
            requested: bytes,
            available: context.available_bytes(),
        })
}

/// One block of distance rows as `f64`, admitted against `context`.
///
/// # Errors
///
/// Memory admission or a row range outside the positions.
pub fn distance_rows_with_context(
    positions: &[[f32; 3]],
    rows: Range<usize>,
    context: &ExecutionContext,
) -> Result<Vec<f64>, MatrixError> {
    let size = positions.len();
    let _held = reserve_block(context, rows.len(), size, size_of::<f64>())?;
    let mut block = vec![0.0_f64; rows.len().saturating_mul(size)];
    distance_rows_into(positions, rows, &mut block)?;
    Ok(block)
}

/// [`distance_rows_with_context`] narrowed to `f32`.
///
/// # Errors
///
/// Same as [`distance_rows_with_context`].
pub fn distance_rows_f32_with_context(
    positions: &[[f32; 3]],
    rows: Range<usize>,
    context: &ExecutionContext,
) -> Result<Vec<f32>, MatrixError> {
    let size = positions.len();
    let _held = reserve_block(context, rows.len(), size, size_of::<f32>())?;
    let mut block = vec![0.0_f32; rows.len().saturating_mul(size)];
    distance_rows_f32_into(positions, rows, &mut block)?;
    Ok(block)
}
