//! How work is divided, independently of how many threads will run it.
//!
//! The whole determinism guarantee rests on one property: block boundaries are
//! a function of the item count and a fixed block size, never of the worker
//! count. Partials are then merged in ascending block index, so the
//! floating-point grouping is identical at one worker and at sixteen.
//!
//! Partitioning by worker count instead — dividing the items into as many
//! contiguous pieces as there are threads — is correct for a pure map, where
//! each item's output is independent and only the concatenation order matters.
//! It is wrong the moment a block accumulates, because the accumulation
//! boundaries then move with the thread count and the sums reassociate. The
//! plan makes the distinction explicit rather than leaving each kernel to
//! rediscover it.

use std::ops::Range;

/// Items per block when a caller expresses no preference.
///
/// Small enough that a block's working set stays cache-resident and that work
/// distributes evenly across threads; large enough that per-block bookkeeping
/// disappears against the work itself.
pub const DEFAULT_BLOCK_ITEMS: usize = 64;

/// A division of `0..count` into fixed-size blocks.
///
/// # Examples
///
/// ```
/// use pdbiox_core::parallel::BlockPlan;
///
/// let plan = BlockPlan::new(10, 4);
/// assert_eq!(plan.blocks(), 3);
/// assert_eq!(plan.range(0), Some(0..4));
/// assert_eq!(plan.range(2), Some(8..10));
/// assert_eq!(plan.range(3), None);
/// ```
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct BlockPlan {
    count: usize,
    block: usize,
}

impl BlockPlan {
    /// Divides `count` items into blocks of at most `block` items.
    ///
    /// A zero block size is raised to one, because a plan that cannot advance
    /// is not a plan.
    #[must_use]
    pub const fn new(count: usize, block: usize) -> Self {
        Self {
            count,
            block: if block == 0 { 1 } else { block },
        }
    }

    /// Divides `count` items using [`DEFAULT_BLOCK_ITEMS`].
    #[must_use]
    pub const fn with_default_block(count: usize) -> Self {
        Self::new(count, DEFAULT_BLOCK_ITEMS)
    }

    /// The number of items covered.
    #[must_use]
    pub const fn count(self) -> usize {
        self.count
    }

    /// The number of items in a full block.
    #[must_use]
    pub const fn block(self) -> usize {
        self.block
    }

    /// The number of blocks.
    #[must_use]
    pub const fn blocks(self) -> usize {
        self.count.div_ceil(self.block)
    }

    /// Returns true when there is nothing to divide.
    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.count == 0
    }

    /// The half-open item range block `index` covers.
    #[must_use]
    pub fn range(self, index: usize) -> Option<Range<usize>> {
        if index >= self.blocks() {
            return None;
        }
        let start = index.checked_mul(self.block)?;
        let end = start.checked_add(self.block)?.min(self.count);
        Some(start..end)
    }

    /// Every block range, in ascending index order.
    pub fn ranges(self) -> impl Iterator<Item = Range<usize>> {
        (0..self.blocks()).filter_map(move |index| self.range(index))
    }

    /// How many workers can be usefully applied to this plan.
    ///
    /// Never more than the number of blocks, because a worker with no block
    /// costs a thread and returns nothing.
    #[must_use]
    pub fn useful_workers(self, requested: usize) -> usize {
        requested.max(1).min(self.blocks().max(1))
    }
}

#[cfg(test)]
#[path = "plan_tests.rs"]
mod tests;
