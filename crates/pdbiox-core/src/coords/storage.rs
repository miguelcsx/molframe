//! Coordinate storage.
//!
//! Positions are `f32` triples in one contiguous, 64-byte-aligned buffer per
//! frame. Three decimal places in ångström is what gets deposited, and `f32`
//! represents that exactly over any box that will ever be modelled, so the
//! wider type would double the bandwidth of every geometric kernel to store
//! digits nobody measured. Kernels that accumulate — centroids, inertia
//! tensors, sums of squares — widen to `f64` internally, where it matters.
//!
//! One buffer per frame, rather than one per chunk, is what makes handing the
//! whole coordinate set to a numeric array a pointer rather than a copy. Chunks
//! address ranges of it.
//!
//! Alignment is obtained by allocating in lanes of sixteen positions. A lane is
//! 192 bytes and demands 64-byte alignment, so the allocator hands back an
//! aligned block and no pointer arithmetic is needed to get one. The cost is at
//! most fifteen unused positions at the end of a frame.

use bytemuck::{Pod, Zeroable};
use std::fmt;
use std::mem::size_of;
use std::ops::Range;
use std::sync::Arc;

/// Positions per lane. Sixteen triples is 192 bytes — three cache lines exactly,
/// which is what lets the lane demand 64-byte alignment without padding.
const LANE: usize = 16;

#[repr(C, align(64))]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CoordLane([[f32; 3]; LANE]);

/// How many times the coordinates of a structure have changed.
///
/// Everything derived from positions — spatial indices, cached distances, views
/// handed to a caller — records the generation it was built against. Using one
/// after the coordinates moved is then a detectable mistake rather than a
/// silently wrong answer.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
#[repr(transparent)]
pub struct CoordinateGeneration(u64);

impl CoordinateGeneration {
    /// The generation of a structure that has not been edited.
    pub const INITIAL: Self = Self(0);

    /// Rebuilds a generation from a value previously produced by this type.
    ///
    /// The counter is opaque to callers that only need freshness checks, but
    /// preserving it at an interop boundary is useful when a cached view is
    /// serialised or handed between language runtimes.  No validation is
    /// required because every `u64` is a representable generation.
    #[must_use]
    pub const fn from_raw(value: u64) -> Self {
        Self(value)
    }

    /// Returns the next generation, or `None` when the counter is exhausted.
    #[must_use]
    pub const fn next(self) -> Option<Self> {
        match self.0.checked_add(1) {
            Some(next) => Some(Self(next)),
            None => None,
        }
    }

    /// The raw counter.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

/// An axis-aligned bounding box.
///
/// Every chunk carries one. A query that wants atoms within some distance of a
/// region can reject a whole chunk by comparing boxes, before reading a single
/// position out of it.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Aabb {
    /// Componentwise minimum.
    pub min: [f32; 3],
    /// Componentwise maximum.
    pub max: [f32; 3],
}

impl Aabb {
    /// An empty box, which absorbs any point it is extended by.
    pub const EMPTY: Self = Self {
        min: [f32::INFINITY; 3],
        max: [f32::NEG_INFINITY; 3],
    };

    /// Returns true when no point has been added.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.min[0] > self.max[0] || self.min[1] > self.max[1] || self.min[2] > self.max[2]
    }

    /// Grows the box to contain `point`.
    ///
    /// Non-finite components are ignored: a coordinate that is not a number
    /// says nothing about where the atoms are, and letting it poison the box
    /// would disable chunk skipping for the whole structure.
    #[inline]
    pub fn extend(&mut self, point: [f32; 3]) {
        extend_axis(&mut self.min[0], &mut self.max[0], point[0]);
        extend_axis(&mut self.min[1], &mut self.max[1], point[1]);
        extend_axis(&mut self.min[2], &mut self.max[2], point[2]);
    }

    /// Grows the box to contain another box.
    #[inline]
    pub fn union(&mut self, other: &Self) {
        if other.is_empty() {
            return;
        }

        self.extend(other.min);
        self.extend(other.max);
    }

    /// Returns true when the two boxes come within `distance` of one another.
    ///
    /// Compares squared separations, so no square root is taken.
    #[must_use]
    #[inline]
    pub fn within(&self, other: &Self, distance: f32) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }

        let maximum_gap_squared = distance * distance;

        let x = axis_gap(self.min[0], self.max[0], other.min[0], other.max[0]);
        let mut gap_squared = x * x;

        if gap_squared > maximum_gap_squared {
            return false;
        }

        let y = axis_gap(self.min[1], self.max[1], other.min[1], other.max[1]);
        gap_squared += y * y;

        if gap_squared > maximum_gap_squared {
            return false;
        }

        let z = axis_gap(self.min[2], self.max[2], other.min[2], other.max[2]);
        gap_squared += z * z;

        gap_squared <= maximum_gap_squared
    }
}

impl Default for Aabb {
    /// Creates an empty bounding box.
    fn default() -> Self {
        Self::EMPTY
    }
}

/// The positions of one frame.
///
/// # Examples
///
/// ```
/// use pdbiox_core::CoordinateBlock;
///
/// let mut block = CoordinateBlock::with_capacity(3);
/// block.push([1.0, 2.0, 3.0]);
/// block.push([4.0, 5.0, 6.0]);
///
/// assert_eq!(block.len(), 2);
/// assert_eq!(block.as_slice()[1], [4.0, 5.0, 6.0]);
/// ```
#[derive(Clone, Default)]
pub struct CoordinateBlock {
    lanes: Arc<Vec<CoordLane>>,
    len: u32,
}

impl CoordinateBlock {
    /// Creates an empty block.
    #[must_use]
    pub fn new() -> Self {
        Self {
            lanes: Arc::new(Vec::new()),
            len: 0,
        }
    }

    /// Creates a block with room for `positions` without reallocating.
    #[must_use]
    pub fn with_capacity(positions: usize) -> Self {
        Self {
            lanes: Arc::new(Vec::with_capacity(lanes_for_positions(positions))),
            len: 0,
        }
    }

    /// The number of positions held.
    #[must_use]
    pub const fn len(&self) -> u32 {
        self.len
    }

    /// Returns true when the block holds no positions.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Appends a position.
    ///
    /// Copy-on-write detachment is performed at most once for the append, and
    /// the visible length changes only after the target slot has been written.
    pub fn push(&mut self, position: [f32; 3]) {
        let Some((index, next_len)) = append_indices(self.len) else {
            return;
        };

        let lanes = Arc::make_mut(&mut self.lanes);

        if write_position(lanes, index, position) {
            self.len = next_len;
        }
    }

    /// All positions, contiguously.
    ///
    /// This slice is the whole point of the layout: it is what an array
    /// interface receives without a copy.
    #[must_use]
    #[inline]
    pub fn as_slice(&self) -> &[[f32; 3]] {
        let all: &[[f32; 3]] = bytemuck::cast_slice(self.lanes.as_slice());
        visible_slice(all, self.len)
    }

    /// All positions, mutably.
    ///
    /// Reaching for this outside a scoped edit skips the generation bump that
    /// invalidates derived state, which is why it is not part of the public
    /// structure API.
    #[inline]
    pub(crate) fn as_mut_slice(&mut self) -> &mut [[f32; 3]] {
        let lanes = Arc::make_mut(&mut self.lanes);
        let all: &mut [[f32; 3]] = bytemuck::cast_slice_mut(lanes.as_mut_slice());

        visible_mut_slice(all, self.len)
    }

    /// The positions a chunk covers, or `None` if the range runs past the end.
    #[must_use]
    #[inline]
    pub fn range(&self, range: Range<u32>) -> Option<&[[f32; 3]]> {
        let range = usize_range(range)?;
        self.as_slice().get(range)
    }

    /// The bounding box of a range of positions.
    #[must_use]
    pub fn bounds(&self, range: Range<u32>) -> Aabb {
        let Some(positions) = self.range(range) else {
            return Aabb::EMPTY;
        };

        let mut bounds = Aabb::EMPTY;

        for &position in positions {
            bounds.extend(position);
        }

        bounds
    }

    /// Bytes the block occupies, including the unused tail of the last lane.
    #[must_use]
    pub fn allocated_bytes(&self) -> usize {
        self.lanes.len() * size_of::<CoordLane>()
    }
}

impl fmt::Debug for CoordinateBlock {
    /// Formats block metadata without dumping the potentially large coordinate buffer.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CoordinateBlock")
            .field("positions", &self.len)
            .finish_non_exhaustive()
    }
}

impl FromIterator<[f32; 3]> for CoordinateBlock {
    /// Collects positions directly into aligned coordinate lanes.
    ///
    /// Construction writes directly into the backing `Vec`, avoiding an
    /// `Arc::make_mut` uniqueness check for every collected position.
    fn from_iter<T: IntoIterator<Item = [f32; 3]>>(iter: T) -> Self {
        let mut iter = iter.into_iter();
        let hinted_positions = iter.size_hint().0;

        let mut lanes = Vec::new();
        let hinted_lanes = lanes_for_positions(hinted_positions);

        // A size hint is advisory and may come from arbitrary iterator code.
        // Failure to honour it must not prevent incremental collection.
        if hinted_lanes != 0 {
            let _reserve_failed = lanes.try_reserve_exact(hinted_lanes).is_err();
        }

        let mut len = 0u32;

        for position in &mut iter {
            let Some((index, next_len)) = append_indices(len) else {
                break;
            };

            if !write_position(&mut lanes, index, position) {
                break;
            }

            len = next_len;
        }

        Self {
            lanes: Arc::new(lanes),
            len,
        }
    }
}

/// Extends one bounding-box axis when `value` is finite.
#[inline]
fn extend_axis(minimum: &mut f32, maximum: &mut f32, value: f32) {
    if !value.is_finite() {
        return;
    }

    *minimum = minimum.min(value);
    *maximum = maximum.max(value);
}

/// Returns the non-overlapping gap between two intervals on one axis.
#[inline]
fn axis_gap(own_min: f32, own_max: f32, other_min: f32, other_max: f32) -> f32 {
    (other_min - own_max).max(own_min - other_max).max(0.0)
}

/// Converts a position count into the number of aligned lanes required.
#[inline]
fn lanes_for_positions(positions: usize) -> usize {
    positions.div_ceil(LANE)
}

/// Computes the backing-buffer index and next visible length for an append.
///
/// Returns `None` when the block has exhausted its `u32` position space or the
/// current index cannot be represented by the platform.
#[inline]
fn append_indices(len: u32) -> Option<(usize, u32)> {
    let next_len = len.checked_add(1)?;
    let index = usize::try_from(len).ok()?;

    Some((index, next_len))
}

/// Writes one position into its lane, growing the lane vector when necessary.
///
/// Returns `false` only when the supplied logical index is inconsistent with
/// the current contiguous lane layout.
#[inline]
fn write_position(lanes: &mut Vec<CoordLane>, index: usize, position: [f32; 3]) -> bool {
    let (lane_index, slot) = lane_position(index);

    if lane_index > lanes.len() {
        return false;
    }

    if lane_index == lanes.len() {
        let mut lane = CoordLane([[0.0; 3]; LANE]);

        let Some(target) = lane.0.get_mut(slot) else {
            return false;
        };

        *target = position;
        lanes.push(lane);
        return true;
    }

    let Some(target) = lanes
        .get_mut(lane_index)
        .and_then(|lane| lane.0.get_mut(slot))
    else {
        return false;
    };

    *target = position;
    true
}

/// Maps a logical position to its lane and in-lane slot.
#[inline]
const fn lane_position(position: usize) -> (usize, usize) {
    (position / LANE, position % LANE)
}

/// Restricts a backing slice to the number of logically visible values.
///
/// An impossible length is treated as an empty view rather than causing an
/// out-of-bounds panic.
#[inline]
fn visible_slice<T>(values: &[T], len: u32) -> &[T] {
    let Ok(len) = usize::try_from(len) else {
        return &[];
    };

    match values.get(..len) {
        Some(values) => values,
        None => &[],
    }
}

/// Restricts a mutable backing slice to the number of logically visible values.
///
/// An impossible length is treated as an empty view rather than exposing
/// storage outside the block's logical extent.
#[inline]
fn visible_mut_slice<T>(values: &mut [T], len: u32) -> &mut [T] {
    let Ok(len) = usize::try_from(len) else {
        return &mut [];
    };

    match values.get_mut(..len) {
        Some(values) => values,
        None => &mut [],
    }
}

/// Converts an atom-position range into platform-native slice indices.
#[inline]
fn usize_range(range: Range<u32>) -> Option<Range<usize>> {
    let start = usize::try_from(range.start).ok()?;
    let end = usize::try_from(range.end).ok()?;

    Some(start..end)
}

#[cfg(test)]
#[path = "storage_tests.rs"]
mod tests;
