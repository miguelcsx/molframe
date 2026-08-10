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
    pub fn is_empty(&self) -> bool {
        self.min.iter().zip(&self.max).any(|(min, max)| min > max)
    }

    /// Grows the box to contain `point`.
    ///
    /// Non-finite components are ignored: a coordinate that is not a number
    /// says nothing about where the atoms are, and letting it poison the box
    /// would disable chunk skipping for the whole structure.
    pub fn extend(&mut self, point: [f32; 3]) {
        for ((min, max), value) in self.min.iter_mut().zip(&mut self.max).zip(point) {
            if !value.is_finite() {
                continue;
            }

            *min = min.min(value);
            *max = max.max(value);
        }
    }

    /// Grows the box to contain another box.
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
    pub fn within(&self, other: &Self, distance: f32) -> bool {
        if self.is_empty() || other.is_empty() {
            return false;
        }

        let maximum_gap_squared = distance * distance;
        let mut gap_squared = 0.0f32;

        for axis in 0..3 {
            let gap = axis_gap(
                self.min[axis],
                self.max[axis],
                other.min[axis],
                other.max[axis],
            );

            gap_squared += gap * gap;

            if gap_squared > maximum_gap_squared {
                return false;
            }
        }

        gap_squared <= maximum_gap_squared
    }
}

impl Default for Aabb {
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
            lanes: Arc::new(Vec::with_capacity(positions.div_ceil(LANE))),
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
    pub fn push(&mut self, position: [f32; 3]) {
        let index = self.len as usize;

        if index == self.lanes.len() * LANE {
            Arc::make_mut(&mut self.lanes).push(CoordLane([[0.0; 3]; LANE]));
        }

        let (lane, slot) = lane_position(index);

        let Some(target) = Arc::make_mut(&mut self.lanes)
            .get_mut(lane)
            .and_then(|lane| lane.0.get_mut(slot))
        else {
            return;
        };

        *target = position;
        self.len += 1;
    }

    /// All positions, contiguously.
    ///
    /// This slice is the whole point of the layout: it is what an array
    /// interface receives without a copy.
    #[must_use]
    pub fn as_slice(&self) -> &[[f32; 3]] {
        let all: &[[f32; 3]] = bytemuck::cast_slice(self.lanes.as_slice());

        visible_slice(all, self.len)
    }

    /// All positions, mutably.
    ///
    /// Reaching for this outside a scoped edit skips the generation bump that
    /// invalidates derived state, which is why it is not part of the public
    /// structure API.
    pub(crate) fn as_mut_slice(&mut self) -> &mut [[f32; 3]] {
        let all: &mut [[f32; 3]] =
            bytemuck::cast_slice_mut(Arc::make_mut(&mut self.lanes).as_mut_slice());

        visible_mut_slice(all, self.len)
    }

    /// The positions a chunk covers, or `None` if the range runs past the end.
    #[must_use]
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

        positions
            .iter()
            .copied()
            .fold(Aabb::EMPTY, |mut bounds, position| {
                bounds.extend(position);
                bounds
            })
    }

    /// Bytes the block occupies, including the unused tail of the last lane.
    #[must_use]
    pub fn allocated_bytes(&self) -> usize {
        self.lanes.len() * size_of::<CoordLane>()
    }
}

impl fmt::Debug for CoordinateBlock {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CoordinateBlock")
            .field("positions", &self.len)
            .finish_non_exhaustive()
    }
}

impl FromIterator<[f32; 3]> for CoordinateBlock {
    fn from_iter<T: IntoIterator<Item = [f32; 3]>>(iter: T) -> Self {
        let iter = iter.into_iter();
        let mut block = Self::with_capacity(iter.size_hint().0);

        for position in iter {
            block.push(position);
        }

        block
    }
}

fn axis_gap(own_min: f32, own_max: f32, other_min: f32, other_max: f32) -> f32 {
    (other_min - own_max).max(own_min - other_max).max(0.0)
}

const fn lane_position(position: usize) -> (usize, usize) {
    (position / LANE, position % LANE)
}

fn visible_slice<T>(values: &[T], len: u32) -> &[T] {
    let Ok(len) = usize::try_from(len) else {
        return &[];
    };

    match values.get(..len) {
        Some(values) => values,
        None => &[],
    }
}

fn visible_mut_slice<T>(values: &mut [T], len: u32) -> &mut [T] {
    let Ok(len) = usize::try_from(len) else {
        return &mut [];
    };

    match values.get_mut(..len) {
        Some(values) => values,
        None => &mut [],
    }
}

fn usize_range(range: Range<u32>) -> Option<Range<usize>> {
    let start = usize::try_from(range.start).ok()?;
    let end = usize::try_from(range.end).ok()?;

    Some(start..end)
}

#[cfg(test)]
#[path = "storage_tests.rs"]
mod tests;
