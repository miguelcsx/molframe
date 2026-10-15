//! Flat cell index over candidate site positions.
//!
//! Entries are sorted by cell once and stored in a single array, so a cell is
//! a slice rather than an allocation. A bounding box small enough to
//! tabulate gets a dense offset table (one array read per cell); a sparse or
//! far-flung set falls back to binary search over the occupied cells, so the
//! memory spent never exceeds a few words per entry.

pub(crate) type Point = [f32; 3];
pub(crate) type Cell = [i64; 3];

/// Dense tables are used while they cost at most this many cells per entry.
const DENSE_CELLS_PER_ENTRY: u128 = 4;
/// Dense tables are always allowed up to this many cells.
const DENSE_CELL_FLOOR: u128 = 4096;

#[derive(Debug)]
enum Index {
    /// `offsets[id]..offsets[id + 1]` for the row-major id inside the box.
    Dense {
        low: Cell,
        dims: [usize; 3],
        offsets: Vec<u32>,
    },
    /// Occupied cells in lexicographic order with their start offsets.
    Sparse { keys: Vec<Cell>, offsets: Vec<u32> },
}

#[derive(Debug)]
pub(crate) struct Grid {
    cell: f64,
    entries: Vec<(u32, Point)>,
    index: Index,
}

impl Grid {
    /// Indexes `items` (candidate ordinal and position) into cells of side `cell`.
    pub(crate) fn build(cell: f64, items: impl Iterator<Item = (u32, Point)>) -> Self {
        let mut grid = Self {
            cell,
            entries: Vec::new(),
            index: Index::Sparse {
                keys: Vec::new(),
                offsets: vec![0],
            },
        };
        let mut keyed: Vec<(Cell, u32, Point)> = items
            .map(|(index, point)| (grid.key(point), index, point))
            .collect();
        keyed.sort_unstable_by_key(|&(key, index, _)| (key, index));
        grid.index = index_of(&keyed);
        grid.entries = keyed
            .into_iter()
            .map(|(_, index, point)| (index, point))
            .collect();
        grid
    }

    pub(crate) fn key(&self, point: Point) -> Cell {
        let scaled = |value: f32| {
            // Saturating: a coordinate this far out cannot be near anything.
            molframe_core::numeric::f64_to_i64((f64::from(value) / self.cell).floor())
        };
        [scaled(point[0]), scaled(point[1]), scaled(point[2])]
    }

    /// Calls `visit` for every entry in the 27 cells around `point`'s cell.
    pub(crate) fn for_each_near(&self, point: Point, mut visit: impl FnMut(u32, Point)) {
        let centre = self.key(point);
        for dx in -1..=1_i64 {
            for dy in -1..=1_i64 {
                for dz in -1..=1_i64 {
                    let (Some(x), Some(y), Some(z)) = (
                        centre[0].checked_add(dx),
                        centre[1].checked_add(dy),
                        centre[2].checked_add(dz),
                    ) else {
                        continue;
                    };
                    let key = [x, y, z];
                    for &(index, other) in self.bucket(key) {
                        visit(index, other);
                    }
                }
            }
        }
    }

    fn bucket(&self, key: Cell) -> &[(u32, Point)] {
        let span = match &self.index {
            Index::Dense { low, dims, offsets } => dense_id(key, *low, *dims)
                .and_then(|id| Some((*offsets.get(id)?, *offsets.get(id + 1)?))),
            Index::Sparse { keys, offsets } => keys
                .binary_search(&key)
                .ok()
                .and_then(|id| Some((*offsets.get(id)?, *offsets.get(id + 1)?))),
        };
        let Some((start, end)) = span else {
            return &[];
        };
        match self.entries.get(start as usize..end as usize) {
            Some(slice) => slice,
            None => &[],
        }
    }
}

fn dense_id(key: Cell, low: Cell, dims: [usize; 3]) -> Option<usize> {
    let mut id = 0usize;
    for axis in 0..3 {
        let offset = usize::try_from(i128::from(key[axis]) - i128::from(low[axis])).ok()?;
        if offset >= dims[axis] {
            return None;
        }
        id = id * dims[axis] + offset;
    }
    Some(id)
}

/// The offset table for entries already sorted by cell.
fn index_of(sorted: &[(Cell, u32, Point)]) -> Index {
    let (Some(first), Some(last)) = (sorted.first(), sorted.last()) else {
        return Index::Sparse {
            keys: Vec::new(),
            offsets: vec![0],
        };
    };
    let mut low = first.0;
    let mut high = first.0;
    for (key, ..) in sorted {
        for axis in 0..3 {
            low[axis] = low[axis].min(key[axis]);
            high[axis] = high[axis].max(key[axis]);
        }
    }
    debug_assert!(first.0 <= last.0);
    let extent = |axis: usize| (i128::from(high[axis]) - i128::from(low[axis]) + 1).unsigned_abs();
    let cells = extent(0)
        .saturating_mul(extent(1))
        .saturating_mul(extent(2));
    let budget = (DENSE_CELLS_PER_ENTRY * sorted.len() as u128).max(DENSE_CELL_FLOOR);
    let dims = [extent(0), extent(1), extent(2)].map(|value| match usize::try_from(value) {
        Ok(dim) => dim,
        Err(_) => 0,
    });
    if cells <= budget && dims.iter().all(|&dim| dim > 0) {
        let mut offsets = vec![0u32; dims[0] * dims[1] * dims[2] + 1];
        for (key, ..) in sorted {
            if let Some(id) = dense_id(*key, low, dims) {
                offsets[id + 1] += 1;
            }
        }
        for id in 1..offsets.len() {
            offsets[id] += offsets[id - 1];
        }
        return Index::Dense { low, dims, offsets };
    }
    let mut keys: Vec<Cell> = Vec::new();
    let mut offsets: Vec<u32> = Vec::new();
    for (position, (key, ..)) in sorted.iter().enumerate() {
        if keys.last() != Some(key) {
            keys.push(*key);
            offsets.push(molframe_core::numeric::usize_to_u32(position));
        }
    }
    offsets.push(molframe_core::numeric::usize_to_u32(sorted.len()));
    Index::Sparse { keys, offsets }
}

#[cfg(test)]
#[path = "cells_tests.rs"]
mod tests;
