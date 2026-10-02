//! A uniform cell grid over points, shared by the neighbour searches.
//!
//! Points are bucketed by the cell that holds them and stored in cell order, so
//! a neighbouring cell is one contiguous run. The cells of the bounding box are
//! tabulated when there are few enough of them; a span too large for that is
//! searched through its occupied cells instead.

use num_traits::ToPrimitive;
use std::ops::Range;

/// The most cells a dense grid may hold.
const DENSE_CELL_LIMIT: usize = 1 << 22;

/// The integer cell a point falls in.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct CellKey(pub(crate) [i32; 3]);

/// The cell of `position` for cells `size` ångström wide.
///
/// `None` for a position that is not finite, a size that is not positive and
/// finite, or a cell index outside `i32`.
pub(crate) fn cell_for(position: [f32; 3], size: f32) -> Option<CellKey> {
    if position.iter().all(|value| value.is_finite()) && size.is_finite() && size > 0.0 {
        Some(CellKey([
            (position[0] / size).floor().to_i32()?,
            (position[1] / size).floor().to_i32()?,
            (position[2] / size).floor().to_i32()?,
        ]))
    } else {
        None
    }
}

/// Items ordered by grid cell; within a cell they keep their input order.
pub(crate) struct CellGrid<T> {
    items: Vec<T>,
    layout: Layout,
}

enum Layout {
    /// One slot per cell of the bounding box: a neighbour is an index away.
    Dense { offsets: Vec<u32>, dims: [usize; 3] },
    /// Only occupied cells, searched by key.
    Sparse { cells: Vec<(CellKey, Range<usize>)> },
}

impl<T: Copy> CellGrid<T> {
    /// Buckets `entries` into a dense grid when the bounding box is small
    /// enough, and a sparse one otherwise.
    pub(crate) fn build(entries: Vec<(CellKey, T)>) -> Self {
        match Self::dense(&entries) {
            Some(grid) => grid,
            None => Self::sparse(entries),
        }
    }

    /// Counting-sorts the entries into the cells of their bounding box.
    ///
    /// Returns `None` when there are no entries or the box would hold more
    /// cells than the limit.
    pub(crate) fn dense(entries: &[(CellKey, T)]) -> Option<Self> {
        let seed = entries.first()?;
        let (mut low, mut high) = (seed.0.0, seed.0.0);
        for (cell, _) in entries {
            for axis in 0..3 {
                low[axis] = low[axis].min(cell.0[axis]);
                high[axis] = high[axis].max(cell.0[axis]);
            }
        }
        let mut dims = [0_usize; 3];
        let mut total = 1_usize;
        for axis in 0..3 {
            let span = i64::from(high[axis]) - i64::from(low[axis]) + 1;
            dims[axis] = usize::try_from(span).ok()?;
            total = total.checked_mul(dims[axis])?;
        }
        if total > DENSE_CELL_LIMIT {
            return None;
        }
        let slot = |cell: &CellKey| -> usize {
            // Every entry lies inside the bounding box, so the offsets from its
            // corner are never negative; a failed conversion cannot occur.
            let at =
                |axis: usize| match usize::try_from(i64::from(cell.0[axis]) - i64::from(low[axis]))
                {
                    Ok(offset) => offset,
                    Err(_) => 0,
                };
            (at(0) * dims[1] + at(1)) * dims[2] + at(2)
        };
        let mut offsets = vec![0_u32; total + 1];
        for (cell, _) in entries {
            offsets[slot(cell) + 1] += 1;
        }
        for index in 1..offsets.len() {
            offsets[index] += offsets[index - 1];
        }
        // Filling in entry order keeps each cell's items in input order.
        let mut fill: Vec<u32> = offsets[..total].to_vec();
        let mut items = vec![seed.1; entries.len()];
        for (cell, item) in entries {
            let target = &mut fill[slot(cell)];
            items[*target as usize] = *item;
            *target += 1;
        }
        Some(Self {
            items,
            layout: Layout::Dense { offsets, dims },
        })
    }

    /// Sorts the entries by cell, stably, and records the occupied cells.
    pub(crate) fn sparse(mut entries: Vec<(CellKey, T)>) -> Self {
        entries.sort_by_key(|(cell, _)| *cell);
        let mut cells: Vec<(CellKey, Range<usize>)> = Vec::new();
        for (position, (cell, _)) in entries.iter().enumerate() {
            match cells.last_mut() {
                Some((last, range)) if last == cell => range.end = position + 1,
                _ => cells.push((*cell, position..position + 1)),
            }
        }
        Self {
            items: entries.into_iter().map(|(_, item)| item).collect(),
            layout: Layout::Sparse { cells },
        }
    }

    /// The items in cell order.
    pub(crate) fn items(&self) -> &[T] {
        &self.items
    }

    /// Whether the cells of the bounding box are tabulated.
    #[cfg(test)]
    pub(crate) fn is_dense(&self) -> bool {
        matches!(self.layout, Layout::Dense { .. })
    }

    /// The number of cells [`Self::visit_cells`] can be asked to visit.
    ///
    /// A dense grid counts every slot of its bounding box, empty ones included;
    /// a sparse grid counts its occupied cells.
    pub(crate) fn cell_count(&self) -> usize {
        match &self.layout {
            Layout::Dense { dims, .. } => dims[0] * dims[1] * dims[2],
            Layout::Sparse { cells } => cells.len(),
        }
    }

    /// Calls `visit` with each occupied cell's items and the items of every
    /// occupied cell around it, itself included, as ranges into [`Self::items`].
    pub(crate) fn for_each_cell(&self, visit: impl FnMut(Range<usize>, &[Range<usize>])) {
        self.visit_cells(0..self.cell_count(), visit);
    }

    /// [`Self::for_each_cell`] restricted to cells `cells` of `0..cell_count()`.
    ///
    /// Cells are numbered in the order `for_each_cell` visits them, so visiting
    /// consecutive ranges in turn reproduces its whole sequence.
    pub(crate) fn visit_cells(
        &self,
        cells: Range<usize>,
        mut visit: impl FnMut(Range<usize>, &[Range<usize>]),
    ) {
        let mut neighbourhood = Vec::with_capacity(27);
        match &self.layout {
            Layout::Dense { offsets, dims } => {
                let range_at = |x: usize, y: usize, z: usize| {
                    let slot = (x * dims[1] + y) * dims[2] + z;
                    offsets[slot] as usize..offsets[slot + 1] as usize
                };
                for cell in cells.start..cells.end.min(self.cell_count()) {
                    if offsets[cell] == offsets[cell + 1] {
                        continue;
                    }
                    let (x, y, z) = (
                        cell / (dims[1] * dims[2]),
                        cell / dims[2] % dims[1],
                        cell % dims[2],
                    );
                    neighbourhood.clear();
                    for nx in x.saturating_sub(1)..=(x + 1).min(dims[0] - 1) {
                        for ny in y.saturating_sub(1)..=(y + 1).min(dims[1] - 1) {
                            for nz in z.saturating_sub(1)..=(z + 1).min(dims[2] - 1) {
                                let range = range_at(nx, ny, nz);
                                if !range.is_empty() {
                                    neighbourhood.push(range);
                                }
                            }
                        }
                    }
                    visit(range_at(x, y, z), &neighbourhood);
                }
            }
            Layout::Sparse { cells: occupied } => {
                let end = cells.end.min(occupied.len());
                let Some(visited) = occupied.get(cells.start.min(end)..end) else {
                    return;
                };
                for (base, own) in visited {
                    neighbourhood.clear();
                    for dx in -1..=1 {
                        for dy in -1..=1 {
                            for dz in -1..=1 {
                                let key = CellKey([
                                    base.0[0].saturating_add(dx),
                                    base.0[1].saturating_add(dy),
                                    base.0[2].saturating_add(dz),
                                ]);
                                if let Ok(index) =
                                    occupied.binary_search_by_key(&key, |(key, _)| *key)
                                {
                                    neighbourhood.push(occupied[index].1.clone());
                                }
                            }
                        }
                    }
                    visit(own.clone(), &neighbourhood);
                }
            }
        }
    }
}
