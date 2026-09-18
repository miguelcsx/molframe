//! Compressed sparse row storage for per-item adjacency.
//!
//! A `Vec<Vec<T>>` adjacency costs one heap allocation and one 24-byte header
//! per row, so a billion-atom neighbour list spends 24 GB and a billion
//! allocator round trips before storing a single edge. `Csr` holds the same
//! information as two flat vectors: a row's items are a contiguous slice, found
//! by two loads and a subtraction.
//!
//! Offsets are `usize` rather than `u32` because edge counts outrun row counts:
//! a billion atoms at fifty neighbours each is 5e10 edges, well past `u32`.

/// Rows of `T` stored contiguously, addressed by a prefix-sum offset table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Csr<T> {
    /// Row `i` occupies `values[offsets[i]..offsets[i + 1]]`.
    ///
    /// Always holds `rows + 1` entries, so the final row needs no special case.
    offsets: Vec<usize>,
    /// Every row's items, concatenated in row order.
    values: Vec<T>,
}

impl<T> Csr<T> {
    /// An adjacency with no rows at all.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            offsets: vec![0],
            values: Vec::new(),
        }
    }

    /// The number of rows.
    #[must_use]
    pub fn rows(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Whether there are no rows.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.rows() == 0
    }

    /// The total number of stored items across every row.
    #[must_use]
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// The items of row `index`, or an empty slice when it is out of range.
    ///
    /// Two loads and a subtraction, with no pointer chase.
    #[must_use]
    pub fn row(&self, index: usize) -> &[T] {
        let (Some(start), Some(end)) = (
            self.offsets.get(index),
            index.checked_add(1).and_then(|next| self.offsets.get(next)),
        ) else {
            return &[];
        };
        match self.values.get(*start..*end) {
            Some(items) => items,
            None => &[],
        }
    }

    /// The items of row `index` as a mutable slice, empty when out of range.
    pub fn row_mut(&mut self, index: usize) -> &mut [T] {
        let (Some(start), Some(end)) = (
            self.offsets.get(index),
            index.checked_add(1).and_then(|next| self.offsets.get(next)),
        ) else {
            return &mut [];
        };
        let (start, end) = (*start, *end);
        match self.values.get_mut(start..end) {
            Some(items) => items,
            None => &mut [],
        }
    }

    /// Every row in order, as slices.
    pub fn iter(&self) -> impl Iterator<Item = &[T]> {
        (0..self.rows()).map(|index| self.row(index))
    }
}

/// Scatters items into their rows before freezing into a [`Csr`].
///
/// The shape is fixed at construction, so scattering is a write to a cursor and
/// never reallocates or moves an earlier row.
#[derive(Debug)]
pub struct CsrBuilder<T> {
    /// The next free slot of each row, starting at that row's offset.
    cursors: Vec<usize>,
    /// The finished prefix sums, moved into the `Csr` unchanged.
    offsets: Vec<usize>,
    /// Value storage, sized to the declared total before any row is read.
    values: Vec<T>,
}

impl<T: Clone + Default> CsrBuilder<T> {
    /// Allocates the offset table for a known shape and fills value storage.
    ///
    /// `degrees[i]` is the exact item count of row `i`. Callers that already
    /// run a counting pass — the usual shape for a neighbour search — pay no
    /// reallocation at all: this reserves the exact total, and rows may then be
    /// scattered into in any order.
    #[must_use]
    pub fn with_degrees(degrees: &[usize]) -> Self {
        let mut offsets = Vec::with_capacity(degrees.len() + 1);
        let mut total = 0usize;
        offsets.push(0);
        for degree in degrees {
            total = total.saturating_add(*degree);
            offsets.push(total);
        }

        Self {
            cursors: offsets[..degrees.len()].to_vec(),
            offsets,
            values: vec![T::default(); total],
        }
    }

    /// Writes `item` into the next free slot of `row`.
    ///
    /// Items beyond a row's declared degree, and rows beyond the declared row
    /// count, are dropped rather than corrupting a neighbouring row.
    pub fn push(&mut self, row: usize, item: T) {
        let (Some(cursor), Some(end)) = (self.cursors.get_mut(row), self.offsets.get(row + 1))
        else {
            return;
        };
        if *cursor >= *end {
            return;
        }
        if let Some(slot) = self.values.get_mut(*cursor) {
            *slot = item;
            *cursor += 1;
        }
    }

    /// Freezes the scattered rows, closing any gap left by an under-filled row.
    ///
    /// Rows given a degree larger than the number of items actually pushed keep
    /// no default padding: the surviving items are shifted down and the offsets
    /// rebuilt, so every row is exactly as long as what reached it.
    #[must_use]
    pub fn finish(mut self) -> Csr<T> {
        let mut offsets = Vec::with_capacity(self.offsets.len());
        let mut write = 0usize;
        offsets.push(0);

        for row in 0..self.cursors.len() {
            let (Some(start), Some(cursor)) = (self.offsets.get(row), self.cursors.get(row)) else {
                break;
            };
            for read in *start..*cursor {
                self.values.swap(write, read);
                write += 1;
            }
            offsets.push(write);
        }

        self.values.truncate(write);
        Csr {
            offsets,
            values: self.values,
        }
    }
}

impl<T: Clone + Default + Ord> CsrBuilder<T> {
    /// Freezes the scattered rows with each row sorted and deduplicated.
    ///
    /// A candidate-pair search may report the same pair from two cells, so an
    /// adjacency built from one is canonicalised before use. Sorting happens
    /// within each row's own slice and compaction is a single linear pass, so
    /// no per-row allocation is involved.
    #[must_use]
    pub fn finish_canonical(mut self) -> Csr<T> {
        let mut offsets = Vec::with_capacity(self.offsets.len());
        let mut write = 0usize;
        offsets.push(0);

        for row in 0..self.cursors.len() {
            let (Some(start), Some(cursor)) = (self.offsets.get(row), self.cursors.get(row)) else {
                break;
            };
            let (start, cursor) = (*start, *cursor);
            if let Some(items) = self.values.get_mut(start..cursor) {
                items.sort_unstable();
            }

            let row_start = write;
            for read in start..cursor {
                if write > row_start && self.values.get(write - 1) == self.values.get(read) {
                    continue;
                }
                self.values.swap(write, read);
                write += 1;
            }
            offsets.push(write);
        }

        self.values.truncate(write);
        Csr {
            offsets,
            values: self.values,
        }
    }
}

#[cfg(test)]
#[path = "csr_tests.rs"]
mod tests;
