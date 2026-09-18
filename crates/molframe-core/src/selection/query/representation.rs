//! Selection shapes, construction, membership and iteration.

use crate::column::{BitVec, Ones};
use smallvec::SmallVec;
use std::ops::Range;

/// A set of atoms, addressed by position in the structure's flat atom order.
///
/// # Examples
///
/// ```
/// use molframe_core::AtomSelection;
///
/// let chain = AtomSelection::range(0..100);
/// let metals = AtomSelection::from_sorted(vec![5, 50, 500]);
/// assert_eq!(chain.intersect(&metals).len(), 2);
/// ```
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub enum AtomSelection {
    /// No atoms.
    #[default]
    Empty,
    /// Every atom of a structure with this many.
    All(u32),
    /// One contiguous run.
    Range(Range<u32>),
    /// Several contiguous runs, ascending and disjoint.
    Ranges(SmallVec<[Range<u32>; 4]>),
    /// Scattered positions, ascending and unique.
    Sparse(Vec<u32>),
    /// A bit per position, for a set dense enough to make that cheaper.
    Dense(BitVec),
}

const DENSE_FRACTION: u32 = 16;

impl AtomSelection {
    /// A selection covering one contiguous run.
    #[must_use]
    pub fn range(atoms: Range<u32>) -> Self {
        if atoms.start >= atoms.end {
            Self::Empty
        } else {
            Self::Range(atoms)
        }
    }

    /// Chooses the cheapest shape for ascending unique positions.
    ///
    /// # Panics
    /// Panics if `positions` is not ascending and unique, or if its storage
    /// exceeds the platform's addressable range. Collect arbitrary positions
    /// into `AtomSelection` to sort and deduplicate them first.
    #[must_use]
    pub fn from_sorted(positions: Vec<u32>) -> Self {
        let (Some(first), Some(last)) = (positions.first().copied(), positions.last().copied())
        else {
            return Self::Empty;
        };
        let span = u64::from(last - first) + 1;
        if last < u32::MAX && span == positions.len() as u64 {
            return Self::Range(first..last + 1);
        }
        if let Some(runs) = inline_runs(&positions) {
            let covered: u64 = runs.iter().map(|run| u64::from(run.end - run.start)).sum();
            if covered == positions.len() as u64 && runs.len() * 2 < positions.len() {
                return Self::Ranges(runs);
            }
        }
        if let Some(mask_len) = last.checked_add(1)
            && positions.len() as u64 * u64::from(DENSE_FRACTION) >= u64::from(mask_len)
        {
            return Self::Dense(mask_from_positions(mask_len, &positions));
        }
        Self::Sparse(positions)
    }

    /// The number of atoms selected.
    ///
    /// # Panics
    /// Panics if manually constructed storage violates normalized `u32` ranges.
    #[must_use]
    pub fn len(&self) -> u64 {
        selection_len(self)
    }

    /// Returns true when nothing is selected.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns true when `atom` is selected.
    #[must_use]
    pub fn contains(&self, atom: u32) -> bool {
        match self {
            Self::Empty => false,
            Self::All(count) => atom < *count,
            Self::Range(run) => run.contains(&atom),
            Self::Ranges(runs) => {
                let candidate = runs.partition_point(|run| run.end <= atom);
                runs.get(candidate).is_some_and(|run| run.contains(&atom))
            }
            Self::Sparse(positions) => positions.binary_search(&atom).is_ok(),
            Self::Dense(mask) => mask.test(atom),
        }
    }

    /// The selected positions, ascending.
    ///
    /// A concrete iterator rather than a boxed one: the boxed form cost a heap
    /// allocation per call and a virtual call per position, which prevented the
    /// consuming loop from inlining anything. Selection walks are the innermost
    /// loop of every spatial and analysis pass.
    #[must_use]
    pub fn iter(&self) -> SelectionIter<'_> {
        match self {
            Self::Empty => SelectionIter::Run(0..0),
            Self::All(count) => SelectionIter::Run(0..*count),
            Self::Range(run) => SelectionIter::Run(run.start..run.end),
            Self::Ranges(runs) => SelectionIter::Runs {
                runs,
                next: 0,
                current: 0..0,
            },
            Self::Sparse(positions) => SelectionIter::Sparse(positions.iter()),
            Self::Dense(mask) => SelectionIter::Dense(mask.ones()),
        }
    }
}

/// The positions of an [`AtomSelection`], ascending.
#[derive(Clone, Debug)]
pub enum SelectionIter<'a> {
    /// One contiguous run, which also covers the empty and whole-structure
    /// selections.
    Run(Range<u32>),
    /// Several runs, walked one after another.
    Runs {
        /// The remaining runs, including the one being walked.
        runs: &'a [Range<u32>],
        /// The index of the next run to start.
        next: usize,
        /// The unfinished part of the current run.
        current: Range<u32>,
    },
    /// Scattered positions read straight from their slice.
    Sparse(std::slice::Iter<'a, u32>),
    /// Set bits of a dense mask.
    Dense(Ones<'a>),
}

impl Iterator for SelectionIter<'_> {
    type Item = u32;

    fn next(&mut self) -> Option<u32> {
        match self {
            Self::Run(run) => run.next(),
            Self::Runs {
                runs,
                next,
                current,
            } => loop {
                if let Some(position) = current.next() {
                    return Some(position);
                }
                let run = runs.get(*next)?;
                *next = next.checked_add(1)?;
                *current = run.start..run.end;
            },
            Self::Sparse(positions) => positions.next().copied(),
            Self::Dense(ones) => ones.next(),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Self::Run(run) => run.size_hint(),
            Self::Runs {
                runs,
                next,
                current,
            } => {
                let remaining = match runs.get(*next..) {
                    Some(tail) => tail.iter().map(ExactSizeIterator::len).sum::<usize>(),
                    None => 0,
                };
                let total = current.len() + remaining;
                (total, Some(total))
            }
            Self::Sparse(positions) => positions.size_hint(),
            Self::Dense(ones) => ones.size_hint(),
        }
    }
}

impl ExactSizeIterator for SelectionIter<'_> {}

impl std::iter::FusedIterator for SelectionIter<'_> {}

impl FromIterator<u32> for AtomSelection {
    fn from_iter<T: IntoIterator<Item = u32>>(iter: T) -> Self {
        let mut positions: Vec<u32> = iter.into_iter().collect();
        positions.sort_unstable();
        positions.dedup();
        Self::from_sorted(positions)
    }
}

impl<'a> IntoIterator for &'a AtomSelection {
    type Item = u32;
    type IntoIter = SelectionIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

fn inline_runs(positions: &[u32]) -> Option<SmallVec<[Range<u32>; 4]>> {
    let mut runs = SmallVec::<[Range<u32>; 4]>::new();
    for position in positions.iter().copied() {
        let end = position.checked_add(1)?;
        if let Some(run) = runs.last_mut()
            && run.end == position
        {
            run.end = end;
            continue;
        }
        if runs.len() == runs.inline_size() {
            return None;
        }
        runs.push(position..end);
    }
    Some(runs)
}

fn mask_from_positions(len: u32, positions: &[u32]) -> BitVec {
    let mut mask = BitVec::repeat(false, len);
    for position in positions {
        mask.set(*position, true);
    }
    mask
}

fn selection_len(selection: &AtomSelection) -> u64 {
    match selection {
        AtomSelection::Empty => 0,
        AtomSelection::All(count) => u64::from(*count),
        AtomSelection::Range(run) => u64::from(run.end - run.start),
        AtomSelection::Ranges(runs) => runs.iter().map(|run| u64::from(run.end - run.start)).sum(),
        AtomSelection::Sparse(positions) => positions.len() as u64,
        AtomSelection::Dense(mask) => u64::from(mask.count_ones()),
    }
}
