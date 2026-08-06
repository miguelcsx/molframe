//! Sets of atoms, in whichever shape suits the set.
//!
//! A selection is never unconditionally a list and never unconditionally a bit
//! set. One chain is a contiguous run and costs two integers. Every protein atom
//! is dense and costs a bit each. Fifteen zinc atoms in a million-atom system
//! are a short list, and storing that as a bit set would spend a hundred and
//! twenty-five kilobytes to record fifteen numbers.
//!
//! Composition is proportional to the smaller operand. Intersecting a run with a
//! short list walks the list and tests the run, never the other way round —
//! which is what keeps a query that narrows quickly from paying for the breadth
//! it started with.
//!
//! Selections carry indices only. Nothing is materialised until a caller asks
//! for bytes.

use crate::column::BitVec;
use smallvec::SmallVec;
use std::ops::Range;

/// A set of atoms, addressed by position in the structure's flat atom order.
///
/// # Examples
///
/// ```
/// use pdbiox_core::AtomSelection;
///
/// let chain = AtomSelection::range(0..100);
/// let metals = AtomSelection::from_sorted(vec![5, 50, 500]);
///
/// assert_eq!(chain.intersect(&metals).len(), 2);
/// assert_eq!(chain.union(&metals).len(), 101);
/// assert_eq!(metals.difference(&chain).len(), 1);
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

/// Above this fraction of the span, a bit set beats a list of positions.
///
/// One member costs four bytes as a position and one bit as a mask, so the
/// crossover is at one member per thirty-two positions of span. The threshold
/// sits above that to leave the cheaper form in place for sets near the line.
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

    /// A selection over positions that are already ascending and unique.
    ///
    /// Chooses the shape by how the positions are laid out: one run stays a run,
    /// a handful of runs stay runs, a dense scatter becomes a mask, and anything
    /// else stays a list.
    #[must_use]
    pub fn from_sorted(positions: Vec<u32>) -> Self {
        let (Some(first), Some(last)) = (positions.first(), positions.last()) else {
            return Self::Empty;
        };
        let (first, last) = (*first, *last);
        let span = last - first + 1;
        if span as usize == positions.len() {
            return Self::Range(first..last + 1);
        }

        let mut runs: SmallVec<[Range<u32>; 4]> = SmallVec::new();
        for position in &positions {
            match runs.last_mut() {
                Some(run) if run.end == *position => run.end += 1,
                _ => {
                    if runs.len() == runs.inline_size() {
                        break;
                    }
                    runs.push(*position..*position + 1);
                }
            }
        }
        // Runs pay off only when they are actually runs. A run costs two
        // integers and a position costs one, so a set of isolated positions
        // stored as single-member runs would be twice the size of the list.
        let covered: u32 = runs.iter().map(|run| run.end - run.start).sum();
        if covered as usize == positions.len() && runs.len() * 2 < positions.len() {
            return Self::Ranges(runs);
        }

        if positions.len() as u32 * DENSE_FRACTION >= span {
            let mut mask = BitVec::repeat(false, last + 1);
            for position in &positions {
                mask.set(*position, true);
            }
            return Self::Dense(mask);
        }
        Self::Sparse(positions)
    }

    /// The number of atoms selected.
    #[must_use]
    pub fn len(&self) -> u32 {
        match self {
            Self::Empty => 0,
            Self::All(count) => *count,
            Self::Range(run) => run.end.saturating_sub(run.start),
            Self::Ranges(runs) => runs.iter().map(|run| run.end - run.start).sum(),
            Self::Sparse(positions) => positions.len() as u32,
            Self::Dense(mask) => mask.count_ones(),
        }
    }

    /// Returns true when nothing is selected.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns true when `atom` is selected.
    ///
    /// Constant time for a run or a mask, logarithmic for a list of positions or
    /// of runs.
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
    pub fn iter(&self) -> Box<dyn Iterator<Item = u32> + '_> {
        match self {
            Self::Empty => Box::new(std::iter::empty()),
            Self::All(count) => Box::new(0..*count),
            Self::Range(run) => Box::new(run.clone()),
            Self::Ranges(runs) => Box::new(runs.iter().flat_map(Clone::clone)),
            Self::Sparse(positions) => Box::new(positions.iter().copied()),
            Self::Dense(mask) => Box::new(mask.ones()),
        }
    }

    /// The atoms in both selections.
    ///
    /// Walks whichever selection is smaller and tests the other, so the cost
    /// follows the smaller operand rather than the larger.
    #[must_use]
    pub fn intersect(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Empty, _) | (_, Self::Empty) => return Self::Empty,
            (Self::All(_), _) => return other.clone(),
            (_, Self::All(_)) => return self.clone(),
            (Self::Range(left), Self::Range(right)) => {
                return Self::range(left.start.max(right.start)..left.end.min(right.end));
            }
            _ => {}
        }
        let (smaller, larger) = if self.len() <= other.len() {
            (self, other)
        } else {
            (other, self)
        };
        let kept: Vec<u32> = smaller
            .iter()
            .filter(|position| larger.contains(*position))
            .collect();
        Self::from_sorted(kept)
    }

    /// The atoms in either selection.
    #[must_use]
    pub fn union(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Empty, _) => return other.clone(),
            (_, Self::Empty) => return self.clone(),
            (all @ Self::All(_), _) | (_, all @ Self::All(_)) => return all.clone(),
            _ => {}
        }
        let mut merged = Vec::with_capacity((self.len() + other.len()) as usize);
        let (mut left, mut right) = (self.iter().peekable(), other.iter().peekable());
        loop {
            let next = match (left.peek(), right.peek()) {
                (None, None) => break,
                // A position present in both is emitted once.
                (Some(l), Some(r)) if l == r => {
                    right.next();
                    left.next()
                }
                (Some(l), Some(r)) if l > r => right.next(),
                (Some(_), _) => left.next(),
                (None, Some(_)) => right.next(),
            };
            if let Some(position) = next {
                merged.push(position);
            }
        }
        Self::from_sorted(merged)
    }

    /// The atoms in this selection but not the other.
    #[must_use]
    pub fn difference(&self, other: &Self) -> Self {
        if other.is_empty() {
            return self.clone();
        }
        let kept: Vec<u32> = self
            .iter()
            .filter(|position| !other.contains(*position))
            .collect();
        Self::from_sorted(kept)
    }

    /// Every atom of a structure of `count` atoms that this selection excludes.
    #[must_use]
    pub fn complement(&self, count: u32) -> Self {
        Self::All(count).difference(self)
    }

    /// The bytes this selection occupies beyond its own header.
    ///
    /// Reported so a caller can see that the shape chosen for a set is the
    /// cheap one; nothing depends on the value.
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        match self {
            Self::Empty | Self::All(_) | Self::Range(_) => 0,
            Self::Ranges(runs) => runs.len() * size_of::<Range<u32>>(),
            Self::Sparse(positions) => positions.len() * size_of::<u32>(),
            Self::Dense(mask) => (mask.len() as usize).div_ceil(8),
        }
    }
}

impl FromIterator<u32> for AtomSelection {
    /// Builds a selection from positions in any order, sorting and deduplicating.
    fn from_iter<T: IntoIterator<Item = u32>>(iter: T) -> Self {
        let mut positions: Vec<u32> = iter.into_iter().collect();
        positions.sort_unstable();
        positions.dedup();
        Self::from_sorted(positions)
    }
}

impl<'a> IntoIterator for &'a AtomSelection {
    type Item = u32;
    type IntoIter = Box<dyn Iterator<Item = u32> + 'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
#[path = "selection_tests.rs"]
mod tests;
