//! Selection shapes, construction, membership and iteration.

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
    #[must_use]
    pub fn iter(&self) -> Box<dyn Iterator<Item = u32> + '_> {
        match self {
            Self::Empty => Box::new(std::iter::empty()),
            Self::All(count) => Box::new(0..*count),
            Self::Range(run) => Box::new(run.start..run.end),
            Self::Ranges(runs) => Box::new(runs.iter().flat_map(|run| run.start..run.end)),
            Self::Sparse(positions) => Box::new(positions.iter().copied()),
            Self::Dense(mask) => Box::new(mask.ones()),
        }
    }
}

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
    type IntoIter = Box<dyn Iterator<Item = u32> + 'a>;

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
