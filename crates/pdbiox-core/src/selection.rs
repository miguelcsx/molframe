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
use std::iter::Peekable;
use std::ops::Range;
use std::slice;

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
        let (Some(first), Some(last)) = (positions.first().copied(), positions.last().copied())
        else {
            return Self::Empty;
        };

        let span = u64::from(last)
            .saturating_sub(u64::from(first))
            .saturating_add(1);

        if last < u32::MAX && span == positions.len() as u64 {
            return Self::Range(first..last + 1);
        }

        if let Some(runs) = inline_runs(&positions) {
            // Runs pay off only when they are actually runs. A run costs two
            // integers and a position costs one, so a set of isolated positions
            // stored as single-member runs would be twice the size of the list.
            let covered = runs.iter().fold(0usize, |covered, run| {
                covered.saturating_add(run.end.saturating_sub(run.start) as usize)
            });

            if covered == positions.len() && runs.len() * 2 < positions.len() {
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
    #[must_use]
    pub fn len(&self) -> u32 {
        match self {
            Self::Empty => 0,
            Self::All(count) => *count,
            Self::Range(run) => run.end.saturating_sub(run.start),
            Self::Ranges(runs) => runs.iter().fold(0u32, |count, run| {
                count.saturating_add(run.end.saturating_sub(run.start))
            }),
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
            Self::Range(run) => Box::new(run.start..run.end),
            Self::Ranges(runs) => Box::new(runs.iter().flat_map(|run| run.start..run.end)),
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
            (Self::Empty, _) | (_, Self::Empty) => {
                return Self::Empty;
            }
            (Self::All(_), _) => {
                return copy_selection(other);
            }
            (_, Self::All(_)) => {
                return copy_selection(self);
            }
            _ => {}
        }

        if let (Some(left), Some(right)) = (RunCursor::new(self), RunCursor::new(other)) {
            return intersect_runs(left, right);
        }

        let (smaller, larger) = if self.len() <= other.len() {
            (self, other)
        } else {
            (other, self)
        };

        let mut kept = Vec::with_capacity(smaller.len() as usize);

        visit_positions(smaller, |position| {
            if larger.contains(position) {
                kept.push(position);
            }
        });

        Self::from_sorted(kept)
    }

    /// The atoms in either selection.
    #[must_use]
    pub fn union(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Empty, _) => {
                return copy_selection(other);
            }
            (_, Self::Empty) => {
                return copy_selection(self);
            }
            (Self::All(count), _) => {
                return Self::All(*count);
            }
            (_, Self::All(count)) => {
                return Self::All(*count);
            }
            _ => {}
        }

        if let (Some(left), Some(right)) = (RunCursor::new(self), RunCursor::new(other)) {
            return union_runs(left, right);
        }

        merge_positions(self, other)
    }

    /// The atoms in this selection but not the other.
    #[must_use]
    pub fn difference(&self, other: &Self) -> Self {
        if self.is_empty() {
            return Self::Empty;
        }

        if other.is_empty() {
            return copy_selection(self);
        }

        if let Self::All(count) = self {
            return complement_within(other, *count);
        }

        if let Self::All(count) = other {
            return retain_after(self, *count);
        }

        if let (Some(left), Some(right)) = (RunCursor::new(self), RunCursor::new(other)) {
            return difference_runs(left, right);
        }

        let mut kept = Vec::with_capacity(self.len() as usize);

        visit_positions(self, |position| {
            if !other.contains(position) {
                kept.push(position);
            }
        });

        Self::from_sorted(kept)
    }

    /// Every atom of a structure of `count` atoms that this selection excludes.
    #[must_use]
    pub fn complement(&self, count: u32) -> Self {
        if self.is_empty() {
            return Self::All(count);
        }

        if count == 0 {
            return Self::Empty;
        }

        complement_within(self, count)
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

enum RunCursor<'a> {
    Single(Option<Range<u32>>),
    Many {
        runs: &'a [Range<u32>],
        position: usize,
    },
}

impl<'a> RunCursor<'a> {
    fn new(selection: &'a AtomSelection) -> Option<Self> {
        match selection {
            AtomSelection::Range(run) => Some(Self::Single(Some(run.start..run.end))),
            AtomSelection::Ranges(runs) => Some(Self::Many { runs, position: 0 }),
            _ => None,
        }
    }
}

impl Iterator for RunCursor<'_> {
    type Item = Range<u32>;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Single(run) => run.take(),
            Self::Many { runs, position } => {
                let run = runs.get(*position)?;
                *position += 1;

                Some(run.start..run.end)
            }
        }
    }
}

enum MergeCursor<'a> {
    Empty,
    Range(Range<u32>),
    Ranges {
        runs: &'a [Range<u32>],
        run: usize,
        next: u32,
    },
    Sparse(slice::Iter<'a, u32>),
    Dense(Box<dyn Iterator<Item = u32> + 'a>),
}

impl<'a> MergeCursor<'a> {
    fn new(selection: &'a AtomSelection) -> Self {
        match selection {
            AtomSelection::Empty => Self::Empty,
            AtomSelection::All(count) => Self::Range(0..*count),
            AtomSelection::Range(run) => Self::Range(run.start..run.end),
            AtomSelection::Ranges(runs) => {
                let next = match runs.first() {
                    Some(run) => run.start,
                    None => 0,
                };

                Self::Ranges { runs, run: 0, next }
            }
            AtomSelection::Sparse(positions) => Self::Sparse(positions.iter()),
            AtomSelection::Dense(mask) => Self::Dense(Box::new(mask.ones())),
        }
    }
}

impl Iterator for MergeCursor<'_> {
    type Item = u32;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Empty => None,
            Self::Range(range) => range.next(),
            Self::Sparse(positions) => positions.next().copied(),
            Self::Dense(positions) => positions.next(),
            Self::Ranges { runs, run, next } => loop {
                let current = runs.get(*run)?;

                if *next < current.start {
                    *next = current.start;
                }

                if *next < current.end {
                    let position = *next;
                    *next = next.saturating_add(1);

                    return Some(position);
                }

                *run += 1;

                let following = runs.get(*run)?;
                *next = following.start;
            },
        }
    }
}

fn inline_runs(positions: &[u32]) -> Option<SmallVec<[Range<u32>; 4]>> {
    let mut runs: SmallVec<[Range<u32>; 4]> = SmallVec::new();

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

fn copy_selection(selection: &AtomSelection) -> AtomSelection {
    match selection {
        AtomSelection::Empty => AtomSelection::Empty,
        AtomSelection::All(count) => AtomSelection::All(*count),
        AtomSelection::Range(run) => AtomSelection::Range(run.start..run.end),
        AtomSelection::Ranges(runs) => {
            AtomSelection::Ranges(runs.iter().map(|run| run.start..run.end).collect())
        }
        AtomSelection::Sparse(positions) => {
            AtomSelection::Sparse(positions.iter().copied().collect())
        }
        AtomSelection::Dense(mask) => {
            let mut copied = BitVec::repeat(false, mask.len());

            for position in mask.ones() {
                copied.set(position, true);
            }

            AtomSelection::Dense(copied)
        }
    }
}

fn visit_positions(selection: &AtomSelection, mut visit: impl FnMut(u32)) {
    match selection {
        AtomSelection::Empty => {}
        AtomSelection::All(count) => {
            for position in 0..*count {
                visit(position);
            }
        }
        AtomSelection::Range(run) => {
            for position in run.start..run.end {
                visit(position);
            }
        }
        AtomSelection::Ranges(runs) => {
            for run in runs {
                for position in run.start..run.end {
                    visit(position);
                }
            }
        }
        AtomSelection::Sparse(positions) => {
            for position in positions {
                visit(*position);
            }
        }
        AtomSelection::Dense(mask) => {
            for position in mask.ones() {
                visit(position);
            }
        }
    }
}

fn merge_positions(left: &AtomSelection, right: &AtomSelection) -> AtomSelection {
    let capacity = left.len().saturating_add(right.len()) as usize;

    let mut merged = Vec::with_capacity(capacity);
    let left = MergeCursor::new(left).peekable();
    let right = MergeCursor::new(right).peekable();

    merge_cursors(left, right, &mut merged);

    AtomSelection::from_sorted(merged)
}

fn merge_cursors(
    mut left: Peekable<MergeCursor<'_>>,
    mut right: Peekable<MergeCursor<'_>>,
    merged: &mut Vec<u32>,
) {
    loop {
        let next = match (left.peek().copied(), right.peek().copied()) {
            (None, None) => break,
            (Some(left_position), Some(right_position)) if left_position == right_position => {
                right.next();
                left.next()
            }
            (Some(left_position), Some(right_position)) if left_position > right_position => {
                right.next()
            }
            (Some(_), _) => left.next(),
            (None, Some(_)) => right.next(),
        };

        if let Some(position) = next {
            merged.push(position);
        }
    }
}

fn intersect_runs(mut left: RunCursor<'_>, mut right: RunCursor<'_>) -> AtomSelection {
    let mut intersections: SmallVec<[Range<u32>; 4]> = SmallVec::new();

    let mut left_run = left.next();
    let mut right_run = right.next();

    loop {
        let (Some(current_left), Some(current_right)) = (left_run.as_ref(), right_run.as_ref())
        else {
            break;
        };

        let start = current_left.start.max(current_right.start);
        let end = current_left.end.min(current_right.end);
        let left_end = current_left.end;
        let right_end = current_right.end;

        if start < end {
            push_run(&mut intersections, start..end);
        }

        if left_end <= right_end {
            left_run = left.next();
        }

        if right_end <= left_end {
            right_run = right.next();
        }
    }

    normalize_runs(intersections)
}

fn union_runs(mut left: RunCursor<'_>, mut right: RunCursor<'_>) -> AtomSelection {
    let mut united: SmallVec<[Range<u32>; 4]> = SmallVec::new();

    let mut left_run = left.next();
    let mut right_run = right.next();

    while let Some(run) = take_next_run(&mut left_run, &mut right_run, &mut left, &mut right) {
        push_run(&mut united, run);
    }

    normalize_runs(united)
}

fn difference_runs(mut left: RunCursor<'_>, mut right: RunCursor<'_>) -> AtomSelection {
    let mut difference: SmallVec<[Range<u32>; 4]> = SmallVec::new();

    let mut right_run = right.next();

    for left_run in &mut left {
        let mut cursor = left_run.start;

        loop {
            let Some(current_right) = right_run.as_ref() else {
                break;
            };

            if current_right.end <= cursor {
                right_run = right.next();
                continue;
            }

            if current_right.start >= left_run.end {
                break;
            }

            if current_right.start > cursor {
                push_run(
                    &mut difference,
                    cursor..current_right.start.min(left_run.end),
                );
            }

            cursor = cursor.max(current_right.end);

            if cursor >= left_run.end {
                break;
            }

            right_run = right.next();
        }

        if cursor < left_run.end {
            push_run(&mut difference, cursor..left_run.end);
        }
    }

    normalize_runs(difference)
}

fn take_next_run(
    left_run: &mut Option<Range<u32>>,
    right_run: &mut Option<Range<u32>>,
    left: &mut RunCursor<'_>,
    right: &mut RunCursor<'_>,
) -> Option<Range<u32>> {
    let take_left = match (left_run.as_ref(), right_run.as_ref()) {
        (Some(left), Some(right)) => left.start <= right.start,
        (Some(_), None) => true,
        (None, Some(_)) => false,
        (None, None) => return None,
    };

    if take_left {
        let selected = left_run.take();
        *left_run = left.next();
        selected
    } else {
        let selected = right_run.take();
        *right_run = right.next();
        selected
    }
}

fn push_run(runs: &mut SmallVec<[Range<u32>; 4]>, run: Range<u32>) {
    if run.start >= run.end {
        return;
    }

    if let Some(last) = runs.last_mut()
        && run.start <= last.end
    {
        last.end = last.end.max(run.end);
        return;
    }

    runs.push(run);
}

fn normalize_runs(mut runs: SmallVec<[Range<u32>; 4]>) -> AtomSelection {
    match runs.len() {
        0 => AtomSelection::Empty,
        1 => match runs.pop() {
            Some(run) => AtomSelection::Range(run),
            None => AtomSelection::Empty,
        },
        _ => AtomSelection::Ranges(runs),
    }
}

fn complement_within(selection: &AtomSelection, count: u32) -> AtomSelection {
    match selection {
        AtomSelection::Empty => AtomSelection::All(count),
        AtomSelection::All(selected) => AtomSelection::range((*selected).min(count)..count),
        AtomSelection::Range(run) => complement_ranges(std::iter::once(run), count),
        AtomSelection::Ranges(runs) => complement_ranges(runs.iter(), count),
        AtomSelection::Sparse(positions) => complement_sparse(positions, count),
        AtomSelection::Dense(mask) => complement_dense(mask, count),
    }
}

fn complement_ranges<'a>(
    runs: impl IntoIterator<Item = &'a Range<u32>>,
    count: u32,
) -> AtomSelection {
    let mut gaps: SmallVec<[Range<u32>; 4]> = SmallVec::new();

    let mut cursor = 0u32;

    for run in runs {
        if cursor >= count {
            break;
        }

        let start = run.start.min(count);
        let end = run.end.min(count);

        if end <= cursor {
            continue;
        }

        if start > cursor {
            push_run(&mut gaps, cursor..start);
        }

        cursor = cursor.max(end);
    }

    if cursor < count {
        push_run(&mut gaps, cursor..count);
    }

    normalize_runs(gaps)
}

fn complement_sparse(positions: &[u32], count: u32) -> AtomSelection {
    let mut gaps: SmallVec<[Range<u32>; 4]> = SmallVec::new();

    let mut cursor = 0u32;

    for position in positions {
        if *position >= count {
            break;
        }

        if *position < cursor {
            continue;
        }

        if *position > cursor {
            push_run(&mut gaps, cursor..*position);
        }

        cursor = position.saturating_add(1);
    }

    if cursor < count {
        push_run(&mut gaps, cursor..count);
    }

    normalize_runs(gaps)
}

fn complement_dense(selected: &BitVec, count: u32) -> AtomSelection {
    let mut complement = BitVec::repeat(true, count);

    for position in selected.ones() {
        if position >= count {
            break;
        }

        complement.set(position, false);
    }

    if complement.none() {
        AtomSelection::Empty
    } else {
        AtomSelection::Dense(complement)
    }
}

fn retain_after(selection: &AtomSelection, first: u32) -> AtomSelection {
    match selection {
        AtomSelection::Empty => AtomSelection::Empty,
        AtomSelection::All(count) => AtomSelection::range(first.min(*count)..*count),
        AtomSelection::Range(run) => AtomSelection::range(run.start.max(first)..run.end),
        AtomSelection::Ranges(runs) => {
            let mut retained: SmallVec<[Range<u32>; 4]> = SmallVec::new();

            for run in runs {
                push_run(&mut retained, run.start.max(first)..run.end);
            }

            normalize_runs(retained)
        }
        AtomSelection::Sparse(positions) => {
            let start = positions.partition_point(|position| *position < first);

            AtomSelection::from_sorted(positions[start..].iter().copied().collect())
        }
        AtomSelection::Dense(mask) => {
            let positions: Vec<u32> = mask.ones().filter(|position| *position >= first).collect();

            AtomSelection::from_sorted(positions)
        }
    }
}

#[cfg(test)]
#[path = "selection_tests.rs"]
mod tests;
