//! Allocation-bounded position iteration and sparse merging.

use super::representation::AtomSelection;
use crate::column::{BitVec, Ones};
use std::iter::Peekable;
use std::ops::Range;
use std::slice;

enum MergeCursor<'a> {
    Empty,
    Range(Range<u32>),
    Ranges {
        runs: &'a [Range<u32>],
        run: usize,
        next: u32,
    },
    Sparse(slice::Iter<'a, u32>),
    Dense(Ones<'a>),
}

impl<'a> MergeCursor<'a> {
    fn new(selection: &'a AtomSelection) -> Self {
        match selection {
            AtomSelection::Empty => Self::Empty,
            AtomSelection::All(count) => Self::Range(0..*count),
            AtomSelection::Range(run) => Self::Range(run.start..run.end),
            AtomSelection::Ranges(runs) => Self::Ranges {
                runs,
                run: 0,
                next: runs.first().map_or(0, |run| run.start),
            },
            AtomSelection::Sparse(positions) => Self::Sparse(positions.iter()),
            AtomSelection::Dense(mask) => Self::Dense(mask.ones()),
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
                *next = (*next).max(current.start);
                if *next < current.end {
                    let position = *next;
                    *next += 1;
                    return Some(position);
                }
                *run += 1;
                *next = runs.get(*run)?.start;
            },
        }
    }
}

pub(super) fn copy_selection(selection: &AtomSelection) -> AtomSelection {
    match selection {
        AtomSelection::Empty => AtomSelection::Empty,
        AtomSelection::All(count) => AtomSelection::All(*count),
        AtomSelection::Range(run) => AtomSelection::Range(run.start..run.end),
        AtomSelection::Ranges(runs) => {
            AtomSelection::Ranges(runs.iter().map(|run| run.start..run.end).collect())
        }
        AtomSelection::Sparse(positions) => AtomSelection::Sparse(positions.clone()),
        AtomSelection::Dense(mask) => AtomSelection::Dense(copy_mask(mask)),
    }
}

fn copy_mask(mask: &BitVec) -> BitVec {
    let mut copied = BitVec::repeat(false, mask.len());
    for position in mask.ones() {
        copied.set(position, true);
    }
    copied
}

pub(super) fn visit_positions(selection: &AtomSelection, mut visit: impl FnMut(u32)) {
    match selection {
        AtomSelection::Empty => {}
        AtomSelection::All(count) => (0..*count).for_each(visit),
        AtomSelection::Range(run) => (run.start..run.end).for_each(visit),
        AtomSelection::Ranges(runs) => {
            for run in runs {
                for position in run.start..run.end {
                    visit(position);
                }
            }
        }
        AtomSelection::Sparse(positions) => positions.iter().copied().for_each(visit),
        AtomSelection::Dense(mask) => mask.ones().for_each(visit),
    }
}

pub(super) fn merge_positions(left: &AtomSelection, right: &AtomSelection) -> AtomSelection {
    let mut merged = Vec::new();
    merge_cursors(
        MergeCursor::new(left).peekable(),
        MergeCursor::new(right).peekable(),
        &mut merged,
    );
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
