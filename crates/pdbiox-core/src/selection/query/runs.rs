//! Linear run algebra and bounded complements.

use super::representation::AtomSelection;
use crate::column::BitVec;
use smallvec::SmallVec;
use std::ops::Range;

pub(super) enum RunCursor<'a> {
    Single(Option<Range<u32>>),
    Many {
        runs: &'a [Range<u32>],
        position: usize,
    },
}

impl<'a> RunCursor<'a> {
    pub(super) fn new(selection: &'a AtomSelection) -> Option<Self> {
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

pub(super) fn intersect_runs(mut left: RunCursor<'_>, mut right: RunCursor<'_>) -> AtomSelection {
    let mut output = SmallVec::<[Range<u32>; 4]>::new();
    let mut left_run = left.next();
    let mut right_run = right.next();
    while let (Some(a), Some(b)) = (left_run.as_ref(), right_run.as_ref()) {
        push_run(&mut output, a.start.max(b.start)..a.end.min(b.end));
        let (a_end, b_end) = (a.end, b.end);
        if a_end <= b_end {
            left_run = left.next();
        }
        if b_end <= a_end {
            right_run = right.next();
        }
    }
    normalize_runs(output)
}

pub(super) fn union_runs(mut left: RunCursor<'_>, mut right: RunCursor<'_>) -> AtomSelection {
    let mut output = SmallVec::<[Range<u32>; 4]>::new();
    let mut left_run = left.next();
    let mut right_run = right.next();
    while let Some(run) = take_next(&mut left_run, &mut right_run, &mut left, &mut right) {
        push_run(&mut output, run);
    }
    normalize_runs(output)
}

pub(super) fn difference_runs(mut left: RunCursor<'_>, mut right: RunCursor<'_>) -> AtomSelection {
    let mut output = SmallVec::<[Range<u32>; 4]>::new();
    let mut right_run = right.next();
    for left_run in &mut left {
        let mut cursor = left_run.start;
        loop {
            let Some(current) = right_run.as_ref() else {
                break;
            };
            if current.end <= cursor {
                right_run = right.next();
                continue;
            }
            if current.start >= left_run.end {
                break;
            }
            if current.start > cursor {
                push_run(&mut output, cursor..current.start.min(left_run.end));
            }
            cursor = cursor.max(current.end);
            if cursor >= left_run.end {
                break;
            }
            right_run = right.next();
        }
        push_run(&mut output, cursor..left_run.end);
    }
    normalize_runs(output)
}

fn take_next(
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

pub(super) fn complement_within(selection: &AtomSelection, count: u32) -> AtomSelection {
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
    let mut gaps = SmallVec::<[Range<u32>; 4]>::new();
    let mut cursor = 0u32;
    for run in runs {
        if cursor >= count {
            break;
        }
        let (start, end) = (run.start.min(count), run.end.min(count));
        if end <= cursor {
            continue;
        }
        push_run(&mut gaps, cursor..start);
        cursor = cursor.max(end);
    }
    push_run(&mut gaps, cursor..count);
    normalize_runs(gaps)
}

fn complement_sparse(positions: &[u32], count: u32) -> AtomSelection {
    let mut gaps = SmallVec::<[Range<u32>; 4]>::new();
    let mut cursor = 0u32;
    for position in positions
        .iter()
        .copied()
        .take_while(|position| *position < count)
    {
        if position >= cursor {
            push_run(&mut gaps, cursor..position);
            cursor = position + 1;
        }
    }
    push_run(&mut gaps, cursor..count);
    normalize_runs(gaps)
}

fn complement_dense(selected: &BitVec, count: u32) -> AtomSelection {
    let mut complement = BitVec::repeat(true, count);
    for position in selected.ones().take_while(|position| *position < count) {
        complement.set(position, false);
    }
    if complement.none() {
        AtomSelection::Empty
    } else {
        AtomSelection::Dense(complement)
    }
}

pub(super) fn retain_after(selection: &AtomSelection, first: u32) -> AtomSelection {
    match selection {
        AtomSelection::Empty => AtomSelection::Empty,
        AtomSelection::All(count) => AtomSelection::range(first.min(*count)..*count),
        AtomSelection::Range(run) => AtomSelection::range(run.start.max(first)..run.end),
        AtomSelection::Ranges(runs) => {
            let mut retained = SmallVec::<[Range<u32>; 4]>::new();
            for run in runs {
                push_run(&mut retained, run.start.max(first)..run.end);
            }
            normalize_runs(retained)
        }
        AtomSelection::Sparse(positions) => {
            let start = positions.partition_point(|position| *position < first);
            AtomSelection::from_sorted(positions[start..].to_vec())
        }
        AtomSelection::Dense(mask) => {
            AtomSelection::from_sorted(mask.ones().filter(|position| *position >= first).collect())
        }
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
    } else {
        runs.push(run);
    }
}

fn normalize_runs(mut runs: SmallVec<[Range<u32>; 4]>) -> AtomSelection {
    match runs.len() {
        0 => AtomSelection::Empty,
        1 => runs
            .pop()
            .map_or(AtomSelection::Empty, AtomSelection::Range),
        _ => AtomSelection::Ranges(runs),
    }
}
