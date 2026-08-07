//! Public set algebra over adaptive selection shapes.

use super::merge::{copy_selection, merge_positions, visit_positions};
use super::representation::AtomSelection;
use super::runs::{
    RunCursor, complement_within, difference_runs, intersect_runs, retain_after, union_runs,
};
use std::ops::Range;

impl AtomSelection {
    /// The atoms in both selections.
    #[must_use]
    pub fn intersect(&self, other: &Self) -> Self {
        match (self, other) {
            (Self::Empty, _) | (_, Self::Empty) => return Self::Empty,
            (Self::All(_), _) => return copy_selection(other),
            (_, Self::All(_)) => return copy_selection(self),
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
            (Self::Empty, _) => return copy_selection(other),
            (_, Self::Empty) => return copy_selection(self),
            (Self::All(count), _) | (_, Self::All(count)) => return Self::All(*count),
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
