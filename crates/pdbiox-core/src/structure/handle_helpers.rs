//! Small shared conversions used by hierarchy handles.

use crate::column::Presence;
use std::ops::Range;

pub(super) fn range_or_empty(range: Option<Range<u32>>) -> Range<u32> {
    match range {
        Some(range) => range,
        None => 0..0,
    }
}

pub(super) fn recorded_value<T>(entry: (T, Presence)) -> Option<T> {
    let (value, presence) = entry;
    presence.is_present().then_some(value)
}

pub(super) fn residue_count_as_u32(count: usize) -> u32 {
    match u32::try_from(count) {
        Ok(count) => count,
        Err(_) => u32::MAX,
    }
}
