//! Explicit compact-index and count conversions for motif evaluation.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_u32(value: usize) -> u32 {
    let Ok(converted) = u32::try_from(value) else {
        return u32::MAX;
    };
    converted
}

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    let Some(converted) = value.to_f64() else {
        return f64::INFINITY;
    };
    converted
}
