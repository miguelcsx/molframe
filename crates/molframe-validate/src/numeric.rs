//! Explicit saturating conversion at validation storage boundaries.

use num_traits::ToPrimitive;

pub(crate) use molframe_core::numeric::{f64_to_f32, usize_to_f64};

pub(crate) fn f64_to_usize(value: f64) -> usize {
    match value.to_usize() {
        Some(converted) => converted,
        None => 0,
    }
}
