//! Explicit saturating conversions at analysis storage boundaries.

use num_traits::ToPrimitive;

pub(crate) use molframe_core::numeric::{
    f64_to_f32, u64_to_f64, usize_to_f32, usize_to_f64, usize_to_u32,
};

pub(crate) fn f32_to_usize(value: f32) -> Option<usize> {
    value.to_usize()
}
