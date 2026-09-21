//! Explicit saturating conversion at validation storage boundaries.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    let Some(converted) = value.to_f64() else {
        return f64::INFINITY;
    };
    converted
}

pub(crate) fn f64_to_f32(value: f64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}
