//! Explicit saturating conversions at fixture-construction boundaries.

use num_traits::ToPrimitive;

pub(crate) fn u64_to_f64(value: u64) -> f64 {
    finite_or_infinity(value.to_f64())
}

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    finite_or_infinity(value.to_f64())
}

fn finite_or_infinity(value: Option<f64>) -> f64 {
    let Some(value) = value else {
        return f64::INFINITY;
    };
    value
}

pub(crate) fn f64_to_f32(value: f64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}

pub(crate) fn f64_to_u64(value: f64) -> u64 {
    match value.to_u64() {
        Some(converted) => converted,
        None if value.is_sign_negative() => 0,
        None => u64::MAX,
    }
}
