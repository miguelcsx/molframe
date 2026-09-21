//! Explicit saturating conversions at analysis storage boundaries.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f32(value: usize) -> f32 {
    let Some(converted) = value.to_f32() else {
        return f32::INFINITY;
    };
    converted
}

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    finite_or_infinity(value.to_f64())
}

pub(crate) fn u64_to_f64(value: u64) -> f64 {
    finite_or_infinity(value.to_f64())
}

pub(crate) fn f64_to_f32(value: f64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}

pub(crate) fn f32_to_usize(value: f32) -> Option<usize> {
    value.to_usize()
}

pub(crate) fn usize_to_u32(value: usize) -> u32 {
    let Ok(converted) = u32::try_from(value) else {
        return u32::MAX;
    };
    converted
}

fn finite_or_infinity(value: Option<f64>) -> f64 {
    let Some(value) = value else {
        return f64::INFINITY;
    };
    value
}
