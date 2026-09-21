//! Explicit saturating conversions at public storage boundaries.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    let Some(converted) = value.to_f64() else {
        return f64::INFINITY;
    };
    converted
}

pub(crate) fn i64_to_f64(value: i64) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None if value.is_negative() => f64::NEG_INFINITY,
        None => f64::INFINITY,
    }
}

pub(crate) fn f64_to_f32(value: f64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}

pub(crate) fn f64_to_i64(value: f64) -> i64 {
    match value.to_i64() {
        Some(converted) => converted,
        None if value.is_sign_negative() => i64::MIN,
        None => i64::MAX,
    }
}

pub(crate) fn f64_to_usize(value: f64) -> usize {
    match value.to_usize() {
        Some(converted) => converted,
        None if value.is_sign_negative() => 0,
        None => usize::MAX,
    }
}

pub(crate) fn f64_to_u16(value: f64) -> u16 {
    match value.to_u16() {
        Some(converted) => converted,
        None if value.is_sign_negative() => 0,
        None => u16::MAX,
    }
}

pub(crate) fn usize_to_u32(value: usize) -> u32 {
    let Ok(converted) = u32::try_from(value) else {
        return u32::MAX;
    };
    converted
}

pub(crate) fn i128_to_i64(value: i128) -> i64 {
    match i64::try_from(value) {
        Ok(converted) => converted,
        Err(_) if value.is_negative() => i64::MIN,
        Err(_) => i64::MAX,
    }
}
