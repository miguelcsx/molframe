//! Explicit saturating conversions at crystallographic storage boundaries.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None => f64::INFINITY,
    }
}

pub(crate) fn usize_to_u32(value: usize) -> u32 {
    match u32::try_from(value) {
        Ok(converted) => converted,
        Err(_) => u32::MAX,
    }
}

pub(crate) fn f64_to_f32(value: f64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}

pub(crate) fn f64_to_i32(value: f64) -> i32 {
    match value.to_i32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => i32::MIN,
        None => i32::MAX,
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

pub(crate) fn f32_to_i64(value: f32) -> i64 {
    match value.to_i64() {
        Some(converted) => converted,
        None if value.is_sign_negative() => i64::MIN,
        None => i64::MAX,
    }
}

pub(crate) fn i64_to_f32(value: i64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}

pub(crate) fn i64_to_f64(value: i64) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None if value.is_negative() => f64::NEG_INFINITY,
        None => f64::INFINITY,
    }
}

pub(crate) fn i32_to_usize(value: i32) -> usize {
    match usize::try_from(value) {
        Ok(converted) => converted,
        Err(_) => 0,
    }
}
