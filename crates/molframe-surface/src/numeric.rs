//! Explicit saturating conversions at public storage boundaries.

use num_traits::ToPrimitive;

pub(crate) use molframe_core::numeric::{
    f32_to_i64, f64_to_f32, i64_to_f64, usize_to_f32, usize_to_f64, usize_to_u32,
};

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

pub(crate) fn i128_to_i64(value: i128) -> i64 {
    match i64::try_from(value) {
        Ok(converted) => converted,
        Err(_) if value.is_negative() => i64::MIN,
        Err(_) => i64::MAX,
    }
}

pub(crate) fn usize_to_i32(value: usize) -> i32 {
    let Ok(converted) = i32::try_from(value) else {
        return i32::MAX;
    };
    converted
}

pub(crate) fn i32_to_usize(value: i32) -> usize {
    let Ok(converted) = usize::try_from(value) else {
        return 0;
    };
    converted
}

pub(crate) fn f32_to_i32(value: f32) -> i32 {
    match value.to_i32() {
        Some(converted) => converted,
        None if value.is_nan() => 0,
        None if value.is_sign_negative() => i32::MIN,
        None => i32::MAX,
    }
}

pub(crate) fn f32_to_usize(value: f32) -> usize {
    match value.to_usize() {
        Some(converted) => converted,
        None if value.is_nan() || value.is_sign_negative() => 0,
        None => usize::MAX,
    }
}
