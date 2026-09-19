//! Explicit saturating conversion at tensor and table storage boundaries.

use num_traits::ToPrimitive;

pub(crate) fn usize_to_f64(value: usize) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None => f64::INFINITY,
    }
}

pub(crate) fn usize_to_u64(value: usize) -> u64 {
    match u64::try_from(value) {
        Ok(converted) => converted,
        Err(_) => u64::MAX,
    }
}

pub(crate) fn usize_to_u32(value: usize) -> u32 {
    match u32::try_from(value) {
        Ok(converted) => converted,
        Err(_) => u32::MAX,
    }
}

pub(crate) fn u64_to_usize(value: u64) -> usize {
    match usize::try_from(value) {
        Ok(converted) => converted,
        Err(_) => usize::MAX,
    }
}

pub(crate) fn f64_to_usize(value: f64) -> usize {
    match value.to_usize() {
        Some(converted) => converted,
        None if value.is_sign_negative() => 0,
        None => usize::MAX,
    }
}

pub(crate) fn f64_to_f32(value: f64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}

pub(crate) fn u32_to_f32(value: u32) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
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
