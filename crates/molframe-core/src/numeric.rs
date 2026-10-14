//! The one place a numeric conversion that can lose information is written.
//!
//! Every function here is total and states what it does at the edge of the
//! target type: integers become the nearest float, a float narrows to the
//! nearest float of the smaller width, and a float becomes an integer by
//! truncation toward zero, saturating at the type's bounds, with NaN as zero.
//! Callers use these in place of `as`, so no conversion in the workspace
//! needs a lint suppression, and each one is a thin wrapper that inlines away.

use num_traits::ToPrimitive;

/// A count as `f64`; exact up to 2^53, rounded to nearest beyond.
#[must_use]
#[inline]
pub fn usize_to_f64(value: usize) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None => f64::INFINITY,
    }
}

/// A count as `f64`; exact up to 2^53, rounded to nearest beyond.
#[must_use]
#[inline]
pub fn u64_to_f64(value: u64) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None => f64::INFINITY,
    }
}

/// A count as `f32`; exact up to 2^24, rounded to nearest beyond.
#[must_use]
#[inline]
pub fn usize_to_f32(value: usize) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None => f32::INFINITY,
    }
}

/// `f64` narrowed to the nearest `f32`; out-of-range magnitudes become
/// infinities and NaN stays NaN.
#[must_use]
#[inline]
pub fn f64_to_f32(value: f64) -> f32 {
    match value.to_f32() {
        Some(converted) => converted,
        None if value.is_sign_negative() => f32::NEG_INFINITY,
        None => f32::INFINITY,
    }
}

/// `f64` truncated toward zero into `i64`, saturating; NaN is zero.
#[must_use]
#[inline]
pub fn f64_to_i64(value: f64) -> i64 {
    match value.to_i64() {
        Some(converted) => converted,
        None if value.is_nan() => 0,
        None if value.is_sign_negative() => i64::MIN,
        None => i64::MAX,
    }
}

/// `f32` truncated toward zero into `i64`, saturating; NaN is zero.
#[must_use]
#[inline]
pub fn f32_to_i64(value: f32) -> i64 {
    f64_to_i64(f64::from(value))
}

/// `f64` truncated toward zero into `usize`, saturating; NaN and negative
/// values are zero.
#[must_use]
#[inline]
pub fn f64_to_usize(value: f64) -> usize {
    match value.to_usize() {
        Some(converted) => converted,
        None if value.is_nan() || value.is_sign_negative() => 0,
        None => usize::MAX,
    }
}

/// An index as `u32`, saturating at `u32::MAX`.
#[must_use]
#[inline]
pub fn usize_to_u32(value: usize) -> u32 {
    match u32::try_from(value) {
        Ok(converted) => converted,
        Err(_) => u32::MAX,
    }
}

/// A `u64` as `u32`, saturating at `u32::MAX`.
#[must_use]
#[inline]
pub fn u64_to_u32(value: u64) -> u32 {
    match u32::try_from(value) {
        Ok(converted) => converted,
        Err(_) => u32::MAX,
    }
}

/// `f64` truncated toward zero into `u32`, saturating; NaN and negative
/// values are zero.
#[must_use]
#[inline]
pub fn f64_to_u32(value: f64) -> u32 {
    match value.to_u32() {
        Some(converted) => converted,
        None if value.is_nan() || value.is_sign_negative() => 0,
        None => u32::MAX,
    }
}

/// `f64` truncated toward zero into `u64`, saturating; NaN and negative
/// values are zero.
#[must_use]
#[inline]
pub fn f64_to_u64(value: f64) -> u64 {
    match value.to_u64() {
        Some(converted) => converted,
        None if value.is_nan() || value.is_sign_negative() => 0,
        None => u64::MAX,
    }
}

/// The low eight bits of `value`; the rest are discarded on purpose.
#[must_use]
#[inline]
pub const fn low_u8(value: u32) -> u8 {
    value.to_le_bytes()[0]
}

/// The low sixteen bits of `value`; the rest are discarded on purpose.
#[must_use]
#[inline]
pub const fn low_u16(value: u32) -> u16 {
    let bytes = value.to_le_bytes();
    u16::from_le_bytes([bytes[0], bytes[1]])
}

/// IEEE equality of two `f32` values: `-0.0` equals `0.0`, NaN equals nothing.
#[must_use]
#[inline]
pub fn f32_eq(left: f32, right: f32) -> bool {
    left.partial_cmp(&right) == Some(std::cmp::Ordering::Equal)
}

/// IEEE equality of two `f64` values: `-0.0` equals `0.0`, NaN equals nothing.
#[must_use]
#[inline]
pub fn f64_eq(left: f64, right: f64) -> bool {
    left.partial_cmp(&right) == Some(std::cmp::Ordering::Equal)
}

/// An `i64` as `f64`; exact up to 2^53 in magnitude, rounded to nearest beyond.
#[must_use]
#[inline]
pub fn i64_to_f64(value: i64) -> f64 {
    match value.to_f64() {
        Some(converted) => converted,
        None => f64::NAN,
    }
}

/// A 24-bit sample `bits` scaled to `[0, 1)`; exact in `f32`.
#[must_use]
#[inline]
pub fn u64_to_f32_unit(bits: u64) -> f32 {
    f64_to_f32(u64_to_f64(bits & 0x00ff_ffff) / 16_777_216.0)
}
