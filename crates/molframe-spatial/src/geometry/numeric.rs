//! Checked numeric conversions shared by spatial kernels.

use num_traits::ToPrimitive;

pub(crate) fn floor_usize(value: f64) -> Option<usize> {
    value.is_finite().then(|| value.floor())?.to_usize()
}

pub(crate) fn ceil_i64(value: f64) -> Option<i64> {
    value.is_finite().then(|| value.ceil())?.to_i64()
}

pub(crate) fn usize_f64(value: usize) -> Option<f64> {
    const MAX_EXACT_F64_INTEGER: usize = if usize::BITS <= f64::MANTISSA_DIGITS {
        usize::MAX
    } else {
        1usize << f64::MANTISSA_DIGITS
    };
    (value <= MAX_EXACT_F64_INTEGER)
        .then_some(value)
        .and_then(|number| number.to_f64())
}

pub(crate) fn i64_f64(value: i64) -> Option<f64> {
    const MAX_EXACT_F64_INTEGER: u64 = 1u64 << f64::MANTISSA_DIGITS;
    (value.unsigned_abs() <= MAX_EXACT_F64_INTEGER)
        .then_some(value)
        .and_then(|number| number.to_f64())
}

pub(crate) fn finite_f32(value: f64) -> Option<f32> {
    let limit = f64::from(f32::MAX);
    (value.is_finite() && value >= -limit && value <= limit)
        .then_some(value)
        .and_then(|number| number.to_f32())
}

pub(crate) fn f64_f32(value: f64) -> f32 {
    if value.is_nan() {
        return f32::NAN;
    }
    let limit = f64::from(f32::MAX);
    if value > limit {
        return f32::INFINITY;
    }
    if value < -limit {
        return f32::NEG_INFINITY;
    }
    match value.to_f32() {
        Some(converted) => converted,
        None => f32::NAN,
    }
}

pub(crate) fn rounded_i64(value: f64) -> i64 {
    match value.to_i64() {
        Some(converted) => converted,
        None if value.is_sign_negative() => i64::MIN,
        None => i64::MAX,
    }
}
