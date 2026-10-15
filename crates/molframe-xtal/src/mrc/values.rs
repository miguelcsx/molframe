//! Numeric summaries and header encoding for MRC maps.

use super::MrcError;
use crate::numeric::{f64_to_f32, usize_to_f64};

pub(super) fn statistics(values: &[f32]) -> (f32, f32, f32, f32) {
    let mut minimum = f32::INFINITY;
    let mut maximum = f32::NEG_INFINITY;
    let mut sum = 0.0_f64;
    for value in values {
        minimum = minimum.min(*value);
        maximum = maximum.max(*value);
        sum += f64::from(*value);
    }
    let mean = sum / usize_to_f64(values.len());
    let variance = values
        .iter()
        .map(|value| (f64::from(*value) - mean).powi(2))
        .sum::<f64>()
        / usize_to_f64(values.len());
    (
        minimum,
        maximum,
        f64_to_f32(mean),
        f64_to_f32(variance.sqrt()),
    )
}

pub(super) fn to_i32(value: usize) -> Result<i32, MrcError> {
    i32::try_from(value).map_err(|_| MrcError::SizeOverflow)
}

fn to_f32(value: f64) -> Result<f32, MrcError> {
    if value.is_finite() && value.abs() <= f64::from(f32::MAX) {
        Ok(f64_to_f32(value))
    } else {
        Err(MrcError::InvalidHeader)
    }
}

pub(super) fn usize_triplet_to_i32(values: [usize; 3]) -> Result<[i32; 3], MrcError> {
    Ok([to_i32(values[0])?, to_i32(values[1])?, to_i32(values[2])?])
}

pub(super) fn f64_triplet_to_f32(values: [f64; 3]) -> Result<[f32; 3], MrcError> {
    Ok([to_f32(values[0])?, to_f32(values[1])?, to_f32(values[2])?])
}

pub(super) fn half_to_f32(value: u16) -> f32 {
    let sign = u32::from(value & 0x8000) << 16;
    let exponent = (value >> 10) & 0x1f;
    let fraction = u32::from(value & 0x03ff);
    let bits = match exponent {
        0 if fraction == 0 => sign,
        0 => {
            let shift = fraction.leading_zeros() - 21;
            sign | ((127 - 15 - shift) << 23) | ((fraction << (shift + 1) & 0x03ff) << 13)
        }
        31 => sign | 0x7f80_0000 | (fraction << 13),
        _ => sign | (u32::from(exponent + 112) << 23) | (fraction << 13),
    };
    f32::from_bits(bits)
}

pub(super) fn write_i32(bytes: &mut [u8], offset: usize, value: i32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

pub(super) fn write_f32(bytes: &mut [u8], offset: usize, value: f32) {
    bytes[offset..offset + 4].copy_from_slice(&value.to_le_bytes());
}

pub(super) fn write_i32_triplet(bytes: &mut [u8], offset: usize, values: [i32; 3]) {
    for (index, value) in values.into_iter().enumerate() {
        write_i32(bytes, offset + index * 4, value);
    }
}

pub(super) fn write_f32_triplet(bytes: &mut [u8], offset: usize, values: [f32; 3]) {
    for (index, value) in values.into_iter().enumerate() {
        write_f32(bytes, offset + index * 4, value);
    }
}
