//! Exact floating-point encodings and explicit interval quantization.

use super::{integer::encode_integers, keep_smaller};
use crate::codec::{DataType, EncodedData, Encoding};
use pdbiox_core::diagnostic::{Code, Diagnostic};

/// Encodes floats without changing their IEEE values.
///
/// The optimiser tries `Float64`, exact `Float32`, and exact fixed-point powers
/// of ten. Lossy interval quantization is available only through
/// [`encode_interval`], where its precision is explicit.
///
/// # Errors
///
/// Returns a codec diagnostic if an exact integer candidate overflows.
pub fn encode_floats(values: &[f64]) -> Result<EncodedData, Diagnostic> {
    let mut best = byte_array_f64(values);
    if values
        .iter()
        .all(|value| f64::from(*value as f32).to_bits() == value.to_bits())
    {
        keep_smaller(&mut best, byte_array_f32(values));
    }
    for factor in [1.0, 10.0, 100.0, 1_000.0, 10_000.0, 100_000.0, 1_000_000.0] {
        if let Some(candidate) = fixed_point(values, factor)? {
            keep_smaller(&mut best, candidate);
        }
    }
    Ok(best)
}

/// Quantizes finite floats over their observed interval into `steps` bins.
///
/// # Errors
///
/// Returns a type diagnostic for non-finite input or fewer than two steps.
pub fn encode_interval(values: &[f64], steps: u32) -> Result<EncodedData, Diagnostic> {
    if steps < 2 || values.iter().any(|value| !value.is_finite()) {
        return Err(type_error(
            "interval quantization requires finite values and two steps",
        ));
    }
    let Some((&first, rest)) = values.split_first() else {
        return Ok(byte_array_f64(values));
    };
    let (mut minimum, mut maximum) = (first, first);
    for value in rest {
        minimum = minimum.min(*value);
        maximum = maximum.max(*value);
    }
    let increment = (maximum - minimum) / f64::from(steps - 1);
    let quantized = if increment == 0.0 {
        vec![0; values.len()]
    } else {
        values
            .iter()
            .map(|value| ((*value - minimum) / increment).round() as i64)
            .collect()
    };
    let mut encoded = encode_integers(&quantized)?;
    encoded.encoding.insert(
        0,
        Encoding::IntervalQuantization {
            min: minimum,
            max: maximum,
            num_steps: steps,
            src_type: DataType::Float64,
        },
    );
    Ok(encoded)
}

fn fixed_point(values: &[f64], factor: f64) -> Result<Option<EncodedData>, Diagnostic> {
    let mut integers = Vec::with_capacity(values.len());
    for value in values {
        let scaled = *value * factor;
        if !scaled.is_finite()
            || scaled.fract().abs().to_bits() != 0.0f64.to_bits()
            || scaled < f64::from(i32::MIN)
            || scaled > f64::from(u32::MAX)
            || (scaled / factor).to_bits() != value.to_bits()
        {
            return Ok(None);
        }
        integers.push(scaled as i64);
    }
    let mut encoded = encode_integers(&integers)?;
    encoded.encoding.insert(
        0,
        Encoding::FixedPoint {
            factor,
            src_type: DataType::Float64,
        },
    );
    Ok(Some(encoded))
}

fn byte_array_f32(values: &[f64]) -> EncodedData {
    let mut data = Vec::with_capacity(values.len().saturating_mul(4));
    for value in values {
        data.extend_from_slice(&(*value as f32).to_le_bytes());
    }
    EncodedData {
        encoding: vec![Encoding::ByteArray {
            r#type: DataType::Float32,
        }],
        data,
    }
}

fn byte_array_f64(values: &[f64]) -> EncodedData {
    let mut data = Vec::with_capacity(values.len().saturating_mul(8));
    for value in values {
        data.extend_from_slice(&value.to_le_bytes());
    }
    EncodedData {
        encoding: vec![Encoding::ByteArray {
            r#type: DataType::Float64,
        }],
        data,
    }
}

fn type_error(message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1403).with_message(message)
}

#[cfg(test)]
#[path = "float_tests.rs"]
mod tests;
