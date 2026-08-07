//! Integer transformations and smallest-candidate selection.

use super::keep_smaller;
use crate::codec::{DataType, EncodedData, Encoding};
use pdbiox_core::diagnostic::{Code, Diagnostic};

/// Encodes integers by trying literal, delta, run-length and packed forms.
///
/// Candidate order is stable and ties keep the earlier representation.
///
/// # Errors
///
/// Returns a length diagnostic if a difference or encoded size overflows.
pub fn encode_integers(values: &[i64]) -> Result<EncodedData, Diagnostic> {
    if !representable(values) {
        return Err(encoding_error());
    }
    let source_type = source_type(values);
    let Some(mut best) = plain(values) else {
        return Err(encoding_error());
    };
    if let Some(candidate) = packed(values, best.data.len())? {
        keep_smaller(&mut best, candidate);
    }

    let differences = differences(values)?;
    if let Some(delta) = transformed(
        &differences,
        Encoding::Delta {
            origin: 0,
            src_type: source_type,
        },
    )? {
        keep_smaller(&mut best, delta);
    }

    let runs = runs(values)?;
    if let Some(run_length) = transformed(
        &runs,
        Encoding::RunLength {
            src_type: source_type,
            src_size: values.len(),
        },
    )? {
        keep_smaller(&mut best, run_length);
    }
    Ok(best)
}

fn transformed(values: &[i64], encoding: Encoding) -> Result<Option<EncodedData>, Diagnostic> {
    let Some(mut candidate) = plain(values) else {
        return Ok(None);
    };
    if let Some(packed) = packed(values, candidate.data.len())? {
        keep_smaller(&mut candidate, packed);
    }
    candidate.encoding.insert(0, encoding);
    Ok(Some(candidate))
}

fn plain(values: &[i64]) -> Option<EncodedData> {
    let data_type = source_type(values);
    if !values_fit(values, data_type) {
        return None;
    }
    let data = match data_type {
        DataType::Int8 => values.iter().map(|value| *value as i8 as u8).collect(),
        DataType::Uint8 => values.iter().map(|value| *value as u8).collect(),
        DataType::Int16 => flatten(values, |value| (*value as i16).to_le_bytes()),
        DataType::Uint16 => flatten(values, |value| (*value as u16).to_le_bytes()),
        DataType::Int32 => flatten(values, |value| (*value as i32).to_le_bytes()),
        DataType::Uint32 => flatten(values, |value| (*value as u32).to_le_bytes()),
        DataType::Float32 | DataType::Float64 => return None,
    };
    Some(EncodedData {
        encoding: vec![Encoding::ByteArray { r#type: data_type }],
        data,
    })
}

fn values_fit(values: &[i64], data_type: DataType) -> bool {
    values.iter().all(|value| match data_type {
        DataType::Int8 => i8::try_from(*value).is_ok(),
        DataType::Int16 => i16::try_from(*value).is_ok(),
        DataType::Int32 => i32::try_from(*value).is_ok(),
        DataType::Uint8 => u8::try_from(*value).is_ok(),
        DataType::Uint16 => u16::try_from(*value).is_ok(),
        DataType::Uint32 => u32::try_from(*value).is_ok(),
        DataType::Float32 | DataType::Float64 => false,
    })
}

fn packed(values: &[i64], maximum_bytes: usize) -> Result<Option<EncodedData>, Diagnostic> {
    let unsigned = values.iter().all(|value| *value >= 0);
    let one_bytes = packed_word_count(values, 1, unsigned)?;
    let two_bytes = packed_word_count(values, 2, unsigned)?.saturating_mul(2);
    let (byte_count, encoded_bytes) = if one_bytes <= two_bytes {
        (1, one_bytes)
    } else {
        (2, two_bytes)
    };
    if encoded_bytes >= maximum_bytes {
        return Ok(None);
    }
    let words = pack(values, byte_count, unsigned)?;
    let data_type = match (byte_count, unsigned) {
        (1, true) => DataType::Uint8,
        (1, false) => DataType::Int8,
        (2, true) => DataType::Uint16,
        _ => DataType::Int16,
    };
    let data = if byte_count == 1 {
        words.iter().map(|value| *value as i8 as u8).collect()
    } else {
        flatten(&words, |value| (*value as i16).to_le_bytes())
    };
    Ok(Some(EncodedData {
        encoding: vec![
            Encoding::IntegerPacking {
                byte_count,
                is_unsigned: unsigned,
                src_size: values.len(),
            },
            Encoding::ByteArray { r#type: data_type },
        ],
        data,
    }))
}

fn packed_word_count(values: &[i64], bytes: u8, unsigned: bool) -> Result<usize, Diagnostic> {
    let upper = packing_upper(bytes, unsigned)?;
    let lower = -upper - 1;
    let mut count = 0usize;
    for value in values {
        let continuations = if *value >= upper {
            value / upper
        } else if !unsigned && *value <= lower {
            value / lower
        } else {
            0
        };
        let Ok(continuations) = usize::try_from(continuations) else {
            return Err(encoding_error());
        };
        count = count
            .checked_add(continuations.saturating_add(1))
            .ok_or_else(encoding_error)?;
    }
    Ok(count)
}

fn pack(values: &[i64], bytes: u8, unsigned: bool) -> Result<Vec<i64>, Diagnostic> {
    let upper = packing_upper(bytes, unsigned)?;
    let lower = -upper - 1;
    let mut output = Vec::new();
    for original in values {
        let mut value = *original;
        while value >= upper {
            output.push(upper);
            value -= upper;
        }
        if !unsigned {
            while value <= lower {
                output.push(lower);
                value -= lower;
            }
        }
        output.push(value);
    }
    Ok(output)
}

fn packing_upper(bytes: u8, unsigned: bool) -> Result<i64, Diagnostic> {
    match (bytes, unsigned) {
        (1, true) => Ok(255),
        (2, true) => Ok(65_535),
        (1, false) => Ok(127),
        (2, false) => Ok(32_767),
        _ => Err(encoding_error()),
    }
}

fn representable(values: &[i64]) -> bool {
    let minimum = values.iter().copied().min();
    let maximum = values.iter().copied().max();
    match (minimum, maximum) {
        (Some(minimum), Some(maximum)) => {
            (minimum >= 0 && maximum <= i64::from(u32::MAX))
                || (minimum >= i64::from(i32::MIN) && maximum <= i64::from(i32::MAX))
        }
        _ => true,
    }
}

fn differences(values: &[i64]) -> Result<Vec<i64>, Diagnostic> {
    let mut previous = 0i64;
    let mut output = Vec::with_capacity(values.len());
    for value in values {
        let Some(delta) = value.checked_sub(previous) else {
            return Err(encoding_error());
        };
        output.push(delta);
        previous = *value;
    }
    Ok(output)
}

fn runs(values: &[i64]) -> Result<Vec<i64>, Diagnostic> {
    let mut output = Vec::new();
    let Some(mut current) = values.first().copied() else {
        return Ok(output);
    };
    let mut count = 0usize;
    for value in values {
        if *value == current {
            count = count.saturating_add(1);
            continue;
        }
        push_run(&mut output, current, count)?;
        current = *value;
        count = 1;
    }
    push_run(&mut output, current, count)?;
    Ok(output)
}

fn push_run(output: &mut Vec<i64>, value: i64, count: usize) -> Result<(), Diagnostic> {
    let Ok(count) = i64::try_from(count) else {
        return Err(encoding_error());
    };
    output.extend([value, count]);
    Ok(())
}

fn source_type(values: &[i64]) -> DataType {
    let minimum = match values.iter().copied().min() {
        Some(value) => value,
        None => 0,
    };
    let maximum = match values.iter().copied().max() {
        Some(value) => value,
        None => 0,
    };
    if minimum >= 0 && maximum <= i64::from(u8::MAX) {
        DataType::Uint8
    } else if minimum >= i64::from(i8::MIN) && maximum <= i64::from(i8::MAX) {
        DataType::Int8
    } else if minimum >= 0 && maximum <= i64::from(u16::MAX) {
        DataType::Uint16
    } else if minimum >= i64::from(i16::MIN) && maximum <= i64::from(i16::MAX) {
        DataType::Int16
    } else if minimum >= 0 && maximum <= i64::from(u32::MAX) {
        DataType::Uint32
    } else {
        DataType::Int32
    }
}

fn flatten<const N: usize>(values: &[i64], encode: impl Fn(&i64) -> [u8; N]) -> Vec<u8> {
    let mut output = Vec::with_capacity(values.len().saturating_mul(N));
    for value in values {
        output.extend_from_slice(&encode(value));
    }
    output
}

fn encoding_error() -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message("integer column cannot be represented")
}

#[cfg(test)]
#[path = "integer_tests.rs"]
mod tests;
