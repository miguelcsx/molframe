//! Integer transformations and smallest-candidate selection.

use super::strategy::{EncodingChoice, encoded_size_lower_bound};
use crate::codec::{DataType, EncodedData, Encoding};
use molframe_core::diagnostic::{Code, Diagnostic};

/// Encodes integers by trying literal, delta, run-length and packed forms.
///
/// Candidate order is stable and ties keep the earlier representation.
///
/// # Errors
///
/// Returns a length diagnostic if a difference or encoded size overflows.
pub fn encode_integers(values: &[i64]) -> Result<EncodedData, Diagnostic> {
    // One pass summarises the column; every type decision below reads that
    // summary instead of walking the values again.
    let range = Range::of(values);
    if !range.is_representable() {
        return Err(encoding_error());
    }
    let source_type = range.source_type();
    let Some(plain) = plain(values, range) else {
        return Err(encoding_error());
    };
    let maximum_bytes = plain.data.len();
    let mut best = EncodingChoice::new(plain);
    if let Some(candidate) = packed(values, range, maximum_bytes)? {
        best.consider(candidate);
    }

    let delta_encoding = Encoding::Delta {
        origin: 0,
        src_type: source_type,
    };
    if can_outperform(&delta_encoding, values.len(), best.size()) {
        let differences = differences(values)?;
        if let Some(delta) = transformed(&differences, delta_encoding)? {
            best.consider(delta);
        }
    }

    let run_encoding = Encoding::RunLength {
        src_type: source_type,
        src_size: values.len(),
    };
    if let Some(runs) = runs_if_competitive(values, &run_encoding, best.size())?
        && let Some(run_length) = transformed(&runs, run_encoding)?
    {
        best.consider(run_length);
    }
    Ok(best.finish())
}

fn transformed(values: &[i64], encoding: Encoding) -> Result<Option<EncodedData>, Diagnostic> {
    let range = Range::of(values);
    let Some(plain) = plain(values, range) else {
        return Ok(None);
    };
    let maximum_bytes = plain.data.len();
    let mut candidate = EncodingChoice::new(plain);
    if let Some(packed) = packed(values, range, maximum_bytes)? {
        candidate.consider(packed);
    }
    let mut candidate = candidate.finish();
    candidate.encoding.insert(0, encoding);
    Ok(Some(candidate))
}

fn plain(values: &[i64], range: Range) -> Option<EncodedData> {
    let data_type = range.source_type();
    if !range.fits(data_type) {
        return None;
    }
    let data = match data_type {
        DataType::Int8 => checked_flatten(values, |value| {
            i8::try_from(value).ok().map(i8::to_le_bytes)
        })?,
        DataType::Uint8 => checked_flatten(values, |value| {
            u8::try_from(value).ok().map(u8::to_le_bytes)
        })?,
        DataType::Int16 => checked_flatten(values, |value| {
            i16::try_from(value).ok().map(i16::to_le_bytes)
        })?,
        DataType::Uint16 => checked_flatten(values, |value| {
            u16::try_from(value).ok().map(u16::to_le_bytes)
        })?,
        DataType::Int32 => checked_flatten(values, |value| {
            i32::try_from(value).ok().map(i32::to_le_bytes)
        })?,
        DataType::Uint32 => checked_flatten(values, |value| {
            u32::try_from(value).ok().map(u32::to_le_bytes)
        })?,
        DataType::Float32 | DataType::Float64 => return None,
    };
    Some(EncodedData {
        encoding: vec![Encoding::ByteArray { r#type: data_type }],
        data,
    })
}

fn packed(
    values: &[i64],
    range: Range,
    maximum_bytes: usize,
) -> Result<Option<EncodedData>, Diagnostic> {
    let unsigned = range.is_unsigned();
    let (one_words, two_words) = packed_word_counts(values, unsigned)?;
    let two_bytes = two_words.checked_mul(2).ok_or_else(encoding_error)?;
    let (byte_count, encoded_bytes) = if one_words <= two_bytes {
        (1, one_words)
    } else {
        (2, two_bytes)
    };
    if encoded_bytes >= maximum_bytes {
        return Ok(None);
    }
    let data_type = match (byte_count, unsigned) {
        (1, true) => DataType::Uint8,
        (1, false) => DataType::Int8,
        (2, true) => DataType::Uint16,
        _ => DataType::Int16,
    };
    let data = pack_bytes(values, byte_count, unsigned, encoded_bytes)?;
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

/// Word counts for one-byte and two-byte packing, in a single pass.
///
/// Both widths are always compared before one is chosen, so counting them
/// together halves the traffic over the column.
fn packed_word_counts(values: &[i64], unsigned: bool) -> Result<(usize, usize), Diagnostic> {
    let one_upper = packing_upper(1, unsigned)?;
    let two_upper = packing_upper(2, unsigned)?;
    let one_lower = -one_upper - 1;
    let two_lower = -two_upper - 1;
    let mut one_count = 0usize;
    let mut two_count = 0usize;

    for value in values {
        one_count = one_count
            .checked_add(words_for(*value, one_upper, one_lower, unsigned)?)
            .ok_or_else(encoding_error)?;
        two_count = two_count
            .checked_add(words_for(*value, two_upper, two_lower, unsigned)?)
            .ok_or_else(encoding_error)?;
    }
    Ok((one_count, two_count))
}

/// Words one value occupies at a given packing width.
fn words_for(value: i64, upper: i64, lower: i64, unsigned: bool) -> Result<usize, Diagnostic> {
    let continuations = if value >= upper {
        value / upper
    } else if !unsigned && value <= lower {
        value / lower
    } else {
        0
    };
    let Ok(continuations) = usize::try_from(continuations) else {
        return Err(encoding_error());
    };
    continuations.checked_add(1).ok_or_else(encoding_error)
}

fn pack_bytes(
    values: &[i64],
    bytes: u8,
    unsigned: bool,
    encoded_bytes: usize,
) -> Result<Vec<u8>, Diagnostic> {
    let upper = packing_upper(bytes, unsigned)?;
    let lower = -upper - 1;
    let mut output = Vec::with_capacity(encoded_bytes);
    for original in values {
        let mut value = *original;
        while value >= upper {
            push_packed_word(&mut output, upper, bytes, unsigned)?;
            value -= upper;
        }
        if !unsigned {
            while value <= lower {
                push_packed_word(&mut output, lower, bytes, unsigned)?;
                value -= lower;
            }
        }
        push_packed_word(&mut output, value, bytes, unsigned)?;
    }
    Ok(output)
}

fn push_packed_word(
    output: &mut Vec<u8>,
    value: i64,
    bytes: u8,
    unsigned: bool,
) -> Result<(), Diagnostic> {
    match (bytes, unsigned) {
        (1, true) => output.extend_from_slice(
            &u8::try_from(value)
                .map_err(|_| encoding_error())?
                .to_le_bytes(),
        ),
        (1, false) => output.extend_from_slice(
            &i8::try_from(value)
                .map_err(|_| encoding_error())?
                .to_le_bytes(),
        ),
        (2, true) => output.extend_from_slice(
            &u16::try_from(value)
                .map_err(|_| encoding_error())?
                .to_le_bytes(),
        ),
        (2, false) => output.extend_from_slice(
            &i16::try_from(value)
                .map_err(|_| encoding_error())?
                .to_le_bytes(),
        ),
        _ => return Err(encoding_error()),
    }
    Ok(())
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

/// The smallest and largest value of a column.
///
/// Representability, the narrowest source type, whether every value fits that
/// type, and whether the column is unsigned are all functions of these two
/// numbers alone. Computing them once replaces four full scans of the column
/// with one, before any encoding candidate is even built.
#[derive(Clone, Copy, Debug)]
struct Range {
    /// The smallest value, or zero for an empty column.
    minimum: i64,
    /// The largest value, or zero for an empty column.
    maximum: i64,
}

impl Range {
    /// Summarises a column in a single pass.
    fn of(values: &[i64]) -> Self {
        let Some(first) = values.first() else {
            return Self {
                minimum: 0,
                maximum: 0,
            };
        };
        let mut minimum = *first;
        let mut maximum = *first;
        for value in values {
            if *value < minimum {
                minimum = *value;
            }
            if *value > maximum {
                maximum = *value;
            }
        }
        Self { minimum, maximum }
    }

    /// Whether one integer type covers the whole column.
    ///
    /// A column spanning a negative minimum and a maximum above `i32::MAX` fits
    /// neither `Int32` nor `Uint32`, and `BinaryCIF` has no wider integer type.
    fn is_representable(self) -> bool {
        (self.minimum >= 0 && self.maximum <= i64::from(u32::MAX))
            || (self.minimum >= i64::from(i32::MIN) && self.maximum <= i64::from(i32::MAX))
    }

    /// Whether every value fits `data_type`.
    fn fits(self, data_type: DataType) -> bool {
        match data_type {
            DataType::Int8 => {
                self.minimum >= i64::from(i8::MIN) && self.maximum <= i64::from(i8::MAX)
            }
            DataType::Int16 => {
                self.minimum >= i64::from(i16::MIN) && self.maximum <= i64::from(i16::MAX)
            }
            DataType::Int32 => {
                self.minimum >= i64::from(i32::MIN) && self.maximum <= i64::from(i32::MAX)
            }
            DataType::Uint8 => self.minimum >= 0 && self.maximum <= i64::from(u8::MAX),
            DataType::Uint16 => self.minimum >= 0 && self.maximum <= i64::from(u16::MAX),
            DataType::Uint32 => self.minimum >= 0 && self.maximum <= i64::from(u32::MAX),
            DataType::Float32 | DataType::Float64 => false,
        }
    }

    /// Whether the column holds no negative value.
    fn is_unsigned(self) -> bool {
        self.minimum >= 0
    }

    /// The narrowest type that covers the column.
    fn source_type(self) -> DataType {
        if self.fits(DataType::Uint8) {
            DataType::Uint8
        } else if self.fits(DataType::Int8) {
            DataType::Int8
        } else if self.fits(DataType::Uint16) {
            DataType::Uint16
        } else if self.fits(DataType::Int16) {
            DataType::Int16
        } else if self.fits(DataType::Uint32) {
            DataType::Uint32
        } else {
            DataType::Int32
        }
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

fn runs_if_competitive(
    values: &[i64],
    encoding: &Encoding,
    current_size: usize,
) -> Result<Option<Vec<i64>>, Diagnostic> {
    let run_count = values
        .windows(2)
        .filter(|pair| pair[0] != pair[1])
        .count()
        .checked_add(usize::from(!values.is_empty()))
        .ok_or_else(encoding_error)?;
    let element_count = run_count.checked_mul(2).ok_or_else(encoding_error)?;
    if !can_outperform(encoding, element_count, current_size) {
        return Ok(None);
    }
    let mut output = Vec::with_capacity(element_count);
    let Some(mut current) = values.first().copied() else {
        return Ok(Some(output));
    };
    let mut count = 0usize;
    for value in values {
        if *value == current {
            count = count.checked_add(1).ok_or_else(encoding_error)?;
            continue;
        }
        push_run(&mut output, current, count)?;
        current = *value;
        count = 1;
    }
    push_run(&mut output, current, count)?;
    Ok(Some(output))
}

fn can_outperform(encoding: &Encoding, element_count: usize, current_size: usize) -> bool {
    let minimum_encoding = vec![
        encoding.clone(),
        Encoding::ByteArray {
            r#type: DataType::Uint8,
        },
    ];
    match encoded_size_lower_bound(minimum_encoding, element_count) {
        Some(lower_bound) => lower_bound < current_size,
        None => true,
    }
}

fn push_run(output: &mut Vec<i64>, value: i64, count: usize) -> Result<(), Diagnostic> {
    let Ok(count) = i64::try_from(count) else {
        return Err(encoding_error());
    };
    output.extend([value, count]);
    Ok(())
}

fn checked_flatten<const N: usize>(
    values: &[i64],
    encode: impl Fn(i64) -> Option<[u8; N]>,
) -> Option<Vec<u8>> {
    let capacity = values.len().checked_mul(N)?;
    let mut output = Vec::with_capacity(capacity);
    for value in values {
        output.extend_from_slice(&encode(*value)?);
    }
    Some(output)
}

fn encoding_error() -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message("integer column cannot be represented")
}

#[cfg(test)]
#[path = "integer_tests.rs"]
mod tests;
