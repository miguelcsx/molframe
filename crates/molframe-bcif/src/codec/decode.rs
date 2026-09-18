//! Complete `BinaryCIF` column codec chains.

use super::DecodedStringColumn;
use indexmap::IndexMap;
use molframe_core::diagnostic::{Code, Diagnostic};
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};
use serde_repr::{Deserialize_repr, Serialize_repr};
use std::borrow::Cow;
use std::sync::Arc;

/// `BinaryCIF` typed-array code.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize_repr, Deserialize_repr)]
#[repr(u8)]
pub enum DataType {
    /// Signed byte.
    Int8 = 1,
    /// Signed little-endian 16-bit integer.
    Int16 = 2,
    /// Signed little-endian 32-bit integer.
    Int32 = 3,
    /// Unsigned byte.
    Uint8 = 4,
    /// Unsigned little-endian 16-bit integer.
    Uint16 = 5,
    /// Unsigned little-endian 32-bit integer.
    Uint32 = 6,
    /// Little-endian IEEE 754 single precision.
    Float32 = 32,
    /// Little-endian IEEE 754 double precision.
    Float64 = 33,
}

/// One transformation in an encoded column's forward chain.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind")]
pub enum Encoding {
    /// Typed values stored as little-endian bytes.
    ByteArray {
        /// Element type.
        r#type: DataType,
    },
    /// Floats multiplied into signed integers.
    FixedPoint {
        /// Scale factor.
        factor: f64,
        /// Original float type.
        #[serde(rename = "srcType")]
        src_type: DataType,
    },
    /// Floats quantized uniformly over an interval.
    IntervalQuantization {
        /// Lower interval boundary.
        min: f64,
        /// Upper interval boundary.
        max: f64,
        /// Number of representable steps.
        #[serde(rename = "numSteps")]
        num_steps: u32,
        /// Original float type.
        #[serde(rename = "srcType")]
        src_type: DataType,
    },
    /// Integer value/count pairs.
    RunLength {
        /// Original integer type.
        #[serde(rename = "srcType")]
        src_type: DataType,
        /// Original element count.
        #[serde(rename = "srcSize")]
        src_size: usize,
    },
    /// Consecutive integer differences.
    Delta {
        /// Value immediately before the first stored delta.
        origin: i64,
        /// Original integer type.
        #[serde(rename = "srcType")]
        src_type: DataType,
    },
    /// Large integers split across signed or unsigned words.
    IntegerPacking {
        /// Packed word width, one or two bytes.
        #[serde(rename = "byteCount")]
        byte_count: u8,
        /// Whether packed words are unsigned.
        #[serde(rename = "isUnsigned")]
        is_unsigned: bool,
        /// Original element count.
        #[serde(rename = "srcSize")]
        src_size: usize,
    },
    /// Dictionary-encoded strings.
    StringArray {
        /// Encoding chain for dictionary indices.
        #[serde(rename = "dataEncoding")]
        data_encoding: Vec<Encoding>,
        /// Concatenated UTF-8 dictionary.
        #[serde(rename = "stringData")]
        string_data: String,
        /// Encoding chain for byte offsets.
        #[serde(rename = "offsetEncoding")]
        offset_encoding: Vec<Encoding>,
        /// Encoded byte offsets.
        #[serde(with = "serde_bytes")]
        offsets: Vec<u8>,
    },
}

/// A payload and the transformations that produced it.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EncodedData {
    /// Forward encoding chain.
    pub encoding: Vec<Encoding>,
    /// Final byte payload.
    #[serde(with = "serde_bytes")]
    pub data: Vec<u8>,
}

/// Fully decoded column values.
#[derive(Clone, Debug, PartialEq)]
pub enum Decoded {
    /// Integer values, widened without changing their value.
    Integers(Vec<i64>),
    /// Floating-point values.
    Floats(Vec<f64>),
    /// Dictionary-shared UTF-8 string values.
    Strings(DecodedStringColumn),
}

enum Stage<'a> {
    Bytes(Cow<'a, [u8]>),
    Integers(Vec<i64>),
    Floats(Vec<f64>),
    Strings(DecodedStringColumn),
}

/// Decodes a full chain in reverse order.
///
/// # Errors
///
/// Returns `E1401` for inconsistent lengths and `E1403` when adjacent codecs
/// disagree on their value type.
pub fn decode(encoded: &EncodedData) -> Result<Decoded, Diagnostic> {
    decode_borrowed(&encoded.encoding, &encoded.data)
}

pub(crate) fn decode_borrowed(encoding: &[Encoding], data: &[u8]) -> Result<Decoded, Diagnostic> {
    decode_stage(Stage::Bytes(Cow::Borrowed(data)), encoding)
}

fn decode_stage(mut stage: Stage<'_>, encoding: &[Encoding]) -> Result<Decoded, Diagnostic> {
    for encoding in encoding.iter().rev() {
        stage = decode_step(stage, encoding)?;
    }
    match stage {
        Stage::Integers(values) => Ok(Decoded::Integers(values)),
        Stage::Floats(values) => Ok(Decoded::Floats(values)),
        Stage::Strings(values) => Ok(Decoded::Strings(values)),
        Stage::Bytes(_) => Err(type_error("codec chain did not assign a type")),
    }
}

fn decode_step<'a>(stage: Stage<'a>, encoding: &Encoding) -> Result<Stage<'a>, Diagnostic> {
    match encoding {
        Encoding::ByteArray { r#type } => byte_array(stage, *r#type),
        Encoding::FixedPoint { factor, src_type } => fixed_point(stage, *factor, *src_type),
        Encoding::IntervalQuantization {
            min,
            max,
            num_steps,
            src_type,
        } => interval(stage, *min, *max, *num_steps, *src_type),
        Encoding::RunLength { src_type, src_size } => run_length(stage, *src_type, *src_size),
        Encoding::Delta { origin, src_type } => delta(stage, *origin, *src_type),
        Encoding::IntegerPacking {
            byte_count,
            is_unsigned,
            src_size,
        } => integer_packing(stage, *byte_count, *is_unsigned, *src_size),
        Encoding::StringArray {
            data_encoding,
            string_data,
            offset_encoding,
            offsets,
        } => string_array(stage, data_encoding, string_data, offset_encoding, offsets),
    }
}

fn byte_array(stage: Stage<'_>, data_type: DataType) -> Result<Stage<'_>, Diagnostic> {
    let Stage::Bytes(bytes) = stage else {
        return Err(type_error("ByteArray requires bytes"));
    };
    match data_type {
        DataType::Int8 => Ok(Stage::Integers(
            bytes
                .iter()
                .copied()
                .map(|value| i64::from(value.cast_signed()))
                .collect(),
        )),
        DataType::Uint8 => Ok(Stage::Integers(
            bytes.iter().copied().map(i64::from).collect(),
        )),
        DataType::Int16 => integers(&bytes, 2, |chunk| {
            i64::from(i16::from_le_bytes([chunk[0], chunk[1]]))
        }),
        DataType::Uint16 => integers(&bytes, 2, |chunk| {
            i64::from(u16::from_le_bytes([chunk[0], chunk[1]]))
        }),
        DataType::Int32 => integers(&bytes, 4, |chunk| {
            i64::from(i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        }),
        DataType::Uint32 => integers(&bytes, 4, |chunk| {
            i64::from(u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        }),
        DataType::Float32 => floats(&bytes, 4, |chunk| {
            f64::from(f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
        }),
        DataType::Float64 => floats(&bytes, 8, |chunk| {
            f64::from_le_bytes([
                chunk[0], chunk[1], chunk[2], chunk[3], chunk[4], chunk[5], chunk[6], chunk[7],
            ])
        }),
    }
}

fn integers<'a>(
    bytes: &[u8],
    width: usize,
    read: impl Fn(&[u8]) -> i64,
) -> Result<Stage<'a>, Diagnostic> {
    if !bytes.len().is_multiple_of(width) {
        return Err(length_error(bytes.len(), width));
    }
    Ok(Stage::Integers(
        bytes.chunks_exact(width).map(read).collect(),
    ))
}

fn floats<'a>(
    bytes: &[u8],
    width: usize,
    read: impl Fn(&[u8]) -> f64,
) -> Result<Stage<'a>, Diagnostic> {
    if !bytes.len().is_multiple_of(width) {
        return Err(length_error(bytes.len(), width));
    }
    Ok(Stage::Floats(bytes.chunks_exact(width).map(read).collect()))
}

fn fixed_point(stage: Stage<'_>, factor: f64, src_type: DataType) -> Result<Stage<'_>, Diagnostic> {
    let Stage::Integers(values) = stage else {
        return Err(type_error("FixedPoint requires integers"));
    };
    if !factor.is_finite() || factor == 0.0 || !is_float(src_type) {
        return Err(type_error("invalid FixedPoint metadata"));
    }
    let values = values
        .into_iter()
        .map(|value| {
            value
                .to_f64()
                .map(|value| value / factor)
                .ok_or_else(|| type_error("FixedPoint integer is outside the floating-point range"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Stage::Floats(values))
}

fn interval(
    stage: Stage<'_>,
    min: f64,
    max: f64,
    steps: u32,
    src_type: DataType,
) -> Result<Stage<'_>, Diagnostic> {
    let Stage::Integers(values) = stage else {
        return Err(type_error("IntervalQuantization requires integers"));
    };
    if steps < 2 || !min.is_finite() || !max.is_finite() || max < min || !is_float(src_type) {
        return Err(type_error("invalid IntervalQuantization metadata"));
    }
    let increment = (max - min) / f64::from(steps - 1);
    let values = values
        .into_iter()
        .map(|value| {
            value
                .to_f64()
                .map(|value| min + increment * value)
                .ok_or_else(|| type_error("quantized integer is outside the floating-point range"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Stage::Floats(values))
}

fn run_length(stage: Stage<'_>, src_type: DataType, size: usize) -> Result<Stage<'_>, Diagnostic> {
    let Stage::Integers(values) = stage else {
        return Err(type_error("RunLength requires integers"));
    };
    if values.len() % 2 != 0 || !is_integer(src_type) {
        return Err(length_error(values.len(), 2));
    }
    let mut output = Vec::with_capacity(size);
    for pair in values.as_chunks::<2>().0 {
        let Ok(count) = usize::try_from(pair[1]) else {
            return Err(length_error(values.len(), size));
        };
        let Some(new_length) = output.len().checked_add(count) else {
            return Err(type_error(
                "decoded run length exceeds the addressable range",
            ));
        };
        if new_length > size {
            return Err(length_error(new_length, size));
        }
        output.resize(new_length, pair[0]);
    }
    if output.len() != size {
        return Err(length_error(output.len(), size));
    }
    Ok(Stage::Integers(output))
}

fn delta(stage: Stage<'_>, origin: i64, src_type: DataType) -> Result<Stage<'_>, Diagnostic> {
    let Stage::Integers(mut values) = stage else {
        return Err(type_error("Delta requires integers"));
    };
    if !is_integer(src_type) {
        return Err(type_error("Delta source type is not integer"));
    }
    let mut previous = origin;
    let size = values.len();
    for value in &mut values {
        previous = previous
            .checked_add(*value)
            .ok_or_else(|| length_error(size, size))?;
        *value = previous;
    }
    Ok(Stage::Integers(values))
}

fn integer_packing(
    stage: Stage<'_>,
    byte_count: u8,
    unsigned: bool,
    size: usize,
) -> Result<Stage<'_>, Diagnostic> {
    let Stage::Integers(values) = stage else {
        return Err(type_error("IntegerPacking requires integers"));
    };
    if values.len() == size {
        return Ok(Stage::Integers(values));
    }
    let upper = match (byte_count, unsigned) {
        (1, true) => 255,
        (2, true) => 65_535,
        (1, false) => 127,
        (2, false) => 32_767,
        _ => return Err(type_error("IntegerPacking byteCount must be one or two")),
    };
    let lower = if unsigned { i64::MIN } else { -upper - 1 };
    let mut output = Vec::with_capacity(size);
    let mut total = 0i64;
    for value in values {
        total = total
            .checked_add(value)
            .ok_or_else(|| length_error(output.len(), size))?;
        if value != upper && value != lower {
            output.push(total);
            total = 0;
            if output.len() > size {
                return Err(length_error(output.len(), size));
            }
        }
    }
    if output.len() != size || total != 0 {
        return Err(length_error(output.len(), size));
    }
    Ok(Stage::Integers(output))
}

fn string_array<'a>(
    stage: Stage<'a>,
    data_encoding: &[Encoding],
    string_data: &str,
    offset_encoding: &[Encoding],
    offsets: &[u8],
) -> Result<Stage<'a>, Diagnostic> {
    let Stage::Bytes(data) = stage else {
        return Err(type_error("StringArray requires bytes"));
    };
    let Decoded::Integers(indices) = decode_stage(Stage::Bytes(data), data_encoding)? else {
        return Err(type_error("StringArray indices are not integers"));
    };
    let Decoded::Integers(offsets) =
        decode_stage(Stage::Bytes(Cow::Borrowed(offsets)), offset_encoding)?
    else {
        return Err(type_error("StringArray offsets are not integers"));
    };
    let source_capacity = match offsets.len().checked_sub(1) {
        Some(capacity) => capacity,
        None => 0,
    };
    let mut lookup = IndexMap::<&str, u32>::with_capacity(source_capacity);
    let mut dictionary = Vec::with_capacity(source_capacity);
    let mut source_to_compact = Vec::with_capacity(source_capacity);
    for pair in offsets.windows(2) {
        let (Ok(start), Ok(end)) = (usize::try_from(pair[0]), usize::try_from(pair[1])) else {
            return Err(length_error(offsets.len(), string_data.len()));
        };
        let Some(value) = string_data.get(start..end) else {
            return Err(length_error(end, string_data.len()));
        };
        let compact = intern_string(value, &mut lookup, &mut dictionary)?;
        source_to_compact.push(compact);
    }
    let mut output = Vec::with_capacity(indices.len());
    for index in indices {
        if index == -1 {
            output.push(intern_string("", &mut lookup, &mut dictionary)?);
            continue;
        }
        let Ok(index) = usize::try_from(index) else {
            return Err(length_error(output.len(), source_to_compact.len()));
        };
        let Some(compact) = source_to_compact.get(index) else {
            return Err(length_error(index, source_to_compact.len()));
        };
        output.push(*compact);
    }
    Ok(Stage::Strings(DecodedStringColumn::from_validated_parts(
        dictionary, output,
    )))
}

fn intern_string<'a>(
    value: &'a str,
    lookup: &mut IndexMap<&'a str, u32>,
    dictionary: &mut Vec<Arc<str>>,
) -> Result<u32, Diagnostic> {
    if let Some(index) = lookup.get(value) {
        return Ok(*index);
    }
    let Ok(index) = u32::try_from(dictionary.len()) else {
        return Err(
            Diagnostic::new(Code::E1401).with_message("string dictionary exceeds u32 indices")
        );
    };
    lookup.insert(value, index);
    dictionary.push(Arc::from(value));
    Ok(index)
}

const fn is_integer(data_type: DataType) -> bool {
    !matches!(data_type, DataType::Float32 | DataType::Float64)
}

const fn is_float(data_type: DataType) -> bool {
    matches!(data_type, DataType::Float32 | DataType::Float64)
}

fn length_error(actual: usize, declared: usize) -> Diagnostic {
    Diagnostic::new(Code::E1401)
        .with_context("actual", actual.to_string())
        .with_context("declared", declared.to_string())
}

fn type_error(message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1403).with_message(message)
}

#[cfg(test)]
#[path = "decode_tests.rs"]
mod tests;
