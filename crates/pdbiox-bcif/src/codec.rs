//! Complete `BinaryCIF` column codec chains.

use pdbiox_core::diagnostic::{Code, Diagnostic};
use serde::{Deserialize, Serialize};
use serde_repr::{Deserialize_repr, Serialize_repr};

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
    /// UTF-8 string values.
    Strings(Vec<String>),
}

enum Stage {
    Bytes(Vec<u8>),
    Integers(Vec<i64>),
    Floats(Vec<f64>),
    Strings(Vec<String>),
}

/// Decodes a full chain in reverse order.
///
/// # Errors
///
/// Returns `E1401` for inconsistent lengths and `E1403` when adjacent codecs
/// disagree on their value type.
pub fn decode(encoded: &EncodedData) -> Result<Decoded, Diagnostic> {
    let mut stage = Stage::Bytes(encoded.data.clone());
    for encoding in encoded.encoding.iter().rev() {
        stage = decode_step(stage, encoding)?;
    }
    match stage {
        Stage::Integers(values) => Ok(Decoded::Integers(values)),
        Stage::Floats(values) => Ok(Decoded::Floats(values)),
        Stage::Strings(values) => Ok(Decoded::Strings(values)),
        Stage::Bytes(_) => Err(type_error("codec chain did not assign a type")),
    }
}

fn decode_step(stage: Stage, encoding: &Encoding) -> Result<Stage, Diagnostic> {
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

fn byte_array(stage: Stage, data_type: DataType) -> Result<Stage, Diagnostic> {
    let Stage::Bytes(bytes) = stage else {
        return Err(type_error("ByteArray requires bytes"));
    };
    match data_type {
        DataType::Int8 => Ok(Stage::Integers(
            bytes
                .into_iter()
                .map(|value| i64::from(value.cast_signed()))
                .collect(),
        )),
        DataType::Uint8 => Ok(Stage::Integers(bytes.into_iter().map(i64::from).collect())),
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

fn integers(bytes: &[u8], width: usize, read: impl Fn(&[u8]) -> i64) -> Result<Stage, Diagnostic> {
    if !bytes.len().is_multiple_of(width) {
        return Err(length_error(bytes.len(), width));
    }
    Ok(Stage::Integers(
        bytes.chunks_exact(width).map(read).collect(),
    ))
}

fn floats(bytes: &[u8], width: usize, read: impl Fn(&[u8]) -> f64) -> Result<Stage, Diagnostic> {
    if !bytes.len().is_multiple_of(width) {
        return Err(length_error(bytes.len(), width));
    }
    Ok(Stage::Floats(bytes.chunks_exact(width).map(read).collect()))
}

fn fixed_point(stage: Stage, factor: f64, src_type: DataType) -> Result<Stage, Diagnostic> {
    let Stage::Integers(values) = stage else {
        return Err(type_error("FixedPoint requires integers"));
    };
    if !factor.is_finite() || factor == 0.0 || !is_float(src_type) {
        return Err(type_error("invalid FixedPoint metadata"));
    }
    Ok(Stage::Floats(
        values
            .into_iter()
            .map(|value| value as f64 / factor)
            .collect(),
    ))
}

fn interval(
    stage: Stage,
    min: f64,
    max: f64,
    steps: u32,
    src_type: DataType,
) -> Result<Stage, Diagnostic> {
    let Stage::Integers(values) = stage else {
        return Err(type_error("IntervalQuantization requires integers"));
    };
    if steps < 2 || !min.is_finite() || !max.is_finite() || max < min || !is_float(src_type) {
        return Err(type_error("invalid IntervalQuantization metadata"));
    }
    let increment = (max - min) / f64::from(steps - 1);
    Ok(Stage::Floats(
        values
            .into_iter()
            .map(|value| min + increment * value as f64)
            .collect(),
    ))
}

fn run_length(stage: Stage, src_type: DataType, size: usize) -> Result<Stage, Diagnostic> {
    let Stage::Integers(values) = stage else {
        return Err(type_error("RunLength requires integers"));
    };
    if values.len() % 2 != 0 || !is_integer(src_type) {
        return Err(length_error(values.len(), 2));
    }
    let mut output = Vec::with_capacity(size);
    for pair in values.chunks_exact(2) {
        let Ok(count) = usize::try_from(pair[1]) else {
            return Err(length_error(values.len(), size));
        };
        if output.len().saturating_add(count) > size {
            return Err(length_error(output.len().saturating_add(count), size));
        }
        output.resize(output.len() + count, pair[0]);
    }
    if output.len() != size {
        return Err(length_error(output.len(), size));
    }
    Ok(Stage::Integers(output))
}

fn delta(stage: Stage, origin: i64, src_type: DataType) -> Result<Stage, Diagnostic> {
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
    stage: Stage,
    byte_count: u8,
    unsigned: bool,
    size: usize,
) -> Result<Stage, Diagnostic> {
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

fn string_array(
    stage: Stage,
    data_encoding: &[Encoding],
    string_data: &str,
    offset_encoding: &[Encoding],
    offsets: &[u8],
) -> Result<Stage, Diagnostic> {
    let Stage::Bytes(data) = stage else {
        return Err(type_error("StringArray requires bytes"));
    };
    let Decoded::Integers(indices) = decode(&EncodedData {
        encoding: data_encoding.to_vec(),
        data,
    })?
    else {
        return Err(type_error("StringArray indices are not integers"));
    };
    let Decoded::Integers(offsets) = decode(&EncodedData {
        encoding: offset_encoding.to_vec(),
        data: offsets.to_vec(),
    })?
    else {
        return Err(type_error("StringArray offsets are not integers"));
    };
    let mut dictionary = Vec::with_capacity(offsets.len().saturating_sub(1));
    for pair in offsets.windows(2) {
        let (Ok(start), Ok(end)) = (usize::try_from(pair[0]), usize::try_from(pair[1])) else {
            return Err(length_error(offsets.len(), string_data.len()));
        };
        let Some(value) = string_data.get(start..end) else {
            return Err(length_error(end, string_data.len()));
        };
        dictionary.push(value);
    }
    let mut output = Vec::with_capacity(indices.len());
    for index in indices {
        if index == -1 {
            output.push(String::new());
            continue;
        }
        let Ok(index) = usize::try_from(index) else {
            return Err(length_error(output.len(), dictionary.len()));
        };
        let Some(value) = dictionary.get(index) else {
            return Err(length_error(index, dictionary.len()));
        };
        output.push((*value).to_owned());
    }
    Ok(Stage::Strings(output))
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
#[path = "codec_tests.rs"]
mod tests;
