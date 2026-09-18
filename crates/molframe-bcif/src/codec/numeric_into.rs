//! Allocation-free numeric decode into caller-owned final columns.

use super::{DataType, EncodedData, Encoding};
use molframe_core::{Code, Diagnostic};
use num_traits::ToPrimitive;

#[cfg(target_arch = "x86_64")]
use wide::f64x4 as F64xN;
#[cfg(target_arch = "x86_64")]
const LANES: usize = 4;
#[cfg(not(target_arch = "x86_64"))]
use wide::f64x2 as F64xN;
#[cfg(not(target_arch = "x86_64"))]
const LANES: usize = 2;

/// Decodes a numeric `BinaryCIF` column directly into caller-owned `f32` storage.
///
/// Literal floats and the integer, packing, delta, run-length, fixed-point and
/// interval chains emitted by molframe are fused into one forward pass. No
/// decoded staging vector is created.
///
/// # Errors
///
/// Returns a codec diagnostic for an unsupported chain, malformed payload or
/// output whose length differs from the decoded column.
pub fn decode_f32_into(encoded: &EncodedData, output: &mut [f32]) -> Result<usize, Diagnostic> {
    decode_f32_borrowed_into(&encoded.encoding, &encoded.data, output)
}

pub(crate) fn decode_f32_borrowed_into(
    encoding: &[Encoding],
    data: &[u8],
    output: &mut [f32],
) -> Result<usize, Diagnostic> {
    let Some((last, prefix)) = encoding.split_last() else {
        return Err(type_error("numeric column has no ByteArray encoding"));
    };
    let Encoding::ByteArray { r#type } = last else {
        return Err(type_error("numeric column does not end in ByteArray"));
    };
    if prefix.is_empty() && matches!(r#type, DataType::Float32 | DataType::Float64) {
        decode_literal_float(data, *r#type, output)?;
        return Ok(output.len());
    }
    if !is_integer(*r#type) {
        return Err(type_error("numeric transform requires integer bytes"));
    }
    let (conversion, integer_chain) = split_conversion(prefix)?;
    let mut reader = IntegerReader::new(data, *r#type, integer_chain)?;
    match reader.expansion {
        Expansion::Plain => decode_plain(&mut reader, conversion, output)?,
        Expansion::Delta(origin) => decode_delta(&mut reader, origin, conversion, output)?,
        Expansion::RunLength(size) => {
            if size != output.len() {
                return Err(length_error("run-length size differs from output"));
            }
            decode_runs(&mut reader, conversion, output)?;
        }
    }
    if reader.next_value()?.is_some() {
        return Err(length_error("encoded payload has trailing numeric values"));
    }
    Ok(output.len())
}

#[derive(Clone, Copy)]
enum Conversion {
    Identity,
    FixedPoint(f64),
    Interval { minimum: f64, increment: f64 },
}

fn split_conversion(prefix: &[Encoding]) -> Result<(Conversion, &[Encoding]), Diagnostic> {
    let Some(first) = prefix.first() else {
        return Ok((Conversion::Identity, prefix));
    };
    match first {
        Encoding::FixedPoint { factor, src_type } => {
            if !factor.is_finite() || *factor == 0.0 || !is_float(*src_type) {
                return Err(type_error("invalid FixedPoint metadata"));
            }
            Ok((Conversion::FixedPoint(*factor), &prefix[1..]))
        }
        Encoding::IntervalQuantization {
            min,
            max,
            num_steps,
            src_type,
        } => {
            if *num_steps < 2
                || !min.is_finite()
                || !max.is_finite()
                || max < min
                || !is_float(*src_type)
            {
                return Err(type_error("invalid IntervalQuantization metadata"));
            }
            Ok((
                Conversion::Interval {
                    minimum: *min,
                    increment: (max - min) / f64::from(*num_steps - 1),
                },
                &prefix[1..],
            ))
        }
        _ => Ok((Conversion::Identity, prefix)),
    }
}

#[derive(Clone, Copy)]
enum Expansion {
    Plain,
    Delta(i64),
    RunLength(usize),
}

struct IntegerReader<'a> {
    bytes: &'a [u8],
    data_type: DataType,
    word: usize,
    packing: Option<(i64, i64)>,
    expansion: Expansion,
}

impl<'a> IntegerReader<'a> {
    fn new(bytes: &'a [u8], data_type: DataType, chain: &[Encoding]) -> Result<Self, Diagnostic> {
        let (expansion, rest) = match chain.first() {
            Some(Encoding::Delta { origin, src_type }) if is_integer(*src_type) => {
                (Expansion::Delta(*origin), &chain[1..])
            }
            Some(Encoding::RunLength { src_type, src_size }) if is_integer(*src_type) => {
                (Expansion::RunLength(*src_size), &chain[1..])
            }
            Some(Encoding::Delta { .. } | Encoding::RunLength { .. }) => {
                return Err(type_error(
                    "integer transform has a non-integer source type",
                ));
            }
            _ => (Expansion::Plain, chain),
        };
        let packing = match rest {
            [] => None,
            [
                Encoding::IntegerPacking {
                    byte_count,
                    is_unsigned,
                    ..
                },
            ] => Some(packing_bounds(*byte_count, *is_unsigned)?),
            _ => return Err(type_error("unsupported numeric codec chain")),
        };
        let width = data_width(data_type);
        if !bytes.len().is_multiple_of(width) {
            return Err(length_error("ByteArray payload has an incomplete value"));
        }
        Ok(Self {
            bytes,
            data_type,
            word: 0,
            packing,
            expansion,
        })
    }

    fn next_value(&mut self) -> Result<Option<i64>, Diagnostic> {
        let Some((upper, lower)) = self.packing else {
            return self.next_word();
        };
        let Some(mut value) = self.next_word()? else {
            return Ok(None);
        };
        let mut total = 0_i64;
        loop {
            total = total
                .checked_add(value)
                .ok_or_else(|| length_error("packed integer overflow"))?;
            if value != upper && value != lower {
                return Ok(Some(total));
            }
            value = self
                .next_word()?
                .ok_or_else(|| length_error("unterminated packed integer"))?;
        }
    }

    fn next_word(&mut self) -> Result<Option<i64>, Diagnostic> {
        let width = data_width(self.data_type);
        let start = self
            .word
            .checked_mul(width)
            .ok_or_else(|| length_error("ByteArray offset overflow"))?;
        let Some(bytes) = self.bytes.get(start..start.saturating_add(width)) else {
            return Ok(None);
        };
        self.word = self.word.saturating_add(1);
        read_word(bytes, self.data_type).map(Some)
    }
}

fn decode_plain(
    reader: &mut IntegerReader<'_>,
    conversion: Conversion,
    output: &mut [f32],
) -> Result<(), Diagnostic> {
    let (blocks, tail) = output.as_chunks_mut::<LANES>();
    for block in blocks {
        let values = read_block(reader, "numeric payload is shorter than output")?;
        block.copy_from_slice(&convert_block(values, conversion)?);
    }
    for target in tail {
        *target = convert(
            reader
                .next_value()?
                .ok_or_else(|| length_error("numeric payload is shorter than output"))?,
            conversion,
        )?;
    }
    Ok(())
}

fn decode_delta(
    reader: &mut IntegerReader<'_>,
    origin: i64,
    conversion: Conversion,
    output: &mut [f32],
) -> Result<(), Diagnostic> {
    let mut previous = origin;
    let (blocks, tail) = output.as_chunks_mut::<LANES>();
    for block in blocks {
        let mut values = [0_i64; LANES];
        for value in &mut values {
            let delta = reader
                .next_value()?
                .ok_or_else(|| length_error("delta payload is shorter than output"))?;
            previous = previous
                .checked_add(delta)
                .ok_or_else(|| length_error("delta integer overflow"))?;
            *value = previous;
        }
        block.copy_from_slice(&convert_block(values, conversion)?);
    }
    for target in tail {
        let delta = reader
            .next_value()?
            .ok_or_else(|| length_error("delta payload is shorter than output"))?;
        previous = previous
            .checked_add(delta)
            .ok_or_else(|| length_error("delta integer overflow"))?;
        *target = convert(previous, conversion)?;
    }
    Ok(())
}

fn read_block(
    reader: &mut IntegerReader<'_>,
    short: &'static str,
) -> Result<[i64; LANES], Diagnostic> {
    let mut values = [0_i64; LANES];
    for value in &mut values {
        *value = reader.next_value()?.ok_or_else(|| length_error(short))?;
    }
    Ok(values)
}

fn convert_block(values: [i64; LANES], conversion: Conversion) -> Result<[f32; LANES], Diagnostic> {
    let mut floats = [0.0_f64; LANES];
    for (target, value) in floats.iter_mut().zip(values) {
        *target = value
            .to_f64()
            .ok_or_else(|| type_error("integer cannot be represented as a float"))?;
    }
    let values = F64xN::from(floats);
    let converted = match conversion {
        Conversion::Identity => values,
        Conversion::FixedPoint(factor) => values / F64xN::splat(factor),
        Conversion::Interval { minimum, increment } => {
            F64xN::splat(minimum) + F64xN::splat(increment) * values
        }
    };
    let mut output = [0.0_f32; LANES];
    for (target, value) in output.iter_mut().zip(converted.to_array()) {
        *target = value
            .to_f32()
            .ok_or_else(|| type_error("decoded value cannot be represented as f32"))?;
    }
    Ok(output)
}

fn decode_runs(
    reader: &mut IntegerReader<'_>,
    conversion: Conversion,
    output: &mut [f32],
) -> Result<(), Diagnostic> {
    let mut written = 0usize;
    while written < output.len() {
        let value = reader
            .next_value()?
            .ok_or_else(|| length_error("run-length payload is shorter than output"))?;
        let count = reader
            .next_value()?
            .and_then(|value| usize::try_from(value).ok())
            .ok_or_else(|| length_error("run length is not addressable"))?;
        let end = written
            .checked_add(count)
            .filter(|end| *end <= output.len())
            .ok_or_else(|| length_error("run length exceeds output"))?;
        output[written..end].fill(convert(value, conversion)?);
        written = end;
    }
    Ok(())
}

fn decode_literal_float(
    bytes: &[u8],
    data_type: DataType,
    output: &mut [f32],
) -> Result<(), Diagnostic> {
    let width = data_width(data_type);
    if bytes.len() != output.len().saturating_mul(width) {
        return Err(length_error("literal float length differs from output"));
    }
    for (source, target) in bytes.chunks_exact(width).zip(output) {
        *target = match data_type {
            DataType::Float32 => f32::from_le_bytes([source[0], source[1], source[2], source[3]]),
            DataType::Float64 => f64::from_le_bytes([
                source[0], source[1], source[2], source[3], source[4], source[5], source[6],
                source[7],
            ])
            .to_f32()
            .ok_or_else(|| type_error("float64 value cannot be represented as f32"))?,
            _ => return Err(type_error("literal numeric payload is not floating point")),
        };
    }
    Ok(())
}

fn convert(value: i64, conversion: Conversion) -> Result<f32, Diagnostic> {
    let value = value
        .to_f64()
        .ok_or_else(|| type_error("integer cannot be represented as a float"))?;
    let value = match conversion {
        Conversion::Identity => value,
        Conversion::FixedPoint(factor) => value / factor,
        Conversion::Interval { minimum, increment } => minimum + increment * value,
    };
    value
        .to_f32()
        .ok_or_else(|| type_error("decoded value cannot be represented as f32"))
}

fn read_word(bytes: &[u8], data_type: DataType) -> Result<i64, Diagnostic> {
    match data_type {
        DataType::Int8 => Ok(i64::from(bytes[0].cast_signed())),
        DataType::Uint8 => Ok(i64::from(bytes[0])),
        DataType::Int16 => Ok(i64::from(i16::from_le_bytes([bytes[0], bytes[1]]))),
        DataType::Uint16 => Ok(i64::from(u16::from_le_bytes([bytes[0], bytes[1]]))),
        DataType::Int32 => Ok(i64::from(i32::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3],
        ]))),
        DataType::Uint32 => Ok(i64::from(u32::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3],
        ]))),
        DataType::Float32 | DataType::Float64 => {
            Err(type_error("integer reader received floating-point bytes"))
        }
    }
}

const fn data_width(data_type: DataType) -> usize {
    match data_type {
        DataType::Int8 | DataType::Uint8 => 1,
        DataType::Int16 | DataType::Uint16 => 2,
        DataType::Int32 | DataType::Uint32 | DataType::Float32 => 4,
        DataType::Float64 => 8,
    }
}

fn packing_bounds(bytes: u8, unsigned: bool) -> Result<(i64, i64), Diagnostic> {
    let upper = match (bytes, unsigned) {
        (1, true) => 255,
        (2, true) => 65_535,
        (1, false) => 127,
        (2, false) => 32_767,
        _ => return Err(type_error("IntegerPacking byteCount must be one or two")),
    };
    let lower = if unsigned { i64::MIN } else { -upper - 1 };
    Ok((upper, lower))
}

const fn is_integer(data_type: DataType) -> bool {
    !matches!(data_type, DataType::Float32 | DataType::Float64)
}

const fn is_float(data_type: DataType) -> bool {
    matches!(data_type, DataType::Float32 | DataType::Float64)
}

fn type_error(message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1403).with_message(message)
}

fn length_error(message: &'static str) -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message(message)
}

#[cfg(test)]
#[path = "numeric_into_tests.rs"]
mod tests;
