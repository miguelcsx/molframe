use num_traits::ToPrimitive;
use pdbiox_core::{Code, Diagnostic};
use serde_bytes::ByteBuf;

const HEADER_BYTES: usize = 12;

pub(super) fn decode_i32(bytes: &[u8]) -> Result<Vec<i32>, Diagnostic> {
    let header = Header::read(bytes)?;
    let payload = &bytes[HEADER_BYTES..];
    let decoded = match header.kind {
        2 => payload
            .iter()
            .map(|value| i32::from(value.cast_signed()))
            .collect(),
        3 => read_i16(payload)?.into_iter().map(i32::from).collect(),
        4 => read_i32(payload)?,
        6 | 7 | 16 => run_length(&read_i32(payload)?, header.length)?,
        8 => delta(run_length(&read_i32(payload)?, header.length)?)?,
        14 => unpack_i16(&read_i16(payload)?, header.length)?,
        15 => unpack_i8(payload, header.length)?,
        _ => return Err(codec_error("integer", header.kind)),
    };
    exact_len(decoded, header.length)
}

pub(super) fn decode_f32(bytes: &[u8]) -> Result<Vec<f32>, Diagnostic> {
    let header = Header::read(bytes)?;
    let payload = &bytes[HEADER_BYTES..];
    let decoded = match header.kind {
        1 => read_f32(payload)?,
        9 => divide(
            run_length(&read_i32(payload)?, header.length)?,
            header.parameter,
        )?,
        10 => divide(
            delta(unpack_i16(&read_i16(payload)?, header.length)?)?,
            header.parameter,
        )?,
        11 => divide(
            read_i16(payload)?.into_iter().map(i32::from).collect(),
            header.parameter,
        )?,
        12 => divide(
            unpack_i16(&read_i16(payload)?, header.length)?,
            header.parameter,
        )?,
        13 => divide(unpack_i8(payload, header.length)?, header.parameter)?,
        _ => return Err(codec_error("floating-point", header.kind)),
    };
    exact_len(decoded, header.length)
}

pub(super) fn decode_chars(bytes: &[u8]) -> Result<Vec<u8>, Diagnostic> {
    decode_i32(bytes)?
        .into_iter()
        .map(|value| {
            u8::try_from(value).map_err(|_| {
                Diagnostic::new(Code::E1102).with_message("MMTF character is outside byte range")
            })
        })
        .collect()
}

pub(super) fn decode_strings(bytes: &[u8]) -> Result<Vec<String>, Diagnostic> {
    let header = Header::read(bytes)?;
    if header.kind != 5 || header.parameter <= 0 {
        return Err(codec_error("fixed-width string", header.kind));
    }
    let width = usize::try_from(header.parameter).map_err(|_| codec_error("string", 5))?;
    let payload = &bytes[HEADER_BYTES..];
    let expected = width
        .checked_mul(header.length)
        .ok_or_else(|| length_error(header.length, payload.len()))?;
    if payload.len() != expected {
        return Err(length_error(header.length, payload.len() / width));
    }
    payload
        .chunks_exact(width)
        .map(|chunk| {
            let end = match chunk.iter().position(|byte| *byte == 0) {
                Some(end) => end,
                None => chunk.len(),
            };
            std::str::from_utf8(&chunk[..end])
                .map(str::to_owned)
                .map_err(|_| Diagnostic::new(Code::E1102).with_message("MMTF string is not UTF-8"))
        })
        .collect()
}

pub(super) fn encode_i32(values: &[i32]) -> Result<ByteBuf, Diagnostic> {
    let mut bytes = header(4, values.len(), 0)?;
    for value in values {
        bytes.extend_from_slice(&value.to_be_bytes());
    }
    Ok(bytes.into())
}

pub(super) fn encode_f32(values: &[f32]) -> Result<ByteBuf, Diagnostic> {
    let mut bytes = header(1, values.len(), 0)?;
    for value in values {
        bytes.extend_from_slice(&value.to_bits().to_be_bytes());
    }
    Ok(bytes.into())
}

pub(super) fn encode_chars(values: &[u8]) -> Result<ByteBuf, Diagnostic> {
    let mut bytes = header(2, values.len(), 0)?;
    bytes.extend(values);
    Ok(bytes.into())
}

pub(super) fn encode_strings(values: &[String], width: usize) -> Result<ByteBuf, Diagnostic> {
    let parameter = i32::try_from(width).map_err(|_| codec_error("string", 5))?;
    let mut bytes = header(5, values.len(), parameter)?;
    for value in values {
        if value.len() > width || !value.is_ascii() {
            return Err(Diagnostic::new(Code::E4102)
                .with_context("identifier", value.as_str())
                .with_context("maximum bytes", width.to_string()));
        }
        bytes.extend_from_slice(value.as_bytes());
        bytes.resize(bytes.len() + width - value.len(), 0);
    }
    Ok(bytes.into())
}

struct Header {
    kind: i32,
    length: usize,
    parameter: i32,
}

impl Header {
    fn read(bytes: &[u8]) -> Result<Self, Diagnostic> {
        let header = bytes.get(..HEADER_BYTES).ok_or_else(|| {
            Diagnostic::new(Code::E1102).with_message("truncated MMTF codec header")
        })?;
        let length = usize::try_from(be_i32(&header[4..8])?).map_err(|_| length_error(0, 0))?;
        Ok(Self {
            kind: be_i32(&header[..4])?,
            length,
            parameter: be_i32(&header[8..12])?,
        })
    }
}

fn header(kind: i32, length: usize, parameter: i32) -> Result<Vec<u8>, Diagnostic> {
    let length = i32::try_from(length).map_err(|_| codec_error("header length", 0))?;
    let mut bytes = Vec::with_capacity(HEADER_BYTES);
    bytes.extend_from_slice(&kind.to_be_bytes());
    bytes.extend_from_slice(&length.to_be_bytes());
    bytes.extend_from_slice(&parameter.to_be_bytes());
    Ok(bytes)
}

fn read_i16(bytes: &[u8]) -> Result<Vec<i16>, Diagnostic> {
    if !bytes.len().is_multiple_of(2) {
        return Err(Diagnostic::new(Code::E1102).with_message("odd MMTF i16 payload length"));
    }
    bytes.chunks_exact(2).map(be_i16).collect()
}

fn read_i32(bytes: &[u8]) -> Result<Vec<i32>, Diagnostic> {
    if !bytes.len().is_multiple_of(4) {
        return Err(Diagnostic::new(Code::E1102).with_message("unaligned MMTF i32 payload"));
    }
    bytes.chunks_exact(4).map(be_i32).collect()
}

fn read_f32(bytes: &[u8]) -> Result<Vec<f32>, Diagnostic> {
    read_i32(bytes).map(|values| {
        values
            .into_iter()
            .map(|value| f32::from_bits(value.cast_unsigned()))
            .collect()
    })
}

fn be_i16(bytes: &[u8]) -> Result<i16, Diagnostic> {
    let array: [u8; 2] = bytes.try_into().map_err(|_| codec_error("i16", 3))?;
    Ok(i16::from_be_bytes(array))
}

fn be_i32(bytes: &[u8]) -> Result<i32, Diagnostic> {
    let array: [u8; 4] = bytes.try_into().map_err(|_| codec_error("i32", 4))?;
    Ok(i32::from_be_bytes(array))
}

fn run_length(values: &[i32], expected: usize) -> Result<Vec<i32>, Diagnostic> {
    if !values.len().is_multiple_of(2) {
        return Err(Diagnostic::new(Code::E1102).with_message("odd MMTF run-length payload"));
    }
    let mut decoded = Vec::with_capacity(expected);
    for pair in values.chunks_exact(2) {
        let count = usize::try_from(pair[1]).map_err(|_| length_error(expected, decoded.len()))?;
        let expanded = decoded
            .len()
            .checked_add(count)
            .ok_or_else(|| length_error(expected, decoded.len()))?;
        if expanded > expected {
            return Err(length_error(expected, expanded));
        }
        decoded.resize(decoded.len() + count, pair[0]);
    }
    Ok(decoded)
}

fn delta(mut values: Vec<i32>) -> Result<Vec<i32>, Diagnostic> {
    for position in 1..values.len() {
        values[position] = values[position - 1]
            .checked_add(values[position])
            .ok_or_else(|| {
                Diagnostic::new(Code::E1102).with_message("MMTF delta integer overflow")
            })?;
    }
    Ok(values)
}

fn unpack_i16(values: &[i16], expected: usize) -> Result<Vec<i32>, Diagnostic> {
    unpack(
        values.iter().map(|value| i32::from(*value)),
        i32::from(i16::MIN),
        i32::from(i16::MAX),
        expected,
    )
}

fn unpack_i8(values: &[u8], expected: usize) -> Result<Vec<i32>, Diagnostic> {
    unpack(
        values.iter().map(|value| i32::from(value.cast_signed())),
        i32::from(i8::MIN),
        i32::from(i8::MAX),
        expected,
    )
}

fn unpack(
    values: impl IntoIterator<Item = i32>,
    minimum: i32,
    maximum: i32,
    expected: usize,
) -> Result<Vec<i32>, Diagnostic> {
    let mut decoded = Vec::with_capacity(expected);
    let mut accumulator = 0i32;
    for value in values {
        accumulator = accumulator.checked_add(value).ok_or_else(|| {
            Diagnostic::new(Code::E1102).with_message("MMTF packed integer overflow")
        })?;
        if value != minimum && value != maximum {
            decoded.push(accumulator);
            accumulator = 0;
        }
    }
    if accumulator != 0 || decoded.len() != expected {
        return Err(length_error(expected, decoded.len()));
    }
    Ok(decoded)
}

fn divide(values: Vec<i32>, divisor: i32) -> Result<Vec<f32>, Diagnostic> {
    if divisor == 0 {
        return Err(Diagnostic::new(Code::E1102).with_message("zero MMTF float divisor"));
    }
    let Some(divisor) = divisor.to_f32() else {
        return Err(codec_error("float divisor", 0));
    };
    values
        .into_iter()
        .map(|value| {
            value
                .to_f32()
                .map(|value| value / divisor)
                .ok_or_else(|| codec_error("scaled float", 0))
        })
        .collect()
}

fn exact_len<T>(values: Vec<T>, expected: usize) -> Result<Vec<T>, Diagnostic> {
    if values.len() == expected {
        Ok(values)
    } else {
        Err(length_error(expected, values.len()))
    }
}

fn codec_error(target: &str, kind: i32) -> Diagnostic {
    Diagnostic::new(Code::E1102)
        .with_message("unsupported or invalid MMTF codec")
        .with_context("target", target)
        .with_context("codec", kind.to_string())
}

fn length_error(expected: usize, actual: usize) -> Diagnostic {
    Diagnostic::new(Code::E1103)
        .with_message("MMTF decoded length disagrees with its header")
        .with_context("expected", expected.to_string())
        .with_context("actual", actual.to_string())
}

#[cfg(test)]
#[path = "codec_tests.rs"]
mod tests;
