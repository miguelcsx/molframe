use molframe_core::{Code, Diagnostic};
use num_traits::ToPrimitive;
use serde_bytes::ByteBuf;

const HEADER_BYTES: usize = 12;

pub(super) fn decode_i32(bytes: &[u8]) -> Result<Vec<i32>, Diagnostic> {
    let header = Header::read(bytes)?;
    let payload = &bytes[HEADER_BYTES..];
    let mut decoded = Vec::with_capacity(header.length);
    match header.kind {
        2 => payload
            .iter()
            .for_each(|value| decoded.push(i32::from(value.cast_signed()))),
        3 => visit_i16(payload, |value| {
            decoded.push(value);
            Ok(())
        })?,
        4 => visit_i32(payload, |value| {
            decoded.push(value);
            Ok(())
        })?,
        6 | 7 | 16 => visit_run_length(payload, header.length, |value| {
            decoded.push(value);
            Ok(())
        })?,
        8 => {
            let mut previous = 0i32;
            visit_run_length(payload, header.length, |difference| {
                previous = previous.checked_add(difference).ok_or_else(|| {
                    Diagnostic::new(Code::E1102).with_message("MMTF delta integer overflow")
                })?;
                decoded.push(previous);
                Ok(())
            })?;
        }
        14 => visit_packed_i16(payload, header.length, |value| {
            decoded.push(value);
            Ok(())
        })?,
        15 => visit_packed_i8(payload, header.length, |value| {
            decoded.push(value);
            Ok(())
        })?,
        _ => return Err(codec_error("integer", header.kind)),
    }
    exact_count(header.length, decoded.len())?;
    Ok(decoded)
}

pub(super) fn decode_f32(bytes: &[u8]) -> Result<Vec<f32>, Diagnostic> {
    let header = Header::read(bytes)?;
    let payload = &bytes[HEADER_BYTES..];
    let mut decoded = Vec::with_capacity(header.length);
    match header.kind {
        1 => visit_i32(payload, |bits| {
            decoded.push(f32::from_bits(bits.cast_unsigned()));
            Ok(())
        })?,
        9 => {
            let divisor = float_divisor(header.parameter)?;
            visit_run_length(payload, header.length, |value| {
                decoded.push(scaled(value, divisor)?);
                Ok(())
            })?;
        }
        10 => {
            let divisor = float_divisor(header.parameter)?;
            let mut previous = 0i32;
            visit_packed_i16(payload, header.length, |difference| {
                previous = previous.checked_add(difference).ok_or_else(|| {
                    Diagnostic::new(Code::E1102).with_message("MMTF delta integer overflow")
                })?;
                decoded.push(scaled(previous, divisor)?);
                Ok(())
            })?;
        }
        11 => {
            let divisor = float_divisor(header.parameter)?;
            visit_i16(payload, |value| {
                decoded.push(scaled(value, divisor)?);
                Ok(())
            })?;
        }
        12 => {
            let divisor = float_divisor(header.parameter)?;
            visit_packed_i16(payload, header.length, |value| {
                decoded.push(scaled(value, divisor)?);
                Ok(())
            })?;
        }
        13 => {
            let divisor = float_divisor(header.parameter)?;
            visit_packed_i8(payload, header.length, |value| {
                decoded.push(scaled(value, divisor)?);
                Ok(())
            })?;
        }
        _ => return Err(codec_error("floating-point", header.kind)),
    }
    exact_count(header.length, decoded.len())?;
    Ok(decoded)
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

fn visit_i16(
    bytes: &[u8],
    mut visit: impl FnMut(i32) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    if !bytes.len().is_multiple_of(2) {
        return Err(Diagnostic::new(Code::E1102).with_message("odd MMTF i16 payload length"));
    }
    for bytes in bytes.as_chunks::<2>().0 {
        visit(i32::from(be_i16(bytes)?))?;
    }
    Ok(())
}

fn visit_i32(
    bytes: &[u8],
    mut visit: impl FnMut(i32) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    if !bytes.len().is_multiple_of(4) {
        return Err(Diagnostic::new(Code::E1102).with_message("unaligned MMTF i32 payload"));
    }
    for bytes in bytes.as_chunks::<4>().0 {
        visit(be_i32(bytes)?)?;
    }
    Ok(())
}

fn be_i16(bytes: &[u8]) -> Result<i16, Diagnostic> {
    let array: [u8; 2] = bytes.try_into().map_err(|_| codec_error("i16", 3))?;
    Ok(i16::from_be_bytes(array))
}

fn be_i32(bytes: &[u8]) -> Result<i32, Diagnostic> {
    let array: [u8; 4] = bytes.try_into().map_err(|_| codec_error("i32", 4))?;
    Ok(i32::from_be_bytes(array))
}

fn visit_run_length(
    bytes: &[u8],
    expected: usize,
    mut visit: impl FnMut(i32) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    if !bytes.len().is_multiple_of(4) {
        return Err(Diagnostic::new(Code::E1102).with_message("unaligned MMTF i32 payload"));
    }
    if !bytes.len().is_multiple_of(8) {
        return Err(Diagnostic::new(Code::E1102).with_message("odd MMTF run-length payload"));
    }
    let mut actual = 0usize;
    for pair in bytes.as_chunks::<8>().0 {
        let value = be_i32(&pair[..4])?;
        let count =
            usize::try_from(be_i32(&pair[4..])?).map_err(|_| length_error(expected, actual))?;
        let expanded = actual
            .checked_add(count)
            .ok_or_else(|| length_error(expected, actual))?;
        if expanded > expected {
            return Err(length_error(expected, expanded));
        }
        for _ in 0..count {
            visit(value)?;
        }
        actual = expanded;
    }
    exact_count(expected, actual)
}

fn visit_packed_i16(
    bytes: &[u8],
    expected: usize,
    visit: impl FnMut(i32) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    if !bytes.len().is_multiple_of(2) {
        return Err(Diagnostic::new(Code::E1102).with_message("odd MMTF i16 payload length"));
    }
    visit_packed(
        bytes
            .as_chunks::<2>()
            .0
            .iter()
            .map(|bytes| i32::from(i16::from_be_bytes([bytes[0], bytes[1]]))),
        i32::from(i16::MIN),
        i32::from(i16::MAX),
        expected,
        visit,
    )
}

fn visit_packed_i8(
    bytes: &[u8],
    expected: usize,
    visit: impl FnMut(i32) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    visit_packed(
        bytes.iter().map(|value| i32::from(value.cast_signed())),
        i32::from(i8::MIN),
        i32::from(i8::MAX),
        expected,
        visit,
    )
}

fn visit_packed(
    values: impl IntoIterator<Item = i32>,
    minimum: i32,
    maximum: i32,
    expected: usize,
    mut visit: impl FnMut(i32) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
    let mut actual = 0usize;
    let mut accumulator = 0i32;
    for value in values {
        accumulator = accumulator.checked_add(value).ok_or_else(|| {
            Diagnostic::new(Code::E1102).with_message("MMTF packed integer overflow")
        })?;
        if value != minimum && value != maximum {
            visit(accumulator)?;
            actual += 1;
            accumulator = 0;
        }
    }
    if accumulator != 0 {
        return Err(length_error(expected, actual));
    }
    exact_count(expected, actual)
}

fn float_divisor(divisor: i32) -> Result<f32, Diagnostic> {
    if divisor == 0 {
        return Err(Diagnostic::new(Code::E1102).with_message("zero MMTF float divisor"));
    }
    divisor
        .to_f32()
        .ok_or_else(|| codec_error("float divisor", 0))
}

fn scaled(value: i32, divisor: f32) -> Result<f32, Diagnostic> {
    value
        .to_f32()
        .map(|value| value / divisor)
        .ok_or_else(|| codec_error("scaled float", 0))
}

fn exact_count(expected: usize, actual: usize) -> Result<(), Diagnostic> {
    if actual == expected {
        Ok(())
    } else {
        Err(length_error(expected, actual))
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
