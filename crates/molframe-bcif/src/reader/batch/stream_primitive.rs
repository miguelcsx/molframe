//! Reusable payload windows and primitive `BinaryCIF` word decoding.

use crate::{DataType, Encoding};
use molframe_core::{Code, Diagnostic, SourceBytes};
use std::ops::Range;

pub(super) const PAYLOAD_BUFFER_BYTES: usize = 16 * 1024;

#[derive(Debug)]
pub(super) struct PayloadInput {
    next: u64,
    end: u64,
    buffer: Vec<u8>,
    position: usize,
}

impl PayloadInput {
    pub(super) fn new(range: Range<u64>) -> Self {
        Self {
            next: range.start,
            end: range.end,
            buffer: Vec::with_capacity(PAYLOAD_BUFFER_BYTES),
            position: 0,
        }
    }

    pub(super) fn read<S: SourceBytes>(
        &mut self,
        source: &mut S,
        width: usize,
    ) -> Result<Option<&[u8]>, Diagnostic> {
        if self.buffer.len().saturating_sub(self.position) < width {
            self.refill(source)?;
        }
        if self.buffer.len().saturating_sub(self.position) < width {
            return Ok(None);
        }
        let start = self.position;
        self.position = self.position.saturating_add(width);
        Ok(self.buffer.get(start..self.position))
    }

    pub(super) fn retained_bytes(&self) -> usize {
        self.buffer.capacity()
    }

    pub(super) fn checkpoint(&self) -> u64 {
        let unread = self.buffer.len().saturating_sub(self.position);
        let unread = match u64::try_from(unread) {
            Ok(unread) => unread,
            Err(_) => u64::MAX,
        };
        self.next.saturating_sub(unread)
    }

    pub(super) fn restore(&mut self, checkpoint: u64) {
        self.next = checkpoint;
        self.buffer.clear();
        self.position = 0;
    }

    fn refill<S: SourceBytes>(&mut self, source: &mut S) -> Result<(), Diagnostic> {
        let remaining = &self.buffer[self.position..];
        let mut carry = [0_u8; 8];
        carry[..remaining.len()].copy_from_slice(remaining);
        let carry_len = remaining.len();
        self.buffer.clear();
        self.buffer.extend_from_slice(&carry[..carry_len]);
        if self.next >= self.end {
            self.position = 0;
            return Ok(());
        }
        let remaining = match usize::try_from(self.end - self.next) {
            Ok(length) => length,
            Err(_) => usize::MAX,
        };
        let demand = PAYLOAD_BUFFER_BYTES
            .saturating_sub(carry_len)
            .min(remaining);
        let window = source.window(self.next, demand)?;
        if window.bytes().is_empty() {
            return Err(short_column());
        }
        self.buffer.extend_from_slice(window.bytes());
        self.next = self
            .next
            .checked_add(u64::try_from(window.bytes().len()).map_err(|_| short_column())?)
            .ok_or_else(short_column)?;
        self.position = 0;
        Ok(())
    }
}

pub(super) fn read_encoding<S: SourceBytes>(
    source: &mut S,
    range: Range<u64>,
    window_bytes: usize,
) -> Result<Vec<Encoding>, Diagnostic> {
    let length = usize::try_from(range.end.saturating_sub(range.start))
        .map_err(|_| error("BinaryCIF encoding metadata exceeds address space"))?;
    let mut bytes = Vec::with_capacity(length);
    let mut offset = range.start;
    while offset < range.end {
        let remaining = match usize::try_from(range.end - offset) {
            Ok(length) => length,
            Err(_) => usize::MAX,
        };
        let window = source.window(offset, remaining.min(window_bytes))?;
        if window.bytes().is_empty() {
            return Err(error("truncated BinaryCIF encoding metadata"));
        }
        bytes.extend_from_slice(window.bytes());
        offset = offset
            .checked_add(u64::try_from(window.bytes().len()).map_err(|_| short_column())?)
            .ok_or_else(short_column)?;
    }
    rmp_serde::from_slice(&bytes).map_err(|failure| error(failure.to_string()))
}

pub(super) fn validate_offsets(dictionary: &str, offsets: &[u32]) -> Result<(), Diagnostic> {
    let Some(first) = offsets.first() else {
        return Err(error("BinaryCIF string dictionary has no offsets"));
    };
    if *first != 0 {
        return Err(error("BinaryCIF string offsets do not start at zero"));
    }
    let mut previous = 0_u32;
    for offset in offsets {
        if *offset < previous {
            return Err(error("BinaryCIF string offsets are not ordered"));
        }
        previous = *offset;
    }
    if usize::try_from(previous).ok() != Some(dictionary.len()) {
        return Err(error("BinaryCIF string offsets do not span dictionary"));
    }
    Ok(())
}

pub(super) fn packing_bounds(bytes: u8, unsigned: bool) -> Result<(i64, i64), Diagnostic> {
    let upper = match (bytes, unsigned) {
        (1, true) => i64::from(u8::MAX),
        (1, false) => i64::from(i8::MAX),
        (2, true) => i64::from(u16::MAX),
        (2, false) => i64::from(i16::MAX),
        _ => return Err(error("BinaryCIF integer packing width is invalid")),
    };
    Ok((upper, -upper - 1))
}

pub(super) const fn data_width(data_type: DataType) -> usize {
    match data_type {
        DataType::Int8 | DataType::Uint8 => 1,
        DataType::Int16 | DataType::Uint16 => 2,
        DataType::Int32 | DataType::Uint32 | DataType::Float32 => 4,
        DataType::Float64 => 8,
    }
}

pub(super) fn read_integer(bytes: &[u8], data_type: DataType) -> Result<i64, Diagnostic> {
    let value = match data_type {
        DataType::Int8 => i64::from(i8::from_le_bytes([bytes[0]])),
        DataType::Uint8 => i64::from(bytes[0]),
        DataType::Int16 => i64::from(i16::from_le_bytes([bytes[0], bytes[1]])),
        DataType::Uint16 => i64::from(u16::from_le_bytes([bytes[0], bytes[1]])),
        DataType::Int32 => i64::from(i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
        DataType::Uint32 => i64::from(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
        DataType::Float32 | DataType::Float64 => {
            return Err(error("floating point word requested as integer"));
        }
    };
    Ok(value)
}

pub(super) fn read_float(bytes: &[u8], data_type: DataType) -> Result<f64, Diagnostic> {
    match data_type {
        DataType::Float32 => Ok(f64::from(f32::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3],
        ]))),
        DataType::Float64 => Ok(f64::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        ])),
        _ => Err(error("integer word requested as floating point")),
    }
}

fn short_column() -> Diagnostic {
    error("BinaryCIF column is shorter than its declared row count")
}

fn error(message: impl Into<Box<str>>) -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message(message)
}
