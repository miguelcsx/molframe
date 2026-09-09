//! Stateful single-pass decoding of one indexed `BinaryCIF` column.

mod checkpoint;
mod metrics;

use super::super::index::IndexedColumn;
use super::decoded::{ColumnChunk, DecodedChunk};
use super::stream_primitive::{
    PayloadInput, data_width, packing_bounds, read_encoding, read_float, read_integer,
    validate_offsets,
};
use crate::{DataType, Decoded, Encoding};
use pdbiox_core::{Code, Diagnostic, SourceBytes};
use std::ops::Range;

pub(super) use super::stream_primitive::PAYLOAD_BUFFER_BYTES;
pub(super) use checkpoint::ColumnCheckpoints;

#[derive(Debug)]
pub(super) struct ColumnStream {
    pub(super) name: String,
    values: ValueStream,
    mask: Option<IntegerStream>,
}

impl ColumnStream {
    pub(super) fn new<S: SourceBytes>(
        source: &mut S,
        indexed: IndexedColumn,
        window_bytes: usize,
    ) -> Result<Self, Diagnostic> {
        let encoding = read_encoding(source, indexed.data.encoding, window_bytes)?;
        let values = ValueStream::new(&encoding, indexed.data.payload)?;
        let mask = match indexed.mask {
            Some(mask) => {
                let encoding = read_encoding(source, mask.encoding, window_bytes)?;
                Some(IntegerStream::new(&encoding, mask.payload)?)
            }
            None => None,
        };
        Ok(Self {
            name: indexed.name,
            values,
            mask,
        })
    }

    pub(super) fn empty_decoded(&self) -> DecodedChunk {
        let values = match &self.values {
            ValueStream::Integer(_) => ColumnChunk::Integer(Vec::new()),
            ValueStream::Float(_) => ColumnChunk::Float(Vec::new()),
            ValueStream::Text(_) => ColumnChunk::Text(Vec::new()),
        };
        DecodedChunk {
            values,
            mask: self.mask.as_ref().map(|_| Vec::new()),
        }
    }

    pub(super) fn decode_into<S: SourceBytes>(
        &mut self,
        source: &mut S,
        rows: usize,
        output: &mut DecodedChunk,
    ) -> Result<(), Diagnostic> {
        output.clear();
        self.values.decode_into(source, rows, &mut output.values)?;
        if let (Some(mask), Some(values)) = (&mut self.mask, &mut output.mask) {
            for _ in 0..rows {
                let value = mask.next(source)?.ok_or_else(short_column)?;
                let value = u8::try_from(value).map_err(|_| invalid_mask(value))?;
                if value > 2 {
                    return Err(invalid_mask(i64::from(value)));
                }
                values.push(value);
            }
        }
        Ok(())
    }

    pub(super) fn text(&self, index: u32) -> Option<&str> {
        match &self.values {
            ValueStream::Text(dictionary) => dictionary.text(index),
            ValueStream::Integer(_) | ValueStream::Float(_) => None,
        }
    }

    pub(super) fn retained_bytes(&self) -> usize {
        self.name
            .capacity()
            .saturating_add(self.values.retained_bytes())
            .saturating_add(self.mask.as_ref().map_or(0, IntegerStream::retained_bytes))
    }
}

#[derive(Debug)]
enum ValueStream {
    Integer(IntegerStream),
    Float(FloatStream),
    Text(TextStream),
}

impl ValueStream {
    fn new(encoding: &[Encoding], payload: Range<u64>) -> Result<Self, Diagnostic> {
        let Some(first) = encoding.first() else {
            return Err(codec_error("BinaryCIF column has no encoding"));
        };
        match first {
            Encoding::StringArray {
                data_encoding,
                string_data,
                offset_encoding,
                offsets,
            } => {
                let Decoded::Integers(decoded_offsets) =
                    crate::codec::decode_borrowed(offset_encoding, offsets)?
                else {
                    return Err(codec_error("BinaryCIF string offsets are not integers"));
                };
                let mut compact = Vec::with_capacity(decoded_offsets.len());
                for offset in decoded_offsets {
                    compact.push(
                        u32::try_from(offset)
                            .map_err(|_| codec_error("BinaryCIF string offset exceeds u32"))?,
                    );
                }
                validate_offsets(string_data, &compact)?;
                Ok(Self::Text(TextStream {
                    indices: IntegerStream::new(data_encoding, payload)?,
                    dictionary: string_data.clone(),
                    offsets: compact,
                }))
            }
            Encoding::FixedPoint { .. }
            | Encoding::IntervalQuantization { .. }
            | Encoding::ByteArray {
                r#type: DataType::Float32 | DataType::Float64,
            } => Ok(Self::Float(FloatStream::new(encoding, payload)?)),
            _ => Ok(Self::Integer(IntegerStream::new(encoding, payload)?)),
        }
    }

    fn decode_into<S: SourceBytes>(
        &mut self,
        source: &mut S,
        rows: usize,
        output: &mut ColumnChunk,
    ) -> Result<(), Diagnostic> {
        match (self, output) {
            (Self::Integer(stream), ColumnChunk::Integer(values)) => {
                for _ in 0..rows {
                    values.push(stream.next(source)?.ok_or_else(short_column)?);
                }
            }
            (Self::Float(stream), ColumnChunk::Float(values)) => {
                for _ in 0..rows {
                    let value = stream.next(source)?.ok_or_else(short_column)?;
                    values.push(num_traits::ToPrimitive::to_f32(&value).ok_or_else(|| {
                        codec_error("BinaryCIF float exceeds the f32 structure range")
                    })?);
                }
            }
            (Self::Text(stream), ColumnChunk::Text(values)) => {
                for _ in 0..rows {
                    let index = stream.indices.next(source)?.ok_or_else(short_column)?;
                    let index = if index == -1 {
                        u32::MAX
                    } else {
                        u32::try_from(index)
                            .map_err(|_| codec_error("BinaryCIF string index exceeds u32"))?
                    };
                    values.push(index);
                }
            }
            _ => return Err(codec_error("BinaryCIF decode buffer type changed")),
        }
        Ok(())
    }

    fn retained_bytes(&self) -> usize {
        match self {
            Self::Integer(stream) => stream.retained_bytes(),
            Self::Float(stream) => stream.retained_bytes(),
            Self::Text(stream) => stream.retained_bytes(),
        }
    }
}

#[derive(Debug)]
struct TextStream {
    indices: IntegerStream,
    dictionary: String,
    offsets: Vec<u32>,
}

impl TextStream {
    fn text(&self, index: u32) -> Option<&str> {
        if index == u32::MAX {
            return Some("");
        }
        let index = index as usize;
        let start = usize::try_from(*self.offsets.get(index)?).ok()?;
        let end = usize::try_from(*self.offsets.get(index.checked_add(1)?)?).ok()?;
        self.dictionary.get(start..end)
    }

    fn retained_bytes(&self) -> usize {
        self.indices
            .retained_bytes()
            .saturating_add(self.dictionary.capacity())
            .saturating_add(
                self.offsets
                    .capacity()
                    .saturating_mul(std::mem::size_of::<u32>()),
            )
    }
}

#[derive(Clone, Copy, Debug)]
enum Expansion {
    Delta(i64),
    RunLength { remaining: usize, value: i64 },
}

#[derive(Debug)]
struct IntegerStream {
    input: PayloadInput,
    data_type: DataType,
    packing: Option<(i64, i64)>,
    expansions: Vec<Expansion>,
}

impl IntegerStream {
    fn new(encoding: &[Encoding], payload: Range<u64>) -> Result<Self, Diagnostic> {
        let Some((last, prefix)) = encoding.split_last() else {
            return Err(codec_error("integer column has no ByteArray encoding"));
        };
        let Encoding::ByteArray { r#type } = last else {
            return Err(codec_error("integer column does not end in ByteArray"));
        };
        if matches!(r#type, DataType::Float32 | DataType::Float64) {
            return Err(codec_error("integer column uses floating-point bytes"));
        }
        let (transforms, packing) = match prefix.split_last() {
            Some((
                Encoding::IntegerPacking {
                    byte_count,
                    is_unsigned,
                    ..
                },
                transforms,
            )) => (transforms, Some(packing_bounds(*byte_count, *is_unsigned)?)),
            _ => (prefix, None),
        };
        let mut expansions = Vec::with_capacity(transforms.len());
        for transform in transforms {
            expansions.push(match transform {
                Encoding::Delta { origin, .. } => Expansion::Delta(*origin),
                Encoding::RunLength { .. } => Expansion::RunLength {
                    remaining: 0,
                    value: 0,
                },
                _ => {
                    return Err(codec_error(format!(
                        "unsupported streaming integer codec: {transform:?}"
                    )));
                }
            });
        }
        Ok(Self {
            input: PayloadInput::new(payload),
            data_type: *r#type,
            packing,
            expansions,
        })
    }

    fn next<S: SourceBytes>(&mut self, source: &mut S) -> Result<Option<i64>, Diagnostic> {
        decode_expansion(
            &mut self.input,
            self.data_type,
            self.packing,
            &mut self.expansions,
            source,
        )
    }

    fn retained_bytes(&self) -> usize {
        self.input.retained_bytes().saturating_add(
            self.expansions
                .capacity()
                .saturating_mul(std::mem::size_of::<Expansion>()),
        )
    }
}

fn decode_expansion<S: SourceBytes>(
    input: &mut PayloadInput,
    data_type: DataType,
    packing: Option<(i64, i64)>,
    expansions: &mut [Expansion],
    source: &mut S,
) -> Result<Option<i64>, Diagnostic> {
    let Some((expansion, inner)) = expansions.split_first_mut() else {
        return next_packed(input, data_type, packing, source);
    };
    match *expansion {
        Expansion::Delta(previous) => {
            let Some(delta) = decode_expansion(input, data_type, packing, inner, source)? else {
                return Ok(None);
            };
            let value = previous
                .checked_add(delta)
                .ok_or_else(|| codec_error("BinaryCIF delta overflow"))?;
            *expansion = Expansion::Delta(value);
            Ok(Some(value))
        }
        Expansion::RunLength { remaining, value } if remaining > 0 => {
            *expansion = Expansion::RunLength {
                remaining: remaining - 1,
                value,
            };
            Ok(Some(value))
        }
        Expansion::RunLength { .. } => {
            let Some(value) = decode_expansion(input, data_type, packing, inner, source)? else {
                return Ok(None);
            };
            let count = decode_expansion(input, data_type, packing, inner, source)?
                .ok_or_else(|| codec_error("BinaryCIF run has no count"))?;
            let count = usize::try_from(count)
                .map_err(|_| codec_error("BinaryCIF run count is invalid"))?;
            if count == 0 {
                return Err(codec_error("BinaryCIF run count is zero"));
            }
            *expansion = Expansion::RunLength {
                remaining: count - 1,
                value,
            };
            Ok(Some(value))
        }
    }
}

fn next_packed<S: SourceBytes>(
    input: &mut PayloadInput,
    data_type: DataType,
    packing: Option<(i64, i64)>,
    source: &mut S,
) -> Result<Option<i64>, Diagnostic> {
    let Some((upper, lower)) = packing else {
        return next_word(input, data_type, source);
    };
    let Some(mut value) = next_word(input, data_type, source)? else {
        return Ok(None);
    };
    let mut total = 0_i64;
    loop {
        total = total
            .checked_add(value)
            .ok_or_else(|| codec_error("BinaryCIF packed integer overflow"))?;
        if value != upper && value != lower {
            return Ok(Some(total));
        }
        value = next_word(input, data_type, source)?
            .ok_or_else(|| codec_error("unterminated BinaryCIF packed integer"))?;
    }
}

fn next_word<S: SourceBytes>(
    input: &mut PayloadInput,
    data_type: DataType,
    source: &mut S,
) -> Result<Option<i64>, Diagnostic> {
    let width = data_width(data_type);
    let Some(bytes) = input.read(source, width)? else {
        return Ok(None);
    };
    read_integer(bytes, data_type).map(Some)
}

#[derive(Clone, Copy, Debug)]
enum Conversion {
    Identity,
    FixedPoint(f64),
    Interval { minimum: f64, increment: f64 },
}

#[derive(Debug)]
enum FloatInner {
    Literal { input: PayloadInput, kind: DataType },
    Integer(IntegerStream),
}

#[derive(Debug)]
struct FloatStream {
    inner: FloatInner,
    conversion: Conversion,
}

impl FloatStream {
    fn new(encoding: &[Encoding], payload: Range<u64>) -> Result<Self, Diagnostic> {
        if let [Encoding::ByteArray { r#type }] = encoding
            && matches!(r#type, DataType::Float32 | DataType::Float64)
        {
            return Ok(Self {
                inner: FloatInner::Literal {
                    input: PayloadInput::new(payload),
                    kind: *r#type,
                },
                conversion: Conversion::Identity,
            });
        }
        let (conversion, integer_encoding) = match encoding.split_first() {
            Some((Encoding::FixedPoint { factor, .. }, rest)) => {
                (Conversion::FixedPoint(*factor), rest.to_vec())
            }
            Some((
                Encoding::IntervalQuantization {
                    min,
                    max,
                    num_steps,
                    ..
                },
                rest,
            )) if *num_steps >= 2 => (
                Conversion::Interval {
                    minimum: *min,
                    increment: (max - min) / f64::from(*num_steps - 1),
                },
                rest.to_vec(),
            ),
            _ => {
                return Err(codec_error(
                    "unsupported streaming floating-point codec chain",
                ));
            }
        };
        Ok(Self {
            inner: FloatInner::Integer(IntegerStream::new(&integer_encoding, payload)?),
            conversion,
        })
    }

    fn next<S: SourceBytes>(&mut self, source: &mut S) -> Result<Option<f64>, Diagnostic> {
        let value = match &mut self.inner {
            FloatInner::Literal { input, kind } => {
                let width = data_width(*kind);
                let Some(bytes) = input.read(source, width)? else {
                    return Ok(None);
                };
                read_float(bytes, *kind)?
            }
            FloatInner::Integer(stream) => {
                let Some(value) = stream.next(source)? else {
                    return Ok(None);
                };
                num_traits::ToPrimitive::to_f64(&value)
                    .ok_or_else(|| codec_error("BinaryCIF integer exceeds f64 range"))?
            }
        };
        let converted = match self.conversion {
            Conversion::Identity => value,
            Conversion::FixedPoint(factor) if factor != 0.0 => value / factor,
            Conversion::Interval { minimum, increment } => minimum + increment * value,
            Conversion::FixedPoint(_) => {
                return Err(codec_error("BinaryCIF fixed-point factor is zero"));
            }
        };
        Ok(Some(converted))
    }

    fn retained_bytes(&self) -> usize {
        match &self.inner {
            FloatInner::Literal { input, .. } => input.retained_bytes(),
            FloatInner::Integer(stream) => stream.retained_bytes(),
        }
    }
}

fn invalid_mask(value: i64) -> Diagnostic {
    codec_error(format!("invalid BinaryCIF mask value {value}"))
}

fn short_column() -> Diagnostic {
    codec_error("BinaryCIF column is shorter than its declared row count")
}

fn codec_error(message: impl Into<Box<str>>) -> Diagnostic {
    Diagnostic::new(Code::E1401).with_message(message)
}
