//! Narrow row-oriented views over decoded `BinaryCIF` values.

use crate::codec::{Decoded, DecodedStringColumn};
use num_traits::ToPrimitive;
use std::sync::Arc;

pub(super) enum ColumnValues {
    Integers(IntegerValues),
    Floats(FloatValues),
    Strings(StringValues),
}

impl ColumnValues {
    pub(super) const fn from_f32(values: Vec<f32>) -> Self {
        Self::Floats(FloatValues::F32(values))
    }

    pub(super) fn new(decoded: Decoded) -> Self {
        match decoded {
            Decoded::Integers(values) => Self::Integers(IntegerValues::new(values)),
            Decoded::Floats(values) => Self::Floats(FloatValues::F64(values)),
            Decoded::Strings(values) => Self::Strings(StringValues::new(values)),
        }
    }

    pub(super) fn len(&self) -> usize {
        match self {
            Self::Integers(values) => values.len(),
            Self::Floats(values) => values.len(),
            Self::Strings(values) => values.len(),
        }
    }

    pub(super) fn value(&self, row: usize) -> Option<ValueRef<'_>> {
        match self {
            Self::Integers(values) => values.get(row).map(ValueRef::Integer),
            Self::Floats(values) => values.get(row).map(ValueRef::Float),
            Self::Strings(values) => values.get(row).map(ValueRef::Text),
        }
    }

    /// The string-dictionary position of a text row.
    pub(super) fn dictionary_slot(&self, row: usize) -> Option<u32> {
        match self {
            Self::Strings(values) => values.slot(row),
            _ => None,
        }
    }

    pub(super) fn compact_float(&mut self) {
        if let Self::Floats(values) = self {
            values.compact();
        }
    }
}

pub(super) enum IntegerValues {
    I8(Vec<i8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
}

impl IntegerValues {
    /// Narrows once: one pass finds the range, one pass builds the result.
    fn new(values: Vec<i64>) -> Self {
        if values.is_empty() {
            return Self::I8(Vec::new());
        }
        let (low, high) = values
            .iter()
            .fold((i64::MAX, i64::MIN), |(low, high), &value| {
                (low.min(value), high.max(value))
            });
        // Truncation below is lossless: every value lies within the checked range.
        if i8::try_from(low).is_ok() && i8::try_from(high).is_ok() {
            return Self::I8(values.iter().map(|&value| value as i8).collect());
        }
        if i16::try_from(low).is_ok() && i16::try_from(high).is_ok() {
            return Self::I16(values.iter().map(|&value| value as i16).collect());
        }
        if i32::try_from(low).is_ok() && i32::try_from(high).is_ok() {
            return Self::I32(values.iter().map(|&value| value as i32).collect());
        }
        Self::I64(values)
    }

    fn len(&self) -> usize {
        match self {
            Self::I8(values) => values.len(),
            Self::I16(values) => values.len(),
            Self::I32(values) => values.len(),
            Self::I64(values) => values.len(),
        }
    }

    fn get(&self, row: usize) -> Option<i64> {
        match self {
            Self::I8(values) => values.get(row).copied().map(i64::from),
            Self::I16(values) => values.get(row).copied().map(i64::from),
            Self::I32(values) => values.get(row).copied().map(i64::from),
            Self::I64(values) => values.get(row).copied(),
        }
    }
}

pub(super) enum FloatValues {
    F32(Vec<f32>),
    F64(Vec<f64>),
}

impl FloatValues {
    fn len(&self) -> usize {
        match self {
            Self::F32(values) => values.len(),
            Self::F64(values) => values.len(),
        }
    }

    fn get(&self, row: usize) -> Option<f64> {
        match self {
            Self::F32(values) => values.get(row).copied().map(f64::from),
            Self::F64(values) => values.get(row).copied(),
        }
    }

    fn compact(&mut self) {
        let Self::F64(values) = self else {
            return;
        };
        if values.iter().any(|value| value.to_f32().is_none()) {
            return;
        }
        let mut compact = Vec::with_capacity(values.len());
        for value in values.iter() {
            let Some(value) = value.to_f32() else {
                return;
            };
            compact.push(value);
        }
        *self = Self::F32(compact);
    }
}

pub(super) struct StringValues {
    dictionary: Arc<[Arc<str>]>,
    indices: StringIndices,
}

impl StringValues {
    fn new(values: DecodedStringColumn) -> Self {
        let (dictionary, indices) = values.into_parts();
        let indices = StringIndices::new(indices, dictionary.len());
        Self {
            dictionary,
            indices,
        }
    }

    fn len(&self) -> usize {
        self.indices.len()
    }

    fn get(&self, row: usize) -> Option<&str> {
        let index = self.indices.get(row)?;
        self.dictionary.get(index).map(AsRef::as_ref)
    }

    /// The dictionary position of a row's text, if the row names an entry.
    fn slot(&self, row: usize) -> Option<u32> {
        let index = self.indices.get(row)?;
        self.dictionary.get(index)?;
        u32::try_from(index).ok()
    }
}

enum StringIndices {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Arc<[u32]>),
}

impl StringIndices {
    /// Picks the width from the dictionary size, in one pass. The widest value
    /// of each narrow width stays out of range for a dictionary that fits, so
    /// an index past the dictionary clamps to a value `get` still rejects.
    fn new(indices: Arc<[u32]>, dictionary_len: usize) -> Self {
        if dictionary_len < usize::from(u8::MAX) {
            return Self::U8(
                indices
                    .iter()
                    .map(|&index| u8::try_from(index).unwrap_or(u8::MAX))
                    .collect(),
            );
        }
        if dictionary_len < usize::from(u16::MAX) {
            return Self::U16(
                indices
                    .iter()
                    .map(|&index| u16::try_from(index).unwrap_or(u16::MAX))
                    .collect(),
            );
        }
        Self::U32(indices)
    }

    fn len(&self) -> usize {
        match self {
            Self::U8(values) => values.len(),
            Self::U16(values) => values.len(),
            Self::U32(values) => values.len(),
        }
    }

    fn get(&self, row: usize) -> Option<usize> {
        match self {
            Self::U8(values) => values.get(row).copied().map(usize::from),
            Self::U16(values) => values.get(row).copied().map(usize::from),
            Self::U32(values) => values
                .get(row)
                .copied()
                .and_then(|value| usize::try_from(value).ok()),
        }
    }
}

#[derive(Clone, Copy)]
pub(super) enum ValueRef<'a> {
    Inapplicable,
    Unknown,
    Text(&'a str),
    Integer(i64),
    Float(f64),
}

#[cfg(test)]
#[path = "values_tests.rs"]
mod tests;
