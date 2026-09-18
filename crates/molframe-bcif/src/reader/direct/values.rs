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
    fn new(values: Vec<i64>) -> Self {
        let i8_values: Result<Vec<_>, _> = values.iter().copied().map(i8::try_from).collect();
        if let Ok(compact) = i8_values {
            return Self::I8(compact);
        }
        let i16_values: Result<Vec<_>, _> = values.iter().copied().map(i16::try_from).collect();
        if let Ok(compact) = i16_values {
            return Self::I16(compact);
        }
        let i32_values: Result<Vec<_>, _> = values.iter().copied().map(i32::try_from).collect();
        if let Ok(compact) = i32_values {
            return Self::I32(compact);
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
        Self {
            dictionary,
            indices: StringIndices::new(indices),
        }
    }

    fn len(&self) -> usize {
        self.indices.len()
    }

    fn get(&self, row: usize) -> Option<&str> {
        let index = self.indices.get(row)?;
        self.dictionary.get(index).map(AsRef::as_ref)
    }
}

enum StringIndices {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Arc<[u32]>),
}

impl StringIndices {
    fn new(indices: Arc<[u32]>) -> Self {
        let u8_indices: Result<Vec<_>, _> = indices.iter().copied().map(u8::try_from).collect();
        if let Ok(compact) = u8_indices {
            return Self::U8(compact);
        }
        let u16_indices: Result<Vec<_>, _> = indices.iter().copied().map(u16::try_from).collect();
        if let Ok(compact) = u16_indices {
            return Self::U16(compact);
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
