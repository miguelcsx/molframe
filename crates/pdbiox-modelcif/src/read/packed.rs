use crate::{ModelCifError, PackedIndices, PackedIntegers};
use std::mem::size_of;

pub(super) enum GrowingIndices {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
}

impl GrowingIndices {
    pub(super) const fn new() -> Self {
        Self::U8(Vec::new())
    }

    pub(super) fn len(&self) -> usize {
        match self {
            Self::U8(values) => values.len(),
            Self::U16(values) => values.len(),
            Self::U32(values) => values.len(),
        }
    }

    pub(super) fn live_bytes(&self) -> usize {
        match self {
            Self::U8(values) => values.capacity(),
            Self::U16(values) => values.capacity().saturating_mul(size_of::<u16>()),
            Self::U32(values) => values.capacity().saturating_mul(size_of::<u32>()),
        }
    }

    pub(super) fn push(&mut self, value: u32) -> Result<(), ModelCifError> {
        match self {
            Self::U8(values) if u8::try_from(value).is_ok() => push(
                values,
                u8::try_from(value).map_err(|_| ModelCifError::Capacity)?,
            ),
            Self::U16(values) if u16::try_from(value).is_ok() => push(
                values,
                u16::try_from(value).map_err(|_| ModelCifError::Capacity)?,
            ),
            Self::U32(values) => push(values, value),
            Self::U8(_) => {
                self.widen_to_u16()?;
                self.push(value)
            }
            Self::U16(_) => {
                self.widen_to_u32()?;
                self.push(value)
            }
        }
    }

    pub(super) fn finish(self) -> PackedIndices {
        match self {
            Self::U8(values) => PackedIndices::U8(values),
            Self::U16(values) => PackedIndices::U16(values),
            Self::U32(values) => PackedIndices::U32(values),
        }
    }

    fn widen_to_u16(&mut self) -> Result<(), ModelCifError> {
        let Self::U8(source) = std::mem::replace(self, Self::U8(Vec::new())) else {
            return Ok(());
        };
        let mut values = reserved(source.len())?;
        values.extend(source.into_iter().map(u16::from));
        *self = Self::U16(values);
        Ok(())
    }

    fn widen_to_u32(&mut self) -> Result<(), ModelCifError> {
        let Self::U16(source) = std::mem::replace(self, Self::U16(Vec::new())) else {
            return Ok(());
        };
        let mut values = reserved(source.len())?;
        values.extend(source.into_iter().map(u32::from));
        *self = Self::U32(values);
        Ok(())
    }
}

pub(super) enum GrowingIntegers {
    I8(Vec<i8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
}

impl GrowingIntegers {
    pub(super) const fn new() -> Self {
        Self::I8(Vec::new())
    }

    pub(super) fn len(&self) -> usize {
        match self {
            Self::I8(values) => values.len(),
            Self::I16(values) => values.len(),
            Self::I32(values) => values.len(),
            Self::I64(values) => values.len(),
        }
    }

    pub(super) fn get(&self, row: usize) -> Option<i64> {
        match self {
            Self::I8(values) => values.get(row).copied().map(i64::from),
            Self::I16(values) => values.get(row).copied().map(i64::from),
            Self::I32(values) => values.get(row).copied().map(i64::from),
            Self::I64(values) => values.get(row).copied(),
        }
    }

    pub(super) fn live_bytes(&self) -> usize {
        match self {
            Self::I8(values) => values.capacity(),
            Self::I16(values) => values.capacity().saturating_mul(size_of::<i16>()),
            Self::I32(values) => values.capacity().saturating_mul(size_of::<i32>()),
            Self::I64(values) => values.capacity().saturating_mul(size_of::<i64>()),
        }
    }

    pub(super) fn push(&mut self, value: i64) -> Result<(), ModelCifError> {
        match self {
            Self::I8(values) if i8::try_from(value).is_ok() => push(
                values,
                i8::try_from(value).map_err(|_| ModelCifError::Capacity)?,
            ),
            Self::I16(values) if i16::try_from(value).is_ok() => push(
                values,
                i16::try_from(value).map_err(|_| ModelCifError::Capacity)?,
            ),
            Self::I32(values) if i32::try_from(value).is_ok() => push(
                values,
                i32::try_from(value).map_err(|_| ModelCifError::Capacity)?,
            ),
            Self::I64(values) => push(values, value),
            Self::I8(_) | Self::I16(_) | Self::I32(_) => {
                self.widen_for(value)?;
                self.push(value)
            }
        }
    }

    pub(super) fn finish(self) -> PackedIntegers {
        match self {
            Self::I8(values) => PackedIntegers::I8(values),
            Self::I16(values) => PackedIntegers::I16(values),
            Self::I32(values) => PackedIntegers::I32(values),
            Self::I64(values) => PackedIntegers::I64(values),
        }
    }

    fn widen_for(&mut self, value: i64) -> Result<(), ModelCifError> {
        let target = width(value);
        loop {
            let current = integer_width(self);
            if current >= target {
                return Ok(());
            }
            match std::mem::replace(self, Self::I8(Vec::new())) {
                Self::I8(source) if target == 2 => {
                    *self = Self::I16(convert(source, i16::from)?);
                }
                Self::I8(source) if target == 4 => {
                    *self = Self::I32(convert(source, i32::from)?);
                }
                Self::I8(source) => {
                    *self = Self::I64(convert(source, i64::from)?);
                }
                Self::I16(source) if target == 4 => {
                    *self = Self::I32(convert(source, i32::from)?);
                }
                Self::I16(source) => {
                    *self = Self::I64(convert(source, i64::from)?);
                }
                Self::I32(source) => {
                    *self = Self::I64(convert(source, i64::from)?);
                }
                Self::I64(source) => {
                    *self = Self::I64(source);
                    return Ok(());
                }
            }
        }
    }
}

fn integer_width(values: &GrowingIntegers) -> u8 {
    match values {
        GrowingIntegers::I8(_) => 1,
        GrowingIntegers::I16(_) => 2,
        GrowingIntegers::I32(_) => 4,
        GrowingIntegers::I64(_) => 8,
    }
}

fn width(value: i64) -> u8 {
    if i8::try_from(value).is_ok() {
        1
    } else if i16::try_from(value).is_ok() {
        2
    } else if i32::try_from(value).is_ok() {
        4
    } else {
        8
    }
}

fn convert<T, U>(source: Vec<T>, widen: impl Fn(T) -> U) -> Result<Vec<U>, ModelCifError> {
    let mut values = reserved(source.len())?;
    values.extend(source.into_iter().map(widen));
    Ok(values)
}

pub(super) fn push<T>(values: &mut Vec<T>, value: T) -> Result<(), ModelCifError> {
    if values.len() == values.capacity() {
        let additional = values.capacity().max(64);
        values
            .try_reserve_exact(additional)
            .map_err(|_| ModelCifError::Allocation)?;
    }
    values.push(value);
    Ok(())
}

pub(super) fn reserved<T>(capacity: usize) -> Result<Vec<T>, ModelCifError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| ModelCifError::Allocation)?;
    Ok(values)
}
