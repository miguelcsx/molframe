//! Checked XDR and post-2020 in-memory serializer cursor.

use super::error::TprError;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Precision {
    Single,
    Double,
}

#[derive(Debug)]
pub(super) struct Decoder<'a> {
    bytes: &'a [u8],
    offset: usize,
    precision: Precision,
    modern: bool,
}

impl<'a> Decoder<'a> {
    pub(super) const fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            offset: 0,
            precision: Precision::Single,
            modern: false,
        }
    }

    pub(super) const fn offset(&self) -> usize {
        self.offset
    }
    pub(super) const fn remaining(&self) -> usize {
        self.bytes.len() - self.offset
    }
    pub(super) fn set_modern(&mut self) {
        self.modern = true;
    }
    pub(super) fn set_precision(&mut self, width: i32) -> Result<(), TprError> {
        self.precision = match width {
            4 => Precision::Single,
            8 => Precision::Double,
            _ => return Err(TprError::UnsupportedPrecision(width)),
        };
        Ok(())
    }

    pub(super) fn i32(&mut self) -> Result<i32, TprError> {
        Ok(i32::from_be_bytes(self.array()?))
    }
    pub(super) fn u32(&mut self) -> Result<u32, TprError> {
        Ok(u32::from_be_bytes(self.array()?))
    }
    pub(super) fn i64(&mut self) -> Result<i64, TprError> {
        Ok(i64::from_be_bytes(self.array()?))
    }
    pub(super) fn u64(&mut self) -> Result<u64, TprError> {
        Ok(u64::from_be_bytes(self.array()?))
    }
    pub(super) fn real(&mut self) -> Result<f64, TprError> {
        match self.precision {
            Precision::Single => Ok(f64::from(f32::from_be_bytes(self.array()?))),
            Precision::Double => Ok(f64::from_be_bytes(self.array()?)),
        }
    }
    pub(super) fn double(&mut self) -> Result<f64, TprError> {
        Ok(f64::from_be_bytes(self.array()?))
    }
    pub(super) fn ushort(&mut self) -> Result<u16, TprError> {
        if self.modern {
            Ok(u16::from_be_bytes(self.array()?))
        } else {
            u16::try_from(self.u32()?).map_err(|_| TprError::SizeOverflow)
        }
    }
    pub(super) fn uchar(&mut self) -> Result<u8, TprError> {
        if self.modern {
            Ok(self.take(1)?[0])
        } else {
            u8::try_from(self.u32()?).map_err(|_| TprError::SizeOverflow)
        }
    }
    pub(super) fn string(&mut self, field: &'static str) -> Result<String, TprError> {
        let bytes = if self.modern {
            let length = usize::try_from(self.u64()?).map_err(|_| TprError::SizeOverflow)?;
            self.take(length)?
        } else {
            let _declared = self.u32()?;
            let length = self.count_u32(field)?;
            let bytes = self.take(length)?;
            let padding = (4 - length % 4) % 4;
            self.take(padding)?;
            bytes
        };
        std::str::from_utf8(bytes)
            .map(str::to_owned)
            .map_err(|_| TprError::InvalidText { field })
    }
    pub(super) fn count(&mut self, field: &'static str) -> Result<usize, TprError> {
        let offset = self.offset;
        let value = self.i32()?;
        if value < 0 {
            return Err(TprError::InvalidCount {
                field,
                value: i64::from(value),
                offset,
            });
        }
        usize::try_from(value).map_err(|_| TprError::SizeOverflow)
    }
    pub(super) fn skip_reals(&mut self, count: usize) -> Result<(), TprError> {
        let width = match self.precision {
            Precision::Single => 4,
            Precision::Double => 8,
        };
        let bytes = count.checked_mul(width).ok_or(TprError::SizeOverflow)?;
        self.take(bytes).map(|_| ())
    }

    fn count_u32(&mut self, _field: &'static str) -> Result<usize, TprError> {
        usize::try_from(self.u32()?).map_err(|_| TprError::SizeOverflow)
    }
    fn take(&mut self, length: usize) -> Result<&'a [u8], TprError> {
        let end = self
            .offset
            .checked_add(length)
            .ok_or(TprError::SizeOverflow)?;
        let bytes = self
            .bytes
            .get(self.offset..end)
            .ok_or(TprError::Truncated {
                offset: self.offset,
            })?;
        self.offset = end;
        Ok(bytes)
    }
    fn array<const N: usize>(&mut self) -> Result<[u8; N], TprError> {
        let bytes = self.take(N)?;
        let mut array = [0; N];
        array.copy_from_slice(bytes);
        Ok(array)
    }
}
