use crate::ModelCifError;
use std::mem::size_of;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PackedIndices {
    U8(Vec<u8>),
    U16(Vec<u16>),
    U32(Vec<u32>),
}

impl PackedIndices {
    pub(crate) fn from_u32(values: Vec<u32>, maximum: u32) -> Result<Self, ModelCifError> {
        if u8::try_from(maximum).is_ok() {
            let mut packed = reserved(values.len())?;
            for value in values {
                packed.push(u8::try_from(value).map_err(|_| ModelCifError::Capacity)?);
            }
            return Ok(Self::U8(packed));
        }
        if u16::try_from(maximum).is_ok() {
            let mut packed = reserved(values.len())?;
            for value in values {
                packed.push(u16::try_from(value).map_err(|_| ModelCifError::Capacity)?);
            }
            return Ok(Self::U16(packed));
        }
        Ok(Self::U32(values))
    }

    pub(crate) fn get(&self, row: usize) -> Option<u32> {
        match self {
            Self::U8(values) => values.get(row).copied().map(u32::from),
            Self::U16(values) => values.get(row).copied().map(u32::from),
            Self::U32(values) => values.get(row).copied(),
        }
    }

    pub(crate) fn bytes(&self) -> usize {
        match self {
            Self::U8(values) => values.len(),
            Self::U16(values) => values.len().saturating_mul(size_of::<u16>()),
            Self::U32(values) => values.len().saturating_mul(size_of::<u32>()),
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum PackedIntegers {
    I8(Vec<i8>),
    I16(Vec<i16>),
    I32(Vec<i32>),
    I64(Vec<i64>),
}

impl PackedIntegers {
    pub(crate) fn from_i64(
        values: &[i64],
        minimum: i64,
        maximum: i64,
    ) -> Result<Self, ModelCifError> {
        if minimum >= i64::from(i8::MIN) && maximum <= i64::from(i8::MAX) {
            let mut packed = reserved(values.len())?;
            for &value in values {
                packed.push(i8::try_from(value).map_err(|_| ModelCifError::Capacity)?);
            }
            return Ok(Self::I8(packed));
        }
        if minimum >= i64::from(i16::MIN) && maximum <= i64::from(i16::MAX) {
            let mut packed = reserved(values.len())?;
            for &value in values {
                packed.push(i16::try_from(value).map_err(|_| ModelCifError::Capacity)?);
            }
            return Ok(Self::I16(packed));
        }
        if minimum >= i64::from(i32::MIN) && maximum <= i64::from(i32::MAX) {
            let mut packed = reserved(values.len())?;
            for &value in values {
                packed.push(i32::try_from(value).map_err(|_| ModelCifError::Capacity)?);
            }
            return Ok(Self::I32(packed));
        }
        let mut packed = reserved(values.len())?;
        packed.extend_from_slice(values);
        Ok(Self::I64(packed))
    }

    pub(crate) fn get(&self, row: usize) -> Option<i64> {
        match self {
            Self::I8(values) => values.get(row).copied().map(i64::from),
            Self::I16(values) => values.get(row).copied().map(i64::from),
            Self::I32(values) => values.get(row).copied().map(i64::from),
            Self::I64(values) => values.get(row).copied(),
        }
    }

    pub(crate) fn bytes(&self) -> usize {
        match self {
            Self::I8(values) => values.len(),
            Self::I16(values) => values.len().saturating_mul(size_of::<i16>()),
            Self::I32(values) => values.len().saturating_mul(size_of::<i32>()),
            Self::I64(values) => values.len().saturating_mul(size_of::<i64>()),
        }
    }
}

pub(crate) fn reserved<T>(capacity: usize) -> Result<Vec<T>, ModelCifError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(capacity)
        .map_err(|_| ModelCifError::Allocation)?;
    Ok(values)
}
