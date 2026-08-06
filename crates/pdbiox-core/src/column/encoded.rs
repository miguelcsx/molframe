//! Per-column physical encodings.
//!
//! Each column is stored in whichever form suits the data it holds. Identifier
//! and metadata columns repeat heavily — a chunk's model number is one value
//! repeated four thousand times, and its alternate-location column is usually
//! blank throughout — so storing them literally wastes memory that the
//! geometric kernels would otherwise have as cache.
//!
//! The rule that governs the choice is not "compress everything":
//!
//! > **Columns a kernel reads on every iteration stay plain.** Encoding pays for
//! > identifiers, flags and metadata. It does not pay for positions, where a
//! > branch in the read path costs more than the bandwidth it saves and defeats
//! > vectorisation outright.

use super::bits::{bit_width, pack, unpack_one};
use crate::symbol::SymbolId;

/// A value that can live in an encoded column.
///
/// Encodings work on the bit pattern, so anything with a lossless 64-bit
/// representation qualifies. Bit-packing a float is legal and never chosen: the
/// encoder's rules only reach for it where the values are small integers.
pub trait ColumnValue: Copy + PartialEq {
    /// Whether the bit pattern is an integer, so that narrowing it or taking
    /// differences between successive values means something.
    ///
    /// False for floating point. Bit-packing a float would sometimes succeed —
    /// small floats have leading zero bits — and would produce a column whose
    /// size depended on the exponents in it, which is not an encoding anyone
    /// chose.
    const IS_INTEGER: bool;

    /// The value as bits.
    fn to_bits(self) -> u64;
    /// The value from bits.
    fn from_bits(bits: u64) -> Self;
}

macro_rules! column_value_via_cast {
    ($($type:ty),*) => {
        $(impl ColumnValue for $type {
            const IS_INTEGER: bool = true;

            #[allow(clippy::cast_lossless, reason = "one body over every width")]
            fn to_bits(self) -> u64 { self as u64 }

            #[allow(clippy::cast_possible_wrap, reason = "reinterpretation, not conversion")]
            fn from_bits(bits: u64) -> Self { bits as Self }
        })*
    };
}

column_value_via_cast!(u8, u16, u32, u64, i8, i16, i32, i64);

impl ColumnValue for f32 {
    const IS_INTEGER: bool = false;

    fn to_bits(self) -> u64 {
        u64::from(self.to_bits())
    }
    fn from_bits(bits: u64) -> Self {
        Self::from_bits(bits as u32)
    }
}

impl ColumnValue for SymbolId {
    const IS_INTEGER: bool = true;

    fn to_bits(self) -> u64 {
        u64::from(self.get())
    }
    fn from_bits(bits: u64) -> Self {
        Self::from_raw(bits as u32)
    }
}

/// One column of a chunk, in whichever physical form fits its data.
///
/// # Examples
///
/// ```
/// use pdbiox_core::column::EncodedColumn;
///
/// let repeated = EncodedColumn::encode(&[7u32; 4096]);
/// assert!(repeated.is_constant());
/// assert_eq!(repeated.get(4095), Some(7));
///
/// let varied = EncodedColumn::encode(&[1u32, 9, 3, 9]);
/// assert_eq!(varied.get(2), Some(3));
/// assert_eq!(varied.len(), 4);
/// ```
#[derive(Clone, PartialEq, Debug)]
pub enum EncodedColumn<T> {
    /// Values as they are. The only form a kernel reads directly.
    Plain(Vec<T>),
    /// One value for the whole column.
    Constant {
        /// The repeated value.
        value: T,
        /// How many times it repeats.
        len: u32,
    },
    /// Values with the position each run ends at.
    RunLength {
        /// One value per run.
        values: Vec<T>,
        /// Position one past the last member of each run.
        run_ends: Vec<u32>,
    },
    /// Narrow integers at a fixed bit width.
    BitPacked {
        /// The packed bits.
        data: Vec<u8>,
        /// Bits per value.
        width: u8,
        /// How many values are packed.
        len: u32,
    },
    /// A first value and the difference to each successive one.
    Delta {
        /// The value at position zero.
        first: T,
        /// Differences between successive values.
        deltas: Vec<i32>,
    },
}

impl<T: ColumnValue> EncodedColumn<T> {
    /// The number of values held.
    #[must_use]
    pub fn len(&self) -> u32 {
        match self {
            Self::Plain(values) => values.len() as u32,
            Self::Constant { len, .. } | Self::BitPacked { len, .. } => *len,
            Self::RunLength { run_ends, .. } => match run_ends.last() {
                Some(end) => *end,
                None => 0,
            },
            Self::Delta { deltas, .. } => deltas.len() as u32 + 1,
        }
    }

    /// Returns true when the column holds no values.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns true when the column is one repeated value.
    #[must_use]
    pub const fn is_constant(&self) -> bool {
        matches!(self, Self::Constant { .. })
    }

    /// The values as a slice, when they are stored literally.
    ///
    /// A kernel that needs a slice asks for one and takes the other path when
    /// it is refused; nothing forces a column into plain form behind the
    /// caller's back.
    #[must_use]
    pub fn as_slice(&self) -> Option<&[T]> {
        match self {
            Self::Plain(values) => Some(values),
            _ => None,
        }
    }

    /// The value at `position`, or `None` if it is past the end.
    #[must_use]
    #[allow(
        clippy::cast_possible_wrap,
        reason = "delta arithmetic runs on bit patterns"
    )]
    pub fn get(&self, position: u32) -> Option<T> {
        if position >= self.len() {
            return None;
        }
        match self {
            Self::Plain(values) => values.get(position as usize).copied(),
            Self::Constant { value, .. } => Some(*value),
            Self::RunLength { values, run_ends } => {
                let run = run_ends.partition_point(|end| *end <= position);
                values.get(run).copied()
            }
            Self::BitPacked { data, width, .. } => {
                unpack_one(data, *width, position).map(T::from_bits)
            }
            Self::Delta { first, deltas } => {
                let mut value = first.to_bits() as i64;
                for delta in deltas.get(..position as usize)? {
                    value += i64::from(*delta);
                }
                Some(T::from_bits(value as u64))
            }
        }
    }

    /// Every value in order.
    ///
    /// Reading sequentially through a delta or run-length column this way is
    /// linear overall, where repeated random access would not be.
    pub fn iter(&self) -> ColumnIter<'_, T> {
        ColumnIter {
            column: self,
            position: 0,
            running: 0,
            run: 0,
        }
    }

    /// Chooses an encoding for `values` by what they contain.
    ///
    /// The rules, in order: one repeated value becomes constant; few runs
    /// become run-length; small non-negative integers become bit-packed when
    /// that is genuinely narrower; a near-monotonic integer sequence becomes
    /// deltas; anything else stays plain.
    #[must_use]
    pub fn encode(values: &[T]) -> Self {
        let Some(first) = values.first() else {
            return Self::Plain(Vec::new());
        };
        if values.iter().all(|value| value == first) {
            return Self::Constant {
                value: *first,
                len: values.len() as u32,
            };
        }
        if let Some(column) = Self::try_run_length(values) {
            return column;
        }
        if let Some(column) = Self::try_bit_packed(values) {
            return column;
        }
        if let Some(column) = Self::try_delta(values) {
            return column;
        }
        Self::Plain(values.to_vec())
    }

    /// Stores values literally, whatever they are.
    ///
    /// This is what a hot column asks for: it declines encoding rather than
    /// letting the rules decide.
    #[must_use]
    pub fn plain(values: &[T]) -> Self {
        Self::Plain(values.to_vec())
    }

    fn try_run_length(values: &[T]) -> Option<Self> {
        let mut runs = 1usize;
        for pair in values.windows(2) {
            if let [earlier, later] = pair
                && earlier != later
            {
                runs += 1;
            }
        }
        // Each run costs a value and an end position, so it must replace more
        // than two values' worth of storage to be worth the indirection.
        if runs * 4 > values.len() {
            return None;
        }
        let mut encoded = Self::RunLength {
            values: Vec::with_capacity(runs),
            run_ends: Vec::with_capacity(runs),
        };
        if let Self::RunLength {
            values: heads,
            run_ends,
        } = &mut encoded
        {
            for (position, value) in values.iter().enumerate() {
                if heads.last() == Some(value) {
                    if let Some(end) = run_ends.last_mut() {
                        *end = position as u32 + 1;
                    }
                } else {
                    heads.push(*value);
                    run_ends.push(position as u32 + 1);
                }
            }
        }
        Some(encoded)
    }

    fn try_bit_packed(values: &[T]) -> Option<Self> {
        if !T::IS_INTEGER {
            return None;
        }
        let max = values.iter().map(|value| value.to_bits()).max()?;
        let width = bit_width(max);
        if usize::from(width) >= size_of::<T>() * 8 {
            return None;
        }
        let raw: Vec<u64> = values.iter().map(|value| value.to_bits()).collect();
        Some(Self::BitPacked {
            data: pack(&raw, width),
            width,
            len: values.len() as u32,
        })
    }

    #[allow(
        clippy::cast_possible_wrap,
        reason = "delta arithmetic runs on bit patterns"
    )]
    fn try_delta(values: &[T]) -> Option<Self> {
        if !T::IS_INTEGER {
            return None;
        }
        let first = *values.first()?;
        let mut deltas = Vec::with_capacity(values.len().saturating_sub(1));
        for pair in values.windows(2) {
            let [earlier, later] = pair else { return None };
            let step = (later.to_bits() as i64).checked_sub(earlier.to_bits() as i64)?;
            deltas.push(i32::try_from(step).ok()?);
        }
        Some(Self::Delta { first, deltas })
    }
}

/// Sequential reader over an encoded column.
///
/// Carries the cursors that make a forward walk linear: the running total for a
/// delta column and the current run for a run-length one. Reading those forms by
/// repeated random access would cost a scan or a search per value.
#[derive(Debug)]
pub struct ColumnIter<'a, T> {
    column: &'a EncodedColumn<T>,
    position: u32,
    running: i64,
    run: usize,
}

impl<T: ColumnValue> Iterator for ColumnIter<'_, T> {
    type Item = T;

    #[allow(
        clippy::cast_possible_wrap,
        reason = "delta arithmetic runs on bit patterns"
    )]
    fn next(&mut self) -> Option<T> {
        if self.position >= self.column.len() {
            return None;
        }
        let value = match self.column {
            EncodedColumn::Delta { first, deltas } => {
                if self.position == 0 {
                    self.running = first.to_bits() as i64;
                } else {
                    self.running += i64::from(*deltas.get(self.position as usize - 1)?);
                }
                T::from_bits(self.running as u64)
            }
            EncodedColumn::RunLength { values, run_ends } => {
                while run_ends
                    .get(self.run)
                    .is_some_and(|end| *end <= self.position)
                {
                    self.run += 1;
                }
                *values.get(self.run)?
            }
            other => other.get(self.position)?,
        };
        self.position += 1;
        Some(value)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.column.len().saturating_sub(self.position) as usize;
        (remaining, Some(remaining))
    }
}

impl<T: ColumnValue> ExactSizeIterator for ColumnIter<'_, T> {}

impl<'a, T: ColumnValue> IntoIterator for &'a EncodedColumn<T> {
    type Item = T;
    type IntoIter = ColumnIter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

#[cfg(test)]
#[path = "encoded_tests.rs"]
mod tests;
