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

use super::super::bits::{bit_width, pack_mapped, unpack_one};
use super::iter::ColumnIter;
use crate::symbol::SymbolId;
use std::mem::size_of;

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
    fn from_bits(bits: u64) -> Option<Self>;
}

macro_rules! unsigned_column_value {
    ($($type:ty),*) => {
        $(impl ColumnValue for $type {
            const IS_INTEGER: bool = true;

            fn to_bits(self) -> u64 {
                u64::from(self)
            }

            fn from_bits(bits: u64) -> Option<Self> {
                Self::try_from(bits).ok()
            }
        })*
    };
}

macro_rules! signed_column_value {
    ($($type:ty),*) => {
        $(impl ColumnValue for $type {
            const IS_INTEGER: bool = true;

            fn to_bits(self) -> u64 {
                i64::from(self).cast_unsigned()
            }

            fn from_bits(bits: u64) -> Option<Self> {
                Self::try_from(bits.cast_signed()).ok()
            }
        })*
    };
}

unsigned_column_value!(u8, u16, u32);
signed_column_value!(i8, i16, i32);

impl ColumnValue for u64 {
    const IS_INTEGER: bool = true;

    fn to_bits(self) -> u64 {
        self
    }

    fn from_bits(bits: u64) -> Option<Self> {
        Some(bits)
    }
}

impl ColumnValue for i64 {
    const IS_INTEGER: bool = true;

    fn to_bits(self) -> u64 {
        self.cast_unsigned()
    }

    fn from_bits(bits: u64) -> Option<Self> {
        Some(bits.cast_signed())
    }
}

impl ColumnValue for f32 {
    const IS_INTEGER: bool = false;

    fn to_bits(self) -> u64 {
        u64::from(self.to_bits())
    }

    fn from_bits(bits: u64) -> Option<Self> {
        u32::try_from(bits).ok().map(Self::from_bits)
    }
}

impl ColumnValue for SymbolId {
    const IS_INTEGER: bool = true;

    fn to_bits(self) -> u64 {
        u64::from(self.get())
    }

    fn from_bits(bits: u64) -> Option<Self> {
        u32::try_from(bits).ok().map(Self::from_raw)
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
        /// The absolute value at every `DELTA_CHECKPOINT_STRIDE`-th position.
        ///
        /// Without these, reading position `p` sums `p` deltas, so walking a
        /// column by index costs `O(n^2)`. Entry `k` holds the value at
        /// position `(k + 1) * DELTA_CHECKPOINT_STRIDE`, which bounds any read
        /// to one stride of additions. At four bytes per delta and eight per
        /// checkpoint the table adds about three percent to the encoding.
        checkpoints: Vec<i64>,
    },
}

impl<T: ColumnValue> EncodedColumn<T> {
    /// The number of values held.
    ///
    #[must_use]
    pub fn len(&self) -> usize {
        match self {
            Self::Plain(values) => values.len(),
            Self::Constant { len, .. } | Self::BitPacked { len, .. } => *len as usize,
            Self::RunLength { run_ends, .. } => match run_ends.last() {
                Some(end) => *end as usize,
                None => 0,
            },
            Self::Delta { deltas, .. } => deltas.len() + 1,
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
    pub fn get(&self, position: u32) -> Option<T> {
        if position as usize >= self.len() {
            return None;
        }

        match self {
            Self::Plain(values) => copied_at(values, position),
            Self::Constant { value, .. } => Some(*value),
            Self::RunLength { values, run_ends } => {
                let run = run_ends.partition_point(|end| *end <= position);

                values.get(run).copied()
            }
            Self::BitPacked { data, width, .. } => {
                unpack_one(data, *width, position).and_then(T::from_bits)
            }
            Self::Delta {
                first,
                deltas,
                checkpoints,
            } => delta_at(*first, deltas, checkpoints, position),
        }
    }

    /// Every value in order.
    ///
    /// Reading sequentially through a delta or run-length column this way is
    /// linear overall, where repeated random access would not be.
    pub fn iter(&self) -> ColumnIter<'_, T> {
        ColumnIter::new(self)
    }

    /// Chooses an encoding for `values` by what they contain.
    ///
    /// The rules, in order: one repeated value becomes constant; few runs
    /// become run-length; small non-negative integers become bit-packed when
    /// that is genuinely narrower; a near-monotonic integer sequence becomes
    /// deltas; anything else stays plain.
    ///
    #[must_use]
    pub fn encode(values: &[T]) -> Self {
        let Some(first) = values.first() else {
            return Self::Plain(Vec::new());
        };

        if values.iter().skip(1).all(|value| value == first)
            && let Ok(len) = u32::try_from(values.len())
        {
            return Self::Constant { value: *first, len };
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

        Self::Plain(copy_values(values))
    }

    /// Stores values literally, whatever they are.
    ///
    /// This is what a hot column asks for: it declines encoding rather than
    /// letting the rules decide.
    #[must_use]
    pub fn plain(values: &[T]) -> Self {
        Self::Plain(copy_values(values))
    }

    fn try_run_length(values: &[T]) -> Option<Self> {
        let runs = run_count(values);

        // Each run costs a value and an end position, so it must replace more
        // than two values' worth of storage to be worth the indirection.
        if runs > values.len() / 4 {
            return None;
        }

        let mut heads = Vec::with_capacity(runs);
        let mut run_ends = Vec::with_capacity(runs);

        for (position, value) in values.iter().enumerate() {
            let end = u32::try_from(position.checked_add(1)?).ok()?;

            if heads.last() == Some(value) {
                if let Some(last_end) = run_ends.last_mut() {
                    *last_end = end;
                }

                continue;
            }

            heads.push(*value);
            run_ends.push(end);
        }

        Some(Self::RunLength {
            values: heads,
            run_ends,
        })
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

        Some(Self::BitPacked {
            data: pack_mapped(values, width, |value| value.to_bits())?,
            width,
            len: u32::try_from(values.len()).ok()?,
        })
    }

    fn try_delta(values: &[T]) -> Option<Self> {
        if !T::IS_INTEGER {
            return None;
        }

        let first = *values.first()?;
        let delta_count = if values.is_empty() {
            0
        } else {
            values.len() - 1
        };
        let mut deltas = Vec::with_capacity(delta_count);

        for pair in values.windows(2) {
            let [earlier, later] = pair else {
                return None;
            };

            let step = later
                .to_bits()
                .cast_signed()
                .checked_sub(earlier.to_bits().cast_signed())?;

            deltas.push(i32::try_from(step).ok()?);
        }

        let mut checkpoints = Vec::with_capacity(deltas.len() / DELTA_CHECKPOINT_STRIDE);
        let mut running = first.to_bits().cast_signed();
        for (index, delta) in deltas.iter().enumerate() {
            running = running.checked_add(i64::from(*delta))?;
            if (index + 1).is_multiple_of(DELTA_CHECKPOINT_STRIDE) {
                checkpoints.push(running);
            }
        }

        Some(Self::Delta {
            first,
            deltas,
            checkpoints,
        })
    }
}

fn run_count<T: PartialEq>(values: &[T]) -> usize {
    if values.is_empty() {
        return 0;
    }

    1 + values
        .windows(2)
        .filter(|pair| pair.first() != pair.get(1))
        .count()
}

fn copy_values<T: Copy>(values: &[T]) -> Vec<T> {
    values.to_vec()
}

fn copied_at<T: Copy>(values: &[T], position: u32) -> Option<T> {
    usize::try_from(position)
        .ok()
        .and_then(|index| values.get(index))
        .copied()
}

/// Positions between stored absolute values in a delta-encoded column.
///
/// Small enough that a read sums at most this many deltas, large enough that
/// the checkpoint table stays a small fraction of the deltas themselves.
const DELTA_CHECKPOINT_STRIDE: usize = 64;

fn delta_at<T: ColumnValue>(
    first: T,
    deltas: &[i32],
    checkpoints: &[i64],
    position: u32,
) -> Option<T> {
    let end = usize::try_from(position).ok()?;
    if end > deltas.len() {
        return None;
    }

    // Resume from the last checkpoint at or before `position` rather than from
    // the head of the column, bounding the walk to one stride.
    let block = end / DELTA_CHECKPOINT_STRIDE;
    let (start, base) = match block.checked_sub(1).and_then(|slot| checkpoints.get(slot)) {
        Some(checkpoint) => (block * DELTA_CHECKPOINT_STRIDE, *checkpoint),
        None => (0, first.to_bits().cast_signed()),
    };

    let value = deltas
        .get(start..end)?
        .iter()
        .fold(base, |value, delta| value + i64::from(*delta));

    T::from_bits(value.cast_unsigned())
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
