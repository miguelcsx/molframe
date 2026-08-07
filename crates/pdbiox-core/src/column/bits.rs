//! Bit-level primitives: a compact bit set, and packing narrow integers.
//!
//! Both exist for the same reason. A per-atom flag stored as a `bool` costs a
//! byte, eight times what it needs; a per-atom value with a range of a hundred
//! costs four bytes as a `u32`, four times what it needs. At a million atoms
//! those factors are the difference between a column that stays in cache during
//! a scan and one that does not.
//!
//! Neither is used for positions. Unpacking costs a shift and a mask per value,
//! which is worth paying for a column read once per atom and not worth paying
//! for one a kernel reads on every iteration.

use std::sync::Arc;

/// A compact set of bits.
///
/// Iteration yields set positions rather than every position, so walking a
/// sparse mask costs one step per member instead of one per candidate.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct BitVec {
    words: Arc<Vec<u64>>,
    len: u32,
}

const BITS: u32 = u64::BITS;

impl BitVec {
    /// Creates an empty set.
    #[must_use]
    pub fn new() -> Self {
        Self {
            words: Arc::new(Vec::new()),
            len: 0,
        }
    }

    /// Creates a set of `len` bits, all equal to `value`.
    #[must_use]
    pub fn repeat(value: bool, len: u32) -> Self {
        let word = if value { u64::MAX } else { 0 };
        let words = Arc::new(vec![word; word_count(len)]);
        let mut bits = Self { words, len };

        bits.clear_tail();
        bits
    }

    /// Creates a set with room for `len` bits without reallocating.
    #[must_use]
    pub fn with_capacity(len: u32) -> Self {
        Self {
            words: Arc::new(Vec::with_capacity(word_count(len))),
            len: 0,
        }
    }

    /// The number of bits held.
    #[must_use]
    pub const fn len(&self) -> u32 {
        self.len
    }

    /// Returns true when no bit is held.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Appends a bit.
    pub fn push(&mut self, value: bool) {
        let position = self.len;

        if position.is_multiple_of(BITS) {
            Arc::make_mut(&mut self.words).push(0);
        }

        self.len += 1;

        if !value {
            return;
        }

        if let Some(word) = Arc::make_mut(&mut self.words).last_mut() {
            *word |= bit_mask(position);
        }
    }

    /// Returns the bit at `position`, or `None` if it is past the end.
    #[must_use]
    pub fn get(&self, position: u32) -> Option<bool> {
        if position >= self.len {
            return None;
        }

        self.words
            .get(word_index(position))
            .map(|word| *word & bit_mask(position) != 0)
    }

    /// Returns the bit at `position`, treating a position past the end as unset.
    #[must_use]
    pub fn test(&self, position: u32) -> bool {
        self.get(position) == Some(true)
    }

    /// Sets the bit at `position`, ignoring a position past the end.
    pub fn set(&mut self, position: u32, value: bool) {
        if position >= self.len {
            return;
        }
        let Some(word) = Arc::make_mut(&mut self.words).get_mut(word_index(position)) else {
            return;
        };

        let mask = bit_mask(position);

        if value {
            *word |= mask;
        } else {
            *word &= !mask;
        }
    }

    /// The number of set bits.
    #[must_use]
    pub fn count_ones(&self) -> u32 {
        self.words.iter().map(|word| word.count_ones()).sum()
    }

    /// Returns true when every bit is set.
    #[must_use]
    pub fn all(&self) -> bool {
        let complete_words = (self.len / BITS) as usize;

        if !self
            .words
            .iter()
            .take(complete_words)
            .all(|word| *word == u64::MAX)
        {
            return false;
        }

        let tail_bits = self.len % BITS;

        tail_bits == 0 || self.words.get(complete_words).copied() == Some(low_mask(tail_bits))
    }

    /// Returns true when no bit is set.
    #[must_use]
    pub fn none(&self) -> bool {
        self.words.iter().all(|word| *word == 0)
    }

    /// The set positions, ascending.
    ///
    /// Skips whole words that hold nothing, so a mask with one member in a
    /// million bits costs a scan of sixteen thousand words rather than a
    /// million tests.
    pub fn ones(&self) -> impl Iterator<Item = u32> + '_ {
        self.words.iter().enumerate().flat_map(|(index, word)| {
            let base = index as u32 * BITS;
            let mut remaining = *word;

            std::iter::from_fn(move || {
                if remaining == 0 {
                    return None;
                }

                let bit = remaining.trailing_zeros();
                remaining &= remaining - 1;

                Some(base + bit)
            })
        })
    }

    /// Keeps only the bits also set in `other`.
    pub fn intersect_with(&mut self, other: &Self) {
        let shared = self.words.len().min(other.words.len());
        let words = Arc::make_mut(&mut self.words);

        for (word, other_word) in words[..shared].iter_mut().zip(&other.words[..shared]) {
            *word &= *other_word;
        }

        for word in &mut words[shared..] {
            *word = 0;
        }
    }

    /// Adds every bit set in `other`.
    ///
    /// Bits of `other` beyond this set's length are ignored; the length is a
    /// property of the structure, not of the operation.
    pub fn union_with(&mut self, other: &Self) {
        for (word, other_word) in Arc::make_mut(&mut self.words)
            .iter_mut()
            .zip(other.words.iter())
        {
            *word |= *other_word;
        }

        self.clear_tail();
    }

    /// Inverts every bit.
    pub fn invert(&mut self) {
        for word in Arc::make_mut(&mut self.words) {
            *word = !*word;
        }

        self.clear_tail();
    }

    /// Zeroes the bits past the end that live in the final word, so that
    /// counting and comparison never see them.
    fn clear_tail(&mut self) {
        let used = self.len % BITS;

        if used == 0 {
            return;
        }
        if let Some(word) = Arc::make_mut(&mut self.words).last_mut() {
            *word &= low_mask(used);
        }
    }
}

impl FromIterator<bool> for BitVec {
    fn from_iter<T: IntoIterator<Item = bool>>(iter: T) -> Self {
        let iter = iter.into_iter();
        let capacity = match u32::try_from(iter.size_hint().0) {
            Ok(capacity) => capacity,
            Err(_) => u32::MAX,
        };
        let mut bits = Self::with_capacity(capacity);

        for value in iter {
            bits.push(value);
        }

        bits
    }
}

/// The number of bits needed to represent every value up to `max`.
#[must_use]
pub fn bit_width(max: u64) -> u8 {
    match max {
        0 => 1,
        _ => (u64::BITS - max.leading_zeros()) as u8,
    }
}

/// Mask covering the low `width` bits.
const fn low_mask(width: u32) -> u64 {
    if width >= u64::BITS {
        u64::MAX
    } else {
        (1u64 << width) - 1
    }
}

fn word_count(len: u32) -> usize {
    len.div_ceil(BITS) as usize
}

const fn word_index(position: u32) -> usize {
    (position / BITS) as usize
}

const fn bit_offset(position: u32) -> u32 {
    position % BITS
}

const fn bit_mask(position: u32) -> u64 {
    1u64 << bit_offset(position)
}

/// Where a packed value sits: which byte it starts in, how far into that byte,
/// and how many bytes it touches.
const fn placement(width: u32, index: u32) -> (usize, u32, usize) {
    let start_bit = index as u64 * width as u64;
    let byte = (start_bit / 8) as usize;
    let offset = (start_bit % 8) as u32;
    let bytes = (offset as usize + width as usize).div_ceil(8);

    (byte, offset, bytes)
}

/// Packs `values` into `width` bits each.
///
/// Values wider than `width` are truncated, which is why the caller derives the
/// width from the data rather than choosing one. Each value is written with a
/// single shifted accumulation rather than a loop over its bits.
#[must_use]
pub fn pack(values: &[u64], width: u8) -> Vec<u8> {
    pack_mapped(values, width, |value| *value)
}

pub(super) fn pack_mapped<T, F>(values: &[T], width: u8, to_bits: F) -> Vec<u8>
where
    F: Fn(&T) -> u64,
{
    let width = u32::from(width.clamp(1, 64));
    let total_bits = values.len() as u64 * u64::from(width);
    let mut packed = vec![0u8; total_bits.div_ceil(8) as usize];

    for (index, value) in values.iter().enumerate() {
        write_packed_value(&mut packed, width, index as u32, to_bits(value));
    }

    packed
}

fn write_packed_value(packed: &mut [u8], width: u32, index: u32, value: u64) {
    let (byte, offset, bytes) = placement(width, index);
    let Some(end) = byte.checked_add(bytes) else {
        return;
    };
    let Some(targets) = packed.get_mut(byte..end) else {
        return;
    };

    let shifted = u128::from(value & low_mask(width)) << offset;

    for (step, target) in targets.iter_mut().enumerate() {
        *target |= (shifted >> (8 * step)) as u8;
    }
}

/// Reads the value at `index` from a packed buffer.
///
/// Returns `None` when the value would run past the end of the buffer, which is
/// the only way a malformed input can reach this code.
#[must_use]
pub fn unpack_one(packed: &[u8], width: u8, index: u32) -> Option<u64> {
    let width = u32::from(width.clamp(1, 64));
    let (byte, offset, bytes) = placement(width, index);
    let end = byte.checked_add(bytes)?;
    let source = packed.get(byte..end)?;

    let accumulator = source
        .iter()
        .enumerate()
        .fold(0u128, |accumulator, (step, value)| {
            accumulator | (u128::from(*value) << (8 * step))
        });

    Some(((accumulator >> offset) as u64) & low_mask(width))
}

#[cfg(test)]
#[path = "bits_tests.rs"]
mod tests;
