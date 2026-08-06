//! Why a value is not there.
//!
//! `Option` collapses three different statements into one. A file can say a
//! field does not apply here, that its value exists but was not recorded, or
//! simply give a value — and those mean different things to an analysis. A
//! B-factor that does not apply is not a B-factor that nobody measured.
//!
//! A fourth case is deliberately *not* represented here: an atom that was never
//! modelled. That atom has no row at all. Inventing a row with an absent
//! position would fabricate an atom the experiment never saw, and the gap
//! between what a component should contain and what the model contains is
//! reported as coverage instead.

use super::bits::BitVec;

/// What a column says about one position.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default)]
pub enum Presence {
    /// A value is recorded.
    #[default]
    Present,
    /// The value exists but was not recorded.
    Unknown,
    /// The field does not apply here.
    Inapplicable,
}

impl Presence {
    /// Returns true when a value is recorded.
    #[must_use]
    pub const fn is_present(self) -> bool {
        matches!(self, Self::Present)
    }
}

/// Which positions of a column carry a value, and why the others do not.
///
/// The overwhelmingly common case is that every position carries one, and that
/// case costs a length rather than two bit sets.
///
/// # Examples
///
/// ```
/// use pdbiox_core::column::{Presence, ValidityMask};
///
/// let mut mask = ValidityMask::all_present(3);
/// assert_eq!(mask.get(0), Presence::Present);
///
/// mask.set(1, Presence::Unknown);
/// assert_eq!(mask.get(1), Presence::Unknown);
/// assert_eq!(mask.present_count(), 2);
/// ```
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum ValidityMask {
    /// Every position carries a value.
    AllPresent(u32),
    /// Positions differ, so both distinctions are tracked.
    Mixed {
        /// Set where a value is recorded.
        present: BitVec,
        /// Set where the value exists but was not recorded. Meaningful only
        /// where `present` is clear; clear there means the field does not apply.
        unknown: BitVec,
    },
}

impl ValidityMask {
    /// Creates a mask over `len` positions that all carry a value.
    #[must_use]
    pub const fn all_present(len: u32) -> Self {
        Self::AllPresent(len)
    }

    /// The number of positions covered.
    #[must_use]
    pub fn len(&self) -> u32 {
        match self {
            Self::AllPresent(len) => *len,
            Self::Mixed { present, .. } => present.len(),
        }
    }

    /// Returns true when no position is covered.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Returns true when every covered position carries a value.
    #[must_use]
    pub fn is_all_present(&self) -> bool {
        match self {
            Self::AllPresent(_) => true,
            Self::Mixed { present, .. } => present.all(),
        }
    }

    /// What the column says about `position`.
    ///
    /// A position past the end reads as inapplicable, because a column that does
    /// not reach that far makes no statement about it.
    #[must_use]
    pub fn get(&self, position: u32) -> Presence {
        match self {
            Self::AllPresent(len) => {
                if position < *len {
                    Presence::Present
                } else {
                    Presence::Inapplicable
                }
            }
            Self::Mixed { present, unknown } => {
                if present.test(position) {
                    Presence::Present
                } else if unknown.test(position) {
                    Presence::Unknown
                } else {
                    Presence::Inapplicable
                }
            }
        }
    }

    /// Records what the column says about `position`.
    ///
    /// Setting anything other than present on an all-present mask materialises
    /// the bit sets, which is the only point at which they are allocated.
    pub fn set(&mut self, position: u32, presence: Presence) {
        if let Self::AllPresent(len) = *self {
            if presence.is_present() || position >= len {
                return;
            }
            *self = Self::Mixed {
                present: BitVec::repeat(true, len),
                unknown: BitVec::repeat(false, len),
            };
        }
        if let Self::Mixed { present, unknown } = self {
            present.set(position, presence == Presence::Present);
            unknown.set(position, presence == Presence::Unknown);
        }
    }

    /// The number of positions carrying a value.
    #[must_use]
    pub fn present_count(&self) -> u32 {
        match self {
            Self::AllPresent(len) => *len,
            Self::Mixed { present, .. } => present.count_ones(),
        }
    }

    /// Collapses back to the compact form when every position carries a value.
    ///
    /// Worth calling once after building a column, so a file that happened to
    /// record everything does not carry two bit sets for the rest of its life.
    pub fn compact(&mut self) {
        if let Self::Mixed { present, .. } = self
            && present.all()
        {
            *self = Self::AllPresent(present.len());
        }
    }
}

impl FromIterator<Presence> for ValidityMask {
    fn from_iter<T: IntoIterator<Item = Presence>>(iter: T) -> Self {
        let mut present = BitVec::new();
        let mut unknown = BitVec::new();
        for presence in iter {
            present.push(presence == Presence::Present);
            unknown.push(presence == Presence::Unknown);
        }
        let mut mask = Self::Mixed { present, unknown };
        mask.compact();
        mask
    }
}

#[cfg(test)]
#[path = "validity_tests.rs"]
mod tests;
