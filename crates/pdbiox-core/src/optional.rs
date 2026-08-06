//! Absent-or-present values that cost no more than a present one.
//!
//! A great many identifier columns are genuinely optional: a file may or may not
//! carry the depositor's own chain label, and a component that is not part of a
//! polymer has no sequence position. `Option` doubles the width of those columns
//! because a 32-bit identifier has no spare bit pattern for the tag.
//!
//! These types reserve one value instead. At a residue table of a hundred
//! thousand rows and three such columns, that is the difference between one and
//! two megabytes — memory that geometric kernels would otherwise have as cache.
//!
//! The reserved values are chosen where no real datum can land: an identifier
//! is never `u32::MAX` because the dictionary refuses to grow that far, and a
//! sequence position is never `i32::MIN`.

use crate::symbol::SymbolId;
use std::fmt;

/// An interned identifier that may be absent.
///
/// # Examples
///
/// ```
/// use pdbiox_core::{Interner, OptionalSymbol};
///
/// let mut interner = Interner::new();
/// let present = OptionalSymbol::some(interner.intern("A")?);
///
/// assert!(present.get().is_some());
/// assert_eq!(OptionalSymbol::NONE.get(), None);
/// assert_eq!(size_of::<OptionalSymbol>(), 4);
/// # Ok::<(), pdbiox_core::DictionaryFull>(())
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct OptionalSymbol(u32);

impl OptionalSymbol {
    /// The value is absent.
    pub const NONE: Self = Self(u32::MAX);

    /// Wraps a present identifier.
    #[must_use]
    pub const fn some(symbol: SymbolId) -> Self {
        Self(symbol.get())
    }

    /// The identifier, if present.
    #[must_use]
    pub const fn get(self) -> Option<SymbolId> {
        if self.0 == u32::MAX {
            None
        } else {
            Some(SymbolId::from_raw(self.0))
        }
    }

    /// Returns true when no identifier is held.
    #[must_use]
    pub const fn is_none(self) -> bool {
        self.0 == u32::MAX
    }

    /// Returns true when an identifier is held.
    #[must_use]
    pub const fn is_some(self) -> bool {
        self.0 != u32::MAX
    }
}

impl Default for OptionalSymbol {
    fn default() -> Self {
        Self::NONE
    }
}

impl From<Option<SymbolId>> for OptionalSymbol {
    fn from(symbol: Option<SymbolId>) -> Self {
        match symbol {
            Some(symbol) => Self::some(symbol),
            None => Self::NONE,
        }
    }
}

impl From<OptionalSymbol> for Option<SymbolId> {
    fn from(optional: OptionalSymbol) -> Self {
        optional.get()
    }
}

impl fmt::Debug for OptionalSymbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.get() {
            Some(symbol) => write!(f, "Some({symbol:?})"),
            None => f.write_str("None"),
        }
    }
}

/// A signed number that may be absent.
///
/// Sequence positions and the depositor's own residue numbers use this: both
/// are genuinely optional, and both may legitimately be negative, so the
/// reserved value is the one no numbering scheme reaches.
///
/// # Examples
///
/// ```
/// use pdbiox_core::OptionalI32;
///
/// assert_eq!(OptionalI32::some(-5).get(), Some(-5));
/// assert_eq!(OptionalI32::NONE.get(), None);
/// assert_eq!(size_of::<OptionalI32>(), 4);
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct OptionalI32(i32);

impl OptionalI32 {
    /// The value is absent.
    pub const NONE: Self = Self(i32::MIN);

    /// Wraps a present value.
    ///
    /// The reserved value maps to absent, which is the only lossy case and is
    /// unreachable for any real numbering.
    #[must_use]
    pub const fn some(value: i32) -> Self {
        Self(value)
    }

    /// The value, if present.
    #[must_use]
    pub const fn get(self) -> Option<i32> {
        if self.0 == i32::MIN {
            None
        } else {
            Some(self.0)
        }
    }

    /// Returns true when no value is held.
    #[must_use]
    pub const fn is_none(self) -> bool {
        self.0 == i32::MIN
    }

    /// Returns true when a value is held.
    #[must_use]
    pub const fn is_some(self) -> bool {
        self.0 != i32::MIN
    }
}

impl Default for OptionalI32 {
    fn default() -> Self {
        Self::NONE
    }
}

impl From<Option<i32>> for OptionalI32 {
    fn from(value: Option<i32>) -> Self {
        match value {
            Some(value) => Self::some(value),
            None => Self::NONE,
        }
    }
}

impl From<OptionalI32> for Option<i32> {
    fn from(optional: OptionalI32) -> Self {
        optional.get()
    }
}

impl fmt::Debug for OptionalI32 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.get() {
            Some(value) => write!(f, "Some({value})"),
            None => f.write_str("None"),
        }
    }
}

#[cfg(test)]
#[path = "optional_tests.rs"]
mod tests;
