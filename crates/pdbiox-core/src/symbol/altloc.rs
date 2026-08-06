//! Alternate-location labels.

use super::interner::SymbolId;
use std::fmt;

/// An alternate-location label.
///
/// The blank label means "this atom belongs to every conformation" — a residue
/// routinely mixes blank-labelled backbone atoms with labelled side-chain ones,
/// which is exactly why picking the highest-occupancy atom one at a time can
/// assemble a conformation that was never modelled.
///
/// A label of any length is representable, because refusing to truncate applies
/// to labels as much as to chain identifiers.
///
/// # Examples
///
/// ```
/// use pdbiox_core::{AltId, Interner};
///
/// let mut interner = Interner::new();
/// let a = AltId::labelled(interner.intern("A")?);
/// assert!(!a.is_blank());
/// assert!(AltId::BLANK.is_blank());
/// # Ok::<(), pdbiox_core::DictionaryFull>(())
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
#[repr(transparent)]
pub struct AltId(u32);

impl AltId {
    /// The label carried by an atom present in every conformation.
    pub const BLANK: Self = Self(0);

    /// Creates a label from an interned string.
    #[must_use]
    pub const fn labelled(symbol: SymbolId) -> Self {
        Self(symbol.get().saturating_add(1))
    }

    /// Returns true for the blank label.
    #[must_use]
    pub const fn is_blank(self) -> bool {
        self.0 == 0
    }

    /// The interned label, or `None` when blank.
    #[must_use]
    pub const fn symbol(self) -> Option<SymbolId> {
        match self.0.checked_sub(1) {
            Some(raw) => Some(SymbolId::from_raw(raw)),
            None => None,
        }
    }

    /// The raw value, for storage in a column.
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0
    }

    /// Rebuilds a label from a stored value.
    #[must_use]
    pub const fn from_raw(raw: u32) -> Self {
        Self(raw)
    }
}

impl fmt::Debug for AltId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.symbol() {
            Some(symbol) => write!(f, "AltId({symbol:?})"),
            None => f.write_str("AltId(blank)"),
        }
    }
}

#[cfg(test)]
#[path = "altloc_tests.rs"]
mod tests;
