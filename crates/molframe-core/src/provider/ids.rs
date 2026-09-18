//! Stable global identities and compact chunk-local rows.

use std::fmt;

macro_rules! global_id {
    ($name:ident, $label:literal) => {
        #[doc = concat!("Stable 64-bit ", $label, ".")]
        #[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
        #[repr(transparent)]
        pub struct $name(u64);

        impl $name {
            #[doc = concat!("Creates a ", $label, " from its stable value.")]
            #[must_use]
            pub const fn new(value: u64) -> Self {
                Self(value)
            }

            #[doc = concat!("Returns the raw ", $label, " value.")]
            #[must_use]
            pub const fn get(self) -> u64 {
                self.0
            }

            #[doc = concat!("Returns the next ", $label, ", if representable.")]
            #[must_use]
            pub const fn next(self) -> Option<Self> {
                match self.0.checked_add(1) {
                    Some(value) => Some(Self(value)),
                    None => None,
                }
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(formatter, "{}({})", $label, self.0)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(formatter, "{}", self.0)
            }
        }
    };
}

global_id!(DatasetId, "dataset");
global_id!(ChunkId, "chunk");
global_id!(LogicalRow, "logical row");

/// Compact row index used only inside one resident chunk.
#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct LocalRow(u32);

impl LocalRow {
    /// Creates a local row from its chunk-relative ordinal.
    #[must_use]
    pub const fn new(value: u32) -> Self {
        Self(value)
    }

    /// Returns the chunk-relative ordinal.
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0
    }

    /// Returns the next local row, if representable.
    #[must_use]
    pub const fn next(self) -> Option<Self> {
        match self.0.checked_add(1) {
            Some(value) => Some(Self(value)),
            None => None,
        }
    }
}

impl fmt::Debug for LocalRow {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "local row({})", self.0)
    }
}

impl fmt::Display for LocalRow {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}", self.0)
    }
}
