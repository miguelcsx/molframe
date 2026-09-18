//! Typed failures at compact-storage representation boundaries.

use std::fmt;

/// A collection cannot represent another row without losing information.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CapacityError {
    resource: &'static str,
}

impl CapacityError {
    /// Creates a capacity error for a named collection.
    #[must_use]
    pub const fn new(resource: &'static str) -> Self {
        Self { resource }
    }

    /// The affected collection.
    #[must_use]
    pub const fn resource(self) -> &'static str {
        self.resource
    }
}

impl fmt::Display for CapacityError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{} exceeds its supported capacity",
            self.resource
        )
    }
}

impl std::error::Error for CapacityError {}

/// A hierarchy table rejected an append operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TableError {
    /// The compact table representation is full.
    Capacity(CapacityError),
    /// The supplied child range ends before it starts.
    ReversedRange {
        /// The affected hierarchy table.
        table: &'static str,
    },
    /// The requested table position does not exist.
    MissingIndex {
        /// The affected hierarchy table.
        table: &'static str,
    },
}

impl fmt::Display for TableError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Capacity(error) => error.fmt(formatter),
            Self::ReversedRange { table } => write!(formatter, "{table} child range is reversed"),
            Self::MissingIndex { table } => write!(formatter, "{table} index does not exist"),
        }
    }
}

impl std::error::Error for TableError {}
