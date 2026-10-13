//! How a computation asks for a table: by profile, by name or with the table in hand.
//!
//! Every option struct that needs a reference table holds one [`TableRequest`] per
//! table instead of a table. `Profile` (the default) lets the store's profile file
//! decide, `Named` is a path or a cache entry name, and `Given` carries a table the
//! caller already decoded (a test fixture, bytes from Python). Resolving never
//! picks silently: an unresolvable request is a typed [`TableError`].
//!
//! Cost: resolving a `Given` request is an `Arc` clone; a `Named` or `Profile`
//! request is one profile read at most, one `stat` and one hash probe once the
//! table has been decoded.

use crate::{Provided, TableError, TableKind, TableStore};
use std::sync::Arc;

/// A request for one table of kind `T`.
pub enum TableRequest<T> {
    /// Whatever the store's profile names for this kind.
    Profile,
    /// A file path or a cache entry name.
    Named(String),
    /// A table already decoded.
    Given(Arc<Provided<T>>),
}

impl<T> Default for TableRequest<T> {
    fn default() -> Self {
        Self::profile()
    }
}

impl<T> Clone for TableRequest<T> {
    fn clone(&self) -> Self {
        match self {
            Self::Profile => Self::Profile,
            Self::Named(name) => Self::Named(name.clone()),
            Self::Given(table) => Self::Given(Arc::clone(table)),
        }
    }
}

impl<T> core::fmt::Debug for TableRequest<T> {
    fn fmt(&self, out: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Profile => out.write_str("Profile"),
            Self::Named(name) => out.debug_tuple("Named").field(name).finish(),
            Self::Given(table) => out.debug_tuple("Given").field(&table.identity).finish(),
        }
    }
}

impl<T> TableRequest<T> {
    /// A request for the profile's choice.
    #[must_use]
    pub const fn profile() -> Self {
        Self::Profile
    }

    /// A request for the table in `provided`.
    #[must_use]
    pub fn given(provided: Provided<T>) -> Self {
        Self::Given(Arc::new(provided))
    }

    /// A request for a path or cache entry name.
    #[must_use]
    pub fn named(name: impl Into<String>) -> Self {
        Self::Named(name.into())
    }
}

impl<T: TableKind + Send + Sync + 'static> TableRequest<T> {
    /// The table this request selects.
    ///
    /// # Errors
    ///
    /// [`TableError::Missing`] when nothing selects a table; otherwise the error
    /// of reading or validating the selected file.
    pub fn resolve(
        &self,
        store: &TableStore,
        option: &'static str,
    ) -> Result<Arc<Provided<T>>, TableError> {
        match self {
            Self::Given(table) => Ok(Arc::clone(table)),
            Self::Named(name) => Ok(store.select::<T>(Some(name), option)?.provided),
            Self::Profile => Ok(store.select::<T>(None, option)?.provided),
        }
    }
}

#[cfg(test)]
#[path = "table_request_tests.rs"]
mod tests;
