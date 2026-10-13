//! Name, version, source and licence of a reference table.
//!
//! A number is only as trustworthy as its provenance, so every table that
//! feeds a result, whether it ships with the library or comes from the caller,
//! is described by the same four facts. Results carry the identity of the
//! table they were computed with, which is what lets a reader reproduce them
//! or tell why two runs disagree.

use core::fmt;
use std::borrow::Cow;

/// The four facts that identify a table.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct TableIdentity {
    /// Stable machine-readable name.
    pub name: Cow<'static, str>,
    /// Version or edition label.
    pub version: Cow<'static, str>,
    /// Where the values come from.
    pub source: Cow<'static, str>,
    /// Licence the values are distributed under.
    pub licence: Cow<'static, str>,
}

impl TableIdentity {
    /// An identity made of static text, for tables compiled into the library.
    #[must_use]
    pub const fn builtin(
        name: &'static str,
        version: &'static str,
        source: &'static str,
        licence: &'static str,
    ) -> Self {
        Self {
            name: Cow::Borrowed(name),
            version: Cow::Borrowed(version),
            source: Cow::Borrowed(source),
            licence: Cow::Borrowed(licence),
        }
    }

    /// An identity read from a caller's table.
    #[must_use]
    pub fn supplied(name: &str, version: &str, source: &str, licence: &str) -> Self {
        Self {
            name: Cow::Owned(name.to_owned()),
            version: Cow::Owned(version.to_owned()),
            source: Cow::Owned(source.to_owned()),
            licence: Cow::Owned(licence.to_owned()),
        }
    }
}

impl fmt::Display for TableIdentity {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            out,
            "{} {} (source: {}; licence: {})",
            self.name, self.version, self.source, self.licence
        )
    }
}

#[cfg(test)]
#[path = "table_identity_tests.rs"]
mod tests;
