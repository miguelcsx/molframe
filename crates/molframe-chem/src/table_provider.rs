//! The one trait through which every reference table reaches a computation.
//!
//! A [`TableProvider`] hands out a table together with the [`TableIdentity`] that
//! describes it, by reference: the table is decoded at most once per provider and
//! every later call borrows it (or shares it as an `Arc`), so a lookup never
//! copies a table. The library compiles in no table. The providers differ only in
//! where the bytes come from: [`UserFile`] (a path, read and validated once),
//! [`Cache`] (a directory file stored by an explicit fetch command) and
//! [`InMemory`] (built by the caller). Choosing among them is the job of the
//! [`TableStore`](crate::TableStore), which never substitutes one source for
//! another silently and never touches the network.
//!
//! Cost: decoding is linear in the table size and happens once per provider;
//! `provide` afterwards is a pointer read.

use crate::{TableIdentity, TableTextError};
use core::fmt;
use std::path::PathBuf;
use std::sync::{Arc, OnceLock};

/// A decoded table and the identity under which it was supplied.
#[derive(Clone, Debug, PartialEq)]
pub struct Provided<T> {
    /// Name, version, source and licence of the table.
    pub identity: TableIdentity,
    /// The table itself.
    pub table: T,
}

/// A kind of table: how its bytes decode and what it is called in messages.
pub trait TableKind: Sized {
    /// Short stable name of the kind, for example `geometry`.
    const KIND: &'static str;

    /// File extension of a cache entry of this kind, without the dot.
    const EXTENSION: &'static str;

    /// Decodes and validates `bytes`. A text kind reads its identity from the
    /// header; a binary kind has none and takes `fallback`.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when the bytes are not a well-formed table of this kind.
    fn decode(bytes: &[u8], fallback: TableIdentity) -> Result<Provided<Self>, TableTextError>;
}

/// Where a provided table came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// A file the caller named by path.
    User,
    /// A cache entry, stored there by an explicit fetch.
    Cache,
    /// A table the caller built in memory.
    Memory,
}

impl Source {
    /// The word used in reports.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::User => "user",
            Self::Cache => "cache",
            Self::Memory => "memory",
        }
    }
}

/// Why a table could not be provided.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TableError {
    /// The provider has no table: a cache file nobody fetched.
    Absent {
        /// The kind of table.
        kind: &'static str,
        /// The path that was looked at.
        path: String,
    },
    /// Nothing selects a table: no request, no profile entry, or an entry the
    /// cache does not hold.
    Missing {
        /// The kind of table.
        kind: &'static str,
        /// What was looked for and where.
        wanted: String,
        /// The option that names a file instead.
        option: &'static str,
    },
    /// A file could not be read.
    Unreadable {
        /// The path.
        path: String,
        /// The operating system's reason.
        reason: String,
    },
    /// The bytes are not a well-formed table.
    Malformed {
        /// The path, or the label of the in-memory or baked data.
        path: String,
        /// What is wrong, with its line.
        error: TableTextError,
    },
}

impl fmt::Display for TableError {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Absent { kind, path } => write!(out, "no {kind} table at {path}"),
            Self::Missing {
                kind,
                wanted,
                option,
            } => write!(
                out,
                "no {kind} table: {wanted}; pass a file with `{option}`, select a cache entry \
                 by name, or run `molframe data fetch <entry> --manifest <file>` to store one \
                 in the cache"
            ),
            Self::Unreadable { path, reason } => write!(out, "cannot read {path}: {reason}"),
            Self::Malformed { path, error } => write!(out, "{path}: {error}"),
        }
    }
}

impl core::error::Error for TableError {}

molframe_core::diagnostic_from!(TableError, |error| match error {
    TableError::Absent { .. } | TableError::Missing { .. } => molframe_core::Code::E7102,
    TableError::Unreadable { .. } => molframe_core::Code::E7101,
    TableError::Malformed { .. } => molframe_core::Code::E2201,
});

/// Something that can hand out a table of kind `T` with its identity.
pub trait TableProvider<T> {
    /// The table and its identity, decoded on first use and borrowed after.
    ///
    /// # Errors
    ///
    /// [`TableError`] when the source is absent, unreadable or malformed.
    fn provide(&self) -> Result<&Provided<T>, TableError>;

    /// The same table as a shared handle, for a caller that outlives the provider.
    ///
    /// # Errors
    ///
    /// As [`TableProvider::provide`].
    fn shared(&self) -> Result<Arc<Provided<T>>, TableError>;

    /// Which source this provider is.
    fn source(&self) -> Source;
}

type Cell<T> = OnceLock<Result<Arc<Provided<T>>, TableError>>;

fn load<T>(
    cell: &Cell<T>,
    init: impl FnOnce() -> Result<Provided<T>, TableError>,
) -> Result<&Arc<Provided<T>>, TableError> {
    match cell.get_or_init(|| init().map(Arc::new)) {
        Ok(provided) => Ok(provided),
        Err(error) => Err(error.clone()),
    }
}

/// A table in a file the caller named, read and validated on first use.
pub struct UserFile<T> {
    path: PathBuf,
    cell: Cell<T>,
}

impl<T> UserFile<T> {
    /// Names the file; nothing is read until the table is first asked for.
    #[must_use]
    pub const fn new(path: PathBuf) -> Self {
        Self {
            path,
            cell: OnceLock::new(),
        }
    }
}

/// Reads `path` whole and decodes it; a missing file is `Absent` when `absent_ok`.
fn decode_file<T: TableKind>(
    path: &std::path::Path,
    label: &str,
    absent_ok: bool,
) -> Result<Provided<T>, TableError> {
    let shown = path.display().to_string();
    let bytes = match std::fs::read(path) {
        Ok(bytes) => bytes,
        Err(error) if absent_ok && error.kind() == std::io::ErrorKind::NotFound => {
            return Err(TableError::Absent {
                kind: T::KIND,
                path: shown,
            });
        }
        Err(error) => {
            return Err(TableError::Unreadable {
                path: shown,
                reason: error.to_string(),
            });
        }
    };
    let fallback = TableIdentity::supplied(
        &format!("{label}-{}", T::KIND),
        "supplied",
        &shown,
        "unspecified",
    );
    T::decode(&bytes, fallback).map_err(|error| TableError::Malformed { path: shown, error })
}

impl<T: TableKind> TableProvider<T> for UserFile<T> {
    fn provide(&self) -> Result<&Provided<T>, TableError> {
        Ok(load(&self.cell, || decode_file(&self.path, "user", false))?)
    }

    fn shared(&self) -> Result<Arc<Provided<T>>, TableError> {
        Ok(Arc::clone(load(&self.cell, || {
            decode_file(&self.path, "user", false)
        })?))
    }

    fn source(&self) -> Source {
        Source::User
    }
}

/// A table in the cache directory, stored there by an explicit fetch.
pub struct Cache<T> {
    path: PathBuf,
    cell: Cell<T>,
}

impl<T> Cache<T> {
    /// The cache file `file` under `directory`; nothing is read yet.
    #[must_use]
    pub fn new(directory: &std::path::Path, file: &str) -> Self {
        Self {
            path: directory.join(file),
            cell: OnceLock::new(),
        }
    }

    /// Where the table is expected.
    #[must_use]
    pub fn path(&self) -> &std::path::Path {
        &self.path
    }
}

impl<T: TableKind> TableProvider<T> for Cache<T> {
    fn provide(&self) -> Result<&Provided<T>, TableError> {
        Ok(load(&self.cell, || decode_file(&self.path, "cache", true))?)
    }

    fn shared(&self) -> Result<Arc<Provided<T>>, TableError> {
        Ok(Arc::clone(load(&self.cell, || {
            decode_file(&self.path, "cache", true)
        })?))
    }

    fn source(&self) -> Source {
        Source::Cache
    }
}

/// A table the caller built.
pub struct InMemory<T> {
    provided: Arc<Provided<T>>,
}

impl<T> InMemory<T> {
    /// Wraps a table with its identity.
    #[must_use]
    pub fn new(provided: Provided<T>) -> Self {
        Self {
            provided: Arc::new(provided),
        }
    }
}

impl<T> TableProvider<T> for InMemory<T> {
    fn provide(&self) -> Result<&Provided<T>, TableError> {
        Ok(&self.provided)
    }

    fn shared(&self) -> Result<Arc<Provided<T>>, TableError> {
        Ok(Arc::clone(&self.provided))
    }

    fn source(&self) -> Source {
        Source::Memory
    }
}

impl<T> fmt::Debug for UserFile<T> {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        out.debug_struct("UserFile")
            .field("path", &self.path)
            .finish_non_exhaustive()
    }
}

impl<T> fmt::Debug for Cache<T> {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        out.debug_struct("Cache")
            .field("path", &self.path)
            .finish_non_exhaustive()
    }
}

impl<T> fmt::Debug for InMemory<T> {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        out.debug_struct("InMemory")
            .field("identity", &self.provided.identity)
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
#[path = "table_provider_tests.rs"]
mod tests;
