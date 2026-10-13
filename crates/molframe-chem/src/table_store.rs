//! Where tables are looked for, and which one a computation gets.
//!
//! A [`TableStore`] owns the cache directory and the profile, and answers one
//! question: given a request for a table of some kind, which file is it? The
//! request is what the caller passed (`--geometry-table`, a Python argument); the
//! profile is a plain text file of `kind: name` lines that says which cache entry a
//! kind uses when the caller passed nothing. Nothing is built in. A request that
//! names a path is read from that path; a request that names an entry is read
//! from `<directory>/<name>.<extension>`; no request and no profile line is a typed
//! [`TableError::Missing`], never a guess.
//!
//! Decoded tables are memoised per (kind, path, size, modification time), so a
//! table is parsed once however many computations ask for it and an edited file
//! is read again. Lookups after the first are one hash probe and an `Arc` clone.
//!
//! Directory resolution order: `MOLFRAME_DATA_DIR`, then `XDG_CACHE_HOME/molframe`,
//! then `~/.cache/molframe` (`%LOCALAPPDATA%\molframe` on Windows). The profile is
//! `MOLFRAME_PROFILE` or `<directory>/profile.txt`.

use crate::{Cache, UserFile};
use crate::{Provided, Source, TableError, TableKind, TableProvider};
use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::{Arc, Mutex};
use std::time::SystemTime;

/// A table chosen for a computation, with where it came from.
pub struct Selected<T> {
    /// The table and its identity, shared with the store's memo.
    pub provided: Arc<Provided<T>>,
    /// Which source supplied it.
    pub source: Source,
}

impl<T> Clone for Selected<T> {
    fn clone(&self) -> Self {
        Self {
            provided: Arc::clone(&self.provided),
            source: self.source,
        }
    }
}

impl<T> core::fmt::Debug for Selected<T> {
    fn fmt(&self, out: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        out.debug_struct("Selected")
            .field("identity", &self.provided.identity)
            .field("source", &self.source)
            .finish()
    }
}

type Stamp = (u64, Option<SystemTime>);
type Memo = HashMap<(TypeId, PathBuf), (Stamp, Arc<dyn Any + Send + Sync>)>;

/// The cache directory, the profile and the memo of decoded tables.
pub struct TableStore {
    directory: PathBuf,
    profile: PathBuf,
    memo: Mutex<Memo>,
}

impl core::fmt::Debug for TableStore {
    fn fmt(&self, out: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        out.debug_struct("TableStore")
            .field("directory", &self.directory)
            .field("profile", &self.profile)
            .finish_non_exhaustive()
    }
}

impl TableStore {
    /// A store over `directory`, with the profile at `directory/profile.txt`.
    #[must_use]
    pub fn new(directory: PathBuf) -> Self {
        let profile = directory.join("profile.txt");
        Self {
            directory,
            profile,
            memo: Mutex::new(HashMap::new()),
        }
    }

    /// The process-wide store over the environment's directory and profile.
    #[must_use]
    pub fn global() -> Arc<Self> {
        static STORE: OnceLock<Arc<TableStore>> = OnceLock::new();
        Arc::clone(STORE.get_or_init(|| Arc::new(Self::from_environment())))
    }

    /// A store over the directory the environment names, with its profile.
    #[must_use]
    pub fn from_environment() -> Self {
        let mut store = Self::new(environment_directory());
        if let Some(path) = std::env::var_os("MOLFRAME_PROFILE") {
            store.profile = PathBuf::from(path);
        }
        store
    }

    /// The cache directory.
    #[must_use]
    pub fn directory(&self) -> &Path {
        &self.directory
    }

    /// Replaces the profile file.
    #[must_use]
    pub fn with_profile(mut self, profile: PathBuf) -> Self {
        self.profile = profile;
        self
    }

    /// The table `request` selects, or the profile's choice for `T` when it is `None`.
    ///
    /// A request containing a path separator, or naming an existing file, is a
    /// path; anything else is a cache entry name.
    ///
    /// # Errors
    ///
    /// [`TableError::Missing`] when nothing selects a table or the cache lacks the
    /// entry; [`TableError::Unreadable`] or [`TableError::Malformed`] for a file
    /// that cannot be read or fails validation.
    pub fn select<T>(
        &self,
        request: Option<&str>,
        option: &'static str,
    ) -> Result<Selected<T>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        let (name, source) = match request {
            Some(text) if is_path(text) => {
                return self.load::<T>(PathBuf::from(text), Source::User, option);
            }
            Some(text) => (text.to_owned(), Source::Cache),
            None => (self.profile_entry::<T>(option)?, Source::Cache),
        };
        let path = self.directory.join(format!("{name}.{}", T::EXTENSION));
        self.load::<T>(path, source, option)
    }

    fn profile_entry<T: TableKind>(&self, option: &'static str) -> Result<String, TableError> {
        let missing = |wanted: String| TableError::Missing {
            kind: T::KIND,
            wanted,
            option,
        };
        let text = match std::fs::read_to_string(&self.profile) {
            Ok(text) => text,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                return Err(missing(format!(
                    "no request was made and there is no profile at {}",
                    self.profile.display()
                )));
            }
            Err(error) => {
                return Err(TableError::Unreadable {
                    path: self.profile.display().to_string(),
                    reason: error.to_string(),
                });
            }
        };
        for line in text.lines() {
            let line = line.trim();
            if let Some((key, value)) = line.split_once(':')
                && !line.starts_with('#')
                && key.trim() == T::KIND
                && !value.trim().is_empty()
            {
                return Ok(value.trim().to_owned());
            }
        }
        Err(missing(format!(
            "no request was made and {} has no `{}:` line",
            self.profile.display(),
            T::KIND
        )))
    }

    fn load<T>(
        &self,
        path: PathBuf,
        source: Source,
        option: &'static str,
    ) -> Result<Selected<T>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        let stamp = match std::fs::metadata(&path) {
            Ok(meta) => (meta.len(), meta.modified().ok()),
            Err(error)
                if error.kind() == std::io::ErrorKind::NotFound && source == Source::Cache =>
            {
                return Err(TableError::Missing {
                    kind: T::KIND,
                    wanted: format!("the cache holds no entry at {}", path.display()),
                    option,
                });
            }
            Err(error) => {
                return Err(TableError::Unreadable {
                    path: path.display().to_string(),
                    reason: error.to_string(),
                });
            }
        };
        let key = (TypeId::of::<T>(), path.clone());
        if let Ok(memo) = self.memo.lock()
            && let Some((seen, any)) = memo.get(&key)
            && *seen == stamp
            && let Ok(provided) = Arc::clone(any).downcast::<Provided<T>>()
        {
            return Ok(Selected { provided, source });
        }
        let provided = match source {
            Source::Cache => Cache::<T>::new(&self.directory, &file_name(&path)).shared()?,
            _ => UserFile::<T>::new(path).shared()?,
        };
        if let Ok(mut memo) = self.memo.lock() {
            memo.insert(
                key,
                (stamp, Arc::clone(&provided) as Arc<dyn Any + Send + Sync>),
            );
        }
        Ok(Selected { provided, source })
    }
}

fn file_name(path: &Path) -> String {
    match path.file_name() {
        Some(name) => name.to_string_lossy().into_owned(),
        None => String::new(),
    }
}

fn is_path(text: &str) -> bool {
    text.contains(['/', '\\']) || Path::new(text).is_file()
}

fn environment_directory() -> PathBuf {
    if let Some(dir) = std::env::var_os("MOLFRAME_DATA_DIR") {
        return PathBuf::from(dir);
    }
    if let Some(dir) = std::env::var_os("XDG_CACHE_HOME") {
        return PathBuf::from(dir).join("molframe");
    }
    if cfg!(windows)
        && let Some(dir) = std::env::var_os("LOCALAPPDATA")
    {
        return PathBuf::from(dir).join("molframe");
    }
    match std::env::var_os("HOME") {
        Some(home) => PathBuf::from(home).join(".cache").join("molframe"),
        None => PathBuf::from(".molframe"),
    }
}

#[cfg(test)]
#[path = "table_store_tests.rs"]
mod tests;
