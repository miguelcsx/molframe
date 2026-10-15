//! Where tables are looked for, and which one a computation gets.
//!
//! A [`TableStore`] owns the cache directory, the profile and the tables held
//! in memory, and answers one question: given a request for a table of some
//! kind, which table is it? The request is what the caller passed
//! (`--geometry-table`, a Python argument); the profile is a plain text of
//! `kind: name` lines that says which entry a kind uses when the caller passed
//! nothing. A request that names a path is read from that path; a request that
//! names an entry, or a profile line, resolves the entry `<name>.<extension>`;
//! no request and no profile line is a typed [`TableError::Missing`], never a
//! guess.
//!
//! Three layers answer, in this order, both for which entry the profile names
//! and for where that entry's bytes are:
//!
//! 1. memory the caller registered ([`TableStore::insert_profile`],
//!    [`TableStore::insert_table`], [`TableStore::insert_entry`]);
//! 2. the disk: the profile file and `<directory>/<name>.<extension>`;
//! 3. the tables compiled into the library, present only with the
//!    `embedded-tables` feature: the bonding reach table and a profile line
//!    naming it.
//!
//! A layer that names nothing for a kind defers to the next. A target without
//! a filesystem (wasm) reads the disk as empty, so memory and the compiled-in
//! tables give it the same answer a native process gets from the same bytes.
//!
//! Decoded tables are memoised: a file per (kind, path, size, modification
//! time), so an edited file is read again; an entry in memory per kind, until
//! it is replaced. Either way a table is parsed once however many computations
//! ask for it, and a lookup after the first is a few hash probes under a read
//! lock and an `Arc` clone. Resolving from memory allocates nothing; the disk
//! layer costs the filesystem calls that notice an edit.
//!
//! Directory resolution order: `MOLFRAME_DATA_DIR`, then `XDG_CACHE_HOME/molframe`,
//! then `~/.cache/molframe` (`%LOCALAPPDATA%\molframe` on Windows). The profile is
//! `MOLFRAME_PROFILE` or `<directory>/profile.txt`.

use crate::table_bundled::bundled;
use crate::table_memory::{MemoryTables, profile_lines, read, write};
use crate::{Cache, UserFile};
use crate::{Provided, Source, TableBytes, TableError, TableKind, TableProvider};
use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};
use std::sync::{Arc, OnceLock, RwLock};
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
type Memo = HashMap<TypeId, HashMap<PathBuf, (Stamp, Arc<dyn Any + Send + Sync>)>>;

/// The cache directory, the profile, the tables in memory and the memo of
/// decoded files.
pub struct TableStore {
    directory: PathBuf,
    profile: PathBuf,
    memo: RwLock<Memo>,
    memory: RwLock<MemoryTables>,
    bundled: Option<Arc<MemoryTables>>,
}

impl core::fmt::Debug for TableStore {
    fn fmt(&self, out: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        out.debug_struct("TableStore")
            .field("directory", &self.directory)
            .field("profile", &self.profile)
            .field("embedded", &self.bundled.is_some())
            .finish_non_exhaustive()
    }
}

impl TableStore {
    /// A store over `directory`, with the profile at `directory/profile.txt`,
    /// nothing registered in memory, and the compiled-in tables when the
    /// `embedded-tables` feature is enabled.
    #[must_use]
    pub fn new(directory: PathBuf) -> Self {
        let profile = directory.join("profile.txt");
        Self {
            directory,
            profile,
            memo: RwLock::new(HashMap::new()),
            memory: RwLock::new(MemoryTables::new("memory")),
            bundled: bundled(),
        }
    }

    /// The process-wide store over the environment's directory and profile.
    ///
    /// Tables registered on it are seen by every later read in the process,
    /// including the bond perception a file read runs.
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

    /// Replaces the compiled-in layer, so a test can stand one in without the
    /// feature.
    #[cfg(test)]
    pub(crate) fn with_bundled(mut self, tables: Option<Arc<MemoryTables>>) -> Self {
        self.bundled = tables;
        self
    }

    /// Whether the tables compiled in by the `embedded-tables` feature back
    /// this store.
    #[must_use]
    pub const fn has_embedded_tables(&self) -> bool {
        self.bundled.is_some()
    }

    /// Registers a profile in memory, replacing any registered before it.
    ///
    /// `text` has the profile file's format, `kind: name` lines with `#`
    /// comments. A kind it names wins over the profile file and the compiled-in
    /// profile; a kind it does not name is still answered by them.
    pub fn insert_profile(&self, text: &str) {
        write(&self.memory).set_profile(text);
    }

    /// Registers the bytes of the entry `name` of kind `T`, as the cache file
    /// `<name>.<T::EXTENSION>` would hold them.
    ///
    /// See [`TableStore::insert_entry`].
    pub fn insert_table<T: TableKind>(&self, name: &str, bytes: impl Into<TableBytes>) {
        self.insert_entry(name, T::EXTENSION, bytes);
    }

    /// Registers the bytes of the entry `<name>.<extension>` in memory.
    ///
    /// The entry wins over a cache file and a compiled-in table of the same name
    /// and replaces an entry registered before it. Nothing is decoded here: the
    /// bytes are decoded once per kind on first selection, so a malformed table
    /// is reported where it is used, as a file's would be.
    pub fn insert_entry(&self, name: &str, extension: &str, bytes: impl Into<TableBytes>) {
        write(&self.memory).insert(name, extension, bytes.into());
    }

    /// The table `request` selects, or the profile's choice for `T` when it is `None`.
    ///
    /// A request containing a path separator, or naming an existing file, is a
    /// path; anything else is an entry name.
    ///
    /// # Errors
    ///
    /// [`TableError::Missing`] when nothing selects a table or no layer holds the
    /// entry; [`TableError::Unreadable`] or [`TableError::Malformed`] for a table
    /// that cannot be read or fails validation.
    pub fn select<T>(
        &self,
        request: Option<&str>,
        option: &'static str,
    ) -> Result<Selected<T>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        match request {
            Some(text) if is_path(text) => self.load::<T>(PathBuf::from(text)),
            Some(name) => self.entry::<T>(name, option),
            None => match self.profile_entry(T::KIND)? {
                Some(name) => self.entry::<T>(&name, option),
                None => Err(TableError::Missing {
                    kind: T::KIND,
                    wanted: self.unnamed(T::KIND),
                    option,
                }),
            },
        }
    }

    /// The profile's choice for `T`, or `None` when no profile names one.
    ///
    /// For a kind that is optional (a component dictionary layered under in-file
    /// definitions, the reach that bond perception needs) rather than required.
    /// An entry a profile names but no layer holds is still an error.
    ///
    /// # Errors
    ///
    /// As [`TableStore::select`], except that an absent profile or an absent line is
    /// `Ok(None)`.
    pub fn select_optional<T>(
        &self,
        option: &'static str,
    ) -> Result<Option<Selected<T>>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        match self.profile_entry(T::KIND)? {
            Some(name) => self.entry::<T>(&name, option).map(Some),
            None => Ok(None),
        }
    }

    /// The entry the first layer that names one names for `kind`.
    fn profile_entry(&self, kind: &str) -> Result<Option<Arc<str>>, TableError> {
        if let Some(name) = read(&self.memory).profile_entry(kind) {
            return Ok(Some(name));
        }
        if let Some(text) = self.profile_text()?
            && let Some((_, name)) = profile_lines(&text).find(|(key, _)| *key == kind)
        {
            return Ok(Some(name.into()));
        }
        Ok(match &self.bundled {
            Some(tables) => tables.profile_entry(kind),
            None => None,
        })
    }

    /// The profile file's text, or `None` when there is no file to read.
    fn profile_text(&self) -> Result<Option<String>, TableError> {
        match std::fs::read_to_string(&self.profile) {
            Ok(text) => Ok(Some(text)),
            Err(error) if absent(&error) => Ok(None),
            Err(error) => Err(TableError::Unreadable {
                path: self.profile.display().to_string(),
                reason: error.to_string(),
            }),
        }
    }

    /// Why no profile names `kind`, for the error that reports it.
    fn unnamed(&self, kind: &str) -> String {
        match self.profile_text() {
            Ok(Some(_)) => format!(
                "no request was made and {} has no `{kind}:` line",
                self.profile.display()
            ),
            _ => format!(
                "no request was made and there is no profile at {}",
                self.profile.display()
            ),
        }
    }

    /// The entry `name` of kind `T`, from the first layer that holds it.
    fn entry<T>(&self, name: &str, option: &'static str) -> Result<Selected<T>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        let registered = read(&self.memory).entry(name, T::EXTENSION);
        if let Some(entry) = registered {
            return from_memory(&entry);
        }
        let path = self.directory.join(format!("{name}.{}", T::EXTENSION));
        match std::fs::metadata(&path) {
            Ok(meta) => return self.stamped::<T>(path, &meta, Source::Cache),
            Err(error) if absent(&error) => {}
            Err(error) => {
                return Err(TableError::Unreadable {
                    path: path.display().to_string(),
                    reason: error.to_string(),
                });
            }
        }
        if let Some(tables) = &self.bundled
            && let Some(entry) = tables.entry(name, T::EXTENSION)
        {
            return from_memory(&entry);
        }
        Err(TableError::Missing {
            kind: T::KIND,
            wanted: format!("the cache holds no entry at {}", path.display()),
            option,
        })
    }

    /// A table file the caller named by path; its absence is an error, not a
    /// fall-through to another layer.
    fn load<T>(&self, path: PathBuf) -> Result<Selected<T>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        match std::fs::metadata(&path) {
            Ok(meta) => self.stamped::<T>(path, &meta, Source::User),
            Err(error) => Err(TableError::Unreadable {
                path: path.display().to_string(),
                reason: error.to_string(),
            }),
        }
    }

    /// The file at `path`, decoded once per size and modification time.
    fn stamped<T>(
        &self,
        path: PathBuf,
        meta: &std::fs::Metadata,
        source: Source,
    ) -> Result<Selected<T>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        let stamp = (meta.len(), meta.modified().ok());
        let kind = TypeId::of::<T>();
        if let Some(files) = read(&self.memo).get(&kind)
            && let Some((seen, any)) = files.get(&path)
            && *seen == stamp
            && let Ok(provided) = Arc::clone(any).downcast::<Provided<T>>()
        {
            return Ok(Selected { provided, source });
        }
        let provided = match source {
            Source::Cache => Cache::<T>::new(&self.directory, &file_name(&path)).shared()?,
            _ => UserFile::<T>::new(path.clone()).shared()?,
        };
        write(&self.memo).entry(kind).or_default().insert(
            path,
            (stamp, Arc::clone(&provided) as Arc<dyn Any + Send + Sync>),
        );
        Ok(Selected { provided, source })
    }
}

fn from_memory<T>(entry: &crate::table_memory::MemoryEntry) -> Result<Selected<T>, TableError>
where
    T: TableKind + Send + Sync + 'static,
{
    Ok(Selected {
        provided: entry.provided::<T>()?,
        source: Source::Memory,
    })
}

/// Whether a filesystem error means the file is not there. A target without a
/// filesystem (wasm) reports `Unsupported`, and a cache path through a plain
/// file reports `NotADirectory`; both are the same absence.
fn absent(error: &std::io::Error) -> bool {
    matches!(
        error.kind(),
        ErrorKind::NotFound | ErrorKind::Unsupported | ErrorKind::NotADirectory
    )
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
