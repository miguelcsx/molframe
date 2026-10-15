//! Tables and a profile held in memory, under the names a cache directory uses.
//!
//! A [`MemoryTables`] is one layer of a [`TableStore`](crate::TableStore): a
//! profile (which entry each kind uses) and a set of entries keyed exactly as
//! cache files are, by name and extension. A target without a filesystem gets
//! its tables here, and so does a host that already holds the bytes.
//!
//! Cost: registering parses a profile once, linear in its length. A lookup is
//! two hash probes and an `Arc` clone with no allocation. Each entry decodes
//! once per table kind, through [`MemoryFile`]; every later selection of that
//! kind shares the same decoded table.

use crate::{MemoryFile, Provided, TableBytes, TableError, TableKind, TableProvider};
use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::sync::{Arc, RwLock, RwLockReadGuard, RwLockWriteGuard};

type Decoded = Result<Arc<dyn Any + Send + Sync>, TableError>;

/// The lock for reading. A writer that panicked left the data whole — every
/// write here is one insertion — so a poisoned lock is still read.
pub(crate) fn read<T>(lock: &RwLock<T>) -> RwLockReadGuard<'_, T> {
    match lock.read() {
        Ok(guard) => guard,
        Err(poisoned) => poisoned.into_inner(),
    }
}

/// The lock for writing, on the same terms as [`read`].
pub(crate) fn write<T>(lock: &RwLock<T>) -> RwLockWriteGuard<'_, T> {
    match lock.write() {
        Ok(guard) => guard,
        Err(poisoned) => poisoned.into_inner(),
    }
}

/// The `kind: name` lines of a profile text, in order, comments skipped.
///
/// The one reading of the profile format, shared by the file on disk and the
/// profile held in memory.
pub(crate) fn profile_lines(text: &str) -> impl Iterator<Item = (&str, &str)> {
    text.lines().filter_map(|line| {
        let line = line.trim();
        if line.starts_with('#') {
            return None;
        }
        let (key, value) = line.split_once(':')?;
        let value = value.trim();
        (!value.is_empty()).then_some((key.trim(), value))
    })
}

/// One table's bytes and the tables decoded from them, one per kind.
pub(crate) struct MemoryEntry {
    label: String,
    bytes: TableBytes,
    decoded: RwLock<Vec<(TypeId, Decoded)>>,
}

impl MemoryEntry {
    fn new(label: String, bytes: TableBytes) -> Self {
        Self {
            label,
            bytes,
            decoded: RwLock::new(Vec::new()),
        }
    }

    /// The entry decoded as `T`, decoding on the first request for that kind.
    ///
    /// The outcome is kept either way, so malformed bytes are reported again
    /// without being decoded again.
    pub(crate) fn provided<T>(&self) -> Result<Arc<Provided<T>>, TableError>
    where
        T: TableKind + Send + Sync + 'static,
    {
        let kind = TypeId::of::<T>();
        if let Some(found) = find::<T>(&read(&self.decoded), kind) {
            return found;
        }
        let mut decoded = write(&self.decoded);
        // Another thread may have decoded it between the two locks; theirs is
        // kept so every caller shares one table.
        if let Some(found) = find::<T>(&decoded, kind) {
            return found;
        }
        let outcome = MemoryFile::<T>::new(&self.label, self.bytes.clone()).shared();
        decoded.push((
            kind,
            match &outcome {
                Ok(provided) => Ok(Arc::clone(provided) as Arc<dyn Any + Send + Sync>),
                Err(error) => Err(error.clone()),
            },
        ));
        outcome
    }
}

fn find<T>(
    decoded: &[(TypeId, Decoded)],
    kind: TypeId,
) -> Option<Result<Arc<Provided<T>>, TableError>>
where
    T: Send + Sync + 'static,
{
    let (_, outcome) = decoded.iter().find(|(seen, _)| *seen == kind)?;
    Some(match outcome {
        Ok(any) => match Arc::clone(any).downcast::<Provided<T>>() {
            Ok(provided) => Ok(provided),
            Err(_) => return None,
        },
        Err(error) => Err(error.clone()),
    })
}

/// A profile and the entries it can name, all in memory.
#[derive(Default)]
pub(crate) struct MemoryTables {
    origin: &'static str,
    profile: HashMap<Box<str>, Arc<str>>,
    entries: HashMap<Box<str>, HashMap<Box<str>, Arc<MemoryEntry>>>,
}

impl MemoryTables {
    /// An empty layer whose entries are labelled `<origin>:<name>.<extension>`.
    pub(crate) fn new(origin: &'static str) -> Self {
        Self {
            origin,
            ..Self::default()
        }
    }

    /// Replaces the profile with the `kind: name` lines of `text`; the first
    /// line for a kind wins, as it does in a profile file.
    pub(crate) fn set_profile(&mut self, text: &str) {
        self.profile.clear();
        for (kind, name) in profile_lines(text) {
            self.profile
                .entry(kind.into())
                .or_insert_with(|| name.into());
        }
    }

    /// Holds `bytes` as the entry `<name>.<extension>`, replacing any before it.
    pub(crate) fn insert(&mut self, name: &str, extension: &str, bytes: TableBytes) {
        let label = format!("{}:{name}.{extension}", self.origin);
        let entry = Arc::new(MemoryEntry::new(label, bytes));
        match self.entries.get_mut(extension) {
            Some(names) => {
                names.insert(name.into(), entry);
            }
            None => {
                self.entries
                    .insert(extension.into(), HashMap::from([(name.into(), entry)]));
            }
        }
    }

    /// The entry the profile names for `kind`, if it names one.
    pub(crate) fn profile_entry(&self, kind: &str) -> Option<Arc<str>> {
        self.profile.get(kind).map(Arc::clone)
    }

    /// The entry `<name>.<extension>`, if it is held.
    pub(crate) fn entry(&self, name: &str, extension: &str) -> Option<Arc<MemoryEntry>> {
        self.entries.get(extension)?.get(name).map(Arc::clone)
    }
}

#[cfg(test)]
#[path = "table_memory_tests.rs"]
mod tests;
