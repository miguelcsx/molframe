//! The identifier dictionary itself.
//!
//! Local strings live in one arena rather than in one allocation each. A file
//! with a hundred thousand atoms has a few hundred distinct identifiers, so
//! interning them costs a few hundred bytes copied once, not an allocation per
//! atom.

use super::canonical;
use hashbrown::HashTable;
use std::fmt;
use std::hash::{BuildHasher, RandomState};

/// An interned identifier.
///
/// Identifiers below the canonical count mean the same string in every
/// structure. Above it they are meaningful only against the interner that
/// issued them.
///
/// # Examples
///
/// ```
/// use pdbiox_core::{Interner, SymbolId};
///
/// let mut interner = Interner::new();
/// let alanine = interner.intern("ALA")?;
/// assert!(alanine.is_canonical());
///
/// let ligand = interner.intern("MY_LIGAND")?;
/// assert!(!ligand.is_canonical());
/// assert_eq!(interner.resolve(ligand), Some("MY_LIGAND"));
/// # Ok::<(), pdbiox_core::DictionaryFull>(())
/// ```
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(transparent)]
pub struct SymbolId(u32);

impl SymbolId {
    /// Returns true when this identifier has the same meaning in every
    /// structure, so comparing it across structures is sound.
    #[must_use]
    pub fn is_canonical(self) -> bool {
        (self.0 as usize) < canonical::CANONICAL.len()
    }

    /// The raw identifier, for storage in a column.
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0
    }

    /// Rebuilds an identifier from a stored value.
    #[must_use]
    pub const fn from_raw(raw: u32) -> Self {
        Self(raw)
    }
}

impl fmt::Debug for SymbolId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match canonical::text_of(self.0) {
            Some(text) => write!(f, "SymbolId({text:?})"),
            None => write!(f, "SymbolId(local {})", self.0),
        }
    }
}

/// The identifier dictionary grew past the limit that guards against a hostile
/// input inventing identifiers without bound.
#[derive(Clone, Copy, PartialEq, Eq, Debug, thiserror::Error)]
#[error("identifier dictionary is full")]
pub struct DictionaryFull;

/// Where a local string sits in the arena.
#[derive(Clone, Copy)]
struct Extent {
    start: u32,
    len: u32,
}

/// A structure's identifier dictionary.
///
/// Interning is idempotent: the same string always yields the same identifier
/// within one interner, and identifiers are issued in first-seen order, so two
/// runs over the same file produce the same numbering.
#[derive(Clone)]
pub struct Interner {
    arena: String,
    extents: Vec<Extent>,
    table: HashTable<u32>,
    hasher: RandomState,
    limit: u32,
}

impl Interner {
    /// The largest number of local identifiers one structure may issue.
    ///
    /// This is a guard, not a design target. A real structure uses a few
    /// hundred; reaching this many means the input is generating identifiers
    /// rather than naming things.
    pub const DEFAULT_LIMIT: u32 = 1_000_000;

    /// Creates an empty dictionary over the canonical one.
    #[must_use]
    pub fn new() -> Self {
        Self {
            arena: String::new(),
            extents: Vec::new(),
            table: HashTable::new(),
            hasher: RandomState::new(),
            limit: Self::DEFAULT_LIMIT,
        }
    }

    /// Creates a dictionary with room reserved for `identifiers` local strings.
    #[must_use]
    pub fn with_capacity(identifiers: usize) -> Self {
        let mut interner = Self::new();
        interner.extents.reserve(identifiers);
        interner.table.reserve(identifiers, |_| 0u64);
        interner
    }

    /// Sets the limit on local identifiers.
    #[must_use]
    pub const fn with_limit(mut self, limit: u32) -> Self {
        self.limit = limit;
        self
    }

    /// Returns the identifier for `text`, interning it if it is new.
    ///
    /// # Errors
    ///
    /// Returns [`DictionaryFull`] when the local dictionary has reached its
    /// limit and `text` is not already in it.
    pub fn intern(&mut self, text: &str) -> Result<SymbolId, DictionaryFull> {
        if let Some(ordinal) = canonical::ordinal_of(text) {
            return Ok(SymbolId(ordinal));
        }
        let Self {
            arena,
            extents,
            table,
            hasher,
            limit,
        } = self;
        let hash = hasher.hash_one(text);
        let equal = |&local: &u32| extent_text(arena, extents, local) == Some(text);

        if let Some(&local) = table.find(hash, equal) {
            return Ok(SymbolId(local_to_id(local)));
        }
        if extents.len() as u32 >= *limit {
            return Err(DictionaryFull);
        }

        let start = arena.len() as u32;
        arena.push_str(text);
        let local = extents.len() as u32;
        extents.push(Extent {
            start,
            len: text.len() as u32,
        });
        table.insert_unique(hash, local, |&other| {
            // Rehashing during growth reads the arena, which the closure may
            // borrow because only the table is being mutated here. An ordinal
            // with no extent cannot occur; hashing the empty string keeps the
            // table consistent rather than aborting if it somehow did.
            match extent_text(arena, extents, other) {
                Some(text) => hasher.hash_one(text),
                None => hasher.hash_one(""),
            }
        });
        Ok(SymbolId(local_to_id(local)))
    }

    /// Returns the identifier for `text` if it is already known, without
    /// interning it.
    ///
    /// This is what a selection uses: a name that no structure contains should
    /// resolve to nothing rather than growing the dictionary.
    #[must_use]
    pub fn get(&self, text: &str) -> Option<SymbolId> {
        if let Some(ordinal) = canonical::ordinal_of(text) {
            return Some(SymbolId(ordinal));
        }
        let hash = self.hasher.hash_one(text);
        self.table
            .find(hash, |&local| {
                extent_text(&self.arena, &self.extents, local) == Some(text)
            })
            .map(|&local| SymbolId(local_to_id(local)))
    }

    /// Returns the string an identifier names.
    #[must_use]
    pub fn resolve(&self, id: SymbolId) -> Option<&str> {
        if let Some(text) = canonical::text_of(id.0) {
            Some(text)
        } else {
            let local = id.0.checked_sub(canonical::CANONICAL.len() as u32)?;
            extent_text(&self.arena, &self.extents, local)
        }
    }

    /// The number of local identifiers issued.
    #[must_use]
    pub fn len(&self) -> usize {
        self.extents.len()
    }

    /// Returns true when no local identifier has been issued.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.extents.is_empty()
    }

    /// The bytes the local strings occupy.
    #[must_use]
    pub fn arena_len(&self) -> usize {
        self.arena.len()
    }
}

impl Default for Interner {
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Debug for Interner {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Interner")
            .field("local", &self.extents.len())
            .field("arena_bytes", &self.arena.len())
            .finish_non_exhaustive()
    }
}

fn local_to_id(local: u32) -> u32 {
    local.saturating_add(canonical::CANONICAL.len() as u32)
}

fn extent_text<'a>(arena: &'a str, extents: &[Extent], local: u32) -> Option<&'a str> {
    let extent = extents.get(local as usize)?;
    let start = extent.start as usize;
    arena.get(start..start + extent.len as usize)
}

#[cfg(test)]
#[path = "interner_tests.rs"]
mod tests;
