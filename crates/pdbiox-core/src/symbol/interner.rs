//! The identifier dictionary itself.
//!
//! Local strings live in one arena rather than in one allocation each. A file
//! with a hundred thousand atoms has a few hundred distinct identifiers, so
//! interning them costs a few hundred bytes copied once, not an allocation per
//! atom.

use super::canonical;
use hashbrown::HashTable;
use hashbrown::hash_table::Entry;
use std::convert::TryFrom;
use std::fmt;
use std::hash::{BuildHasher, RandomState};
use std::sync::Arc;

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
    /// Largest identifier emitted by an [`Interner`].
    ///
    /// The final `u32` value is reserved because [`crate::AltId`] encodes a
    /// labelled symbol by adding one while retaining zero for the blank label.
    const MAX_VALID_RAW: u32 = u32::MAX - 1;

    /// Returns true when this identifier has the same meaning in every
    /// structure, so comparing it across structures is sound.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn is_canonical(self) -> bool {
        canonical::text_of(self.0).is_some()
    }

    /// The raw identifier, for storage in a column.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub const fn get(self) -> u32 {
        self.0
    }

    /// Rebuilds an identifier from a stored value.
    ///
    /// `u32::MAX` is reserved and is not emitted by [`Interner`]. Callers
    /// should pass values previously produced by this crate.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub const fn from_raw(raw: u32) -> Self {
        Self(raw)
    }
}

impl fmt::Debug for SymbolId {
    /// Formats canonical identifiers by text and local identifiers by raw value.
    ///
    /// Runs in `O(L)` formatting time for canonical text of length `L` and
    /// allocates no intermediate heap storage.
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
///
/// Storing the exclusive end offset avoids repeated addition during resolution
/// without increasing the eight-byte extent representation.
#[derive(Clone, Copy)]
struct Extent {
    start: u32,
    end: u32,
}

/// A type-safe index into the local extent array.
///
/// The transparent representation keeps each hash-table entry at four bytes
/// while preventing local ordinals from being confused with global IDs.
#[derive(Clone, Copy)]
#[repr(transparent)]
struct LocalOrdinal(u32);

/// Copy-on-write storage for the local dictionary.
#[derive(Clone)]
struct Storage {
    arena: String,
    extents: Vec<Extent>,
    table: HashTable<LocalOrdinal>,
}

/// A structure's identifier dictionary.
///
/// Interning is idempotent: the same string always yields the same identifier
/// within one interner, and identifiers are issued in first-seen order, so two
/// runs over the same file produce the same numbering.
#[derive(Clone)]
pub struct Interner {
    storage: Arc<Storage>,
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
    ///
    /// Runs in `O(1)` time and creates no local-string allocation.
    #[must_use]
    pub fn new() -> Self {
        Self {
            storage: Arc::new(Storage {
                arena: String::new(),
                extents: Vec::new(),
                table: HashTable::new(),
            }),
            hasher: RandomState::new(),
            limit: Self::DEFAULT_LIMIT,
        }
    }

    /// Creates a dictionary with room reserved for `identifiers` local strings.
    ///
    /// The extent vector and hash table reserve capacity up front. The arena is
    /// left empty because the total byte length cannot be inferred from the
    /// number of identifiers.
    ///
    /// Uses `O(identifiers)` reserved memory.
    #[must_use]
    pub fn with_capacity(identifiers: usize) -> Self {
        let mut interner = Self::new();
        let storage = Arc::make_mut(&mut interner.storage);
        storage.extents.reserve(identifiers);
        storage.table.reserve(identifiers, |_| 0u64);
        interner
    }

    /// Sets the limit on local identifiers.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    pub const fn with_limit(mut self, limit: u32) -> Self {
        self.limit = limit;
        self
    }

    /// Returns the identifier for `text`, interning it if it is new.
    ///
    /// Existing strings require no allocation. A new local string is copied
    /// exactly once into the arena and adds one extent and one table entry.
    ///
    /// Expected running time is `O(L)`, where `L` is `text.len()`. Insertion is
    /// amortized over hash-table and arena growth.
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
            storage,
            hasher,
            limit,
        } = self;
        let Storage {
            arena,
            extents,
            table,
        } = Arc::make_mut(storage);
        let hash = hasher.hash_one(text);

        let pending = pending_local(extents.len(), *limit, arena.len(), text.len());

        let (local, id, extent) = match pending {
            Ok(pending) => pending,
            Err(error) => {
                return match find_local(table, hash, arena, extents, text) {
                    Some(local) => local_to_id(local).ok_or(error),
                    None => Err(error),
                };
            }
        };

        match table.entry(
            hash,
            |&other| extent_text(arena, extents, other) == Some(text),
            |&other| {
                // Rehashing during growth reads the arena, which the closure may
                // borrow because only the table is being mutated here. An ordinal
                // with no extent cannot occur; hashing the empty string keeps the
                // table consistent rather than aborting if it somehow did.
                match extent_text(arena, extents, other) {
                    Some(existing) => hasher.hash_one(existing),
                    None => hasher.hash_one(""),
                }
            },
        ) {
            Entry::Occupied(occupied) => local_to_id(*occupied.get()).ok_or(DictionaryFull),
            Entry::Vacant(vacant) => {
                arena.push_str(text);
                extents.push(extent);
                let _ = vacant.insert(local);
                Ok(id)
            }
        }
    }

    /// Returns the identifier for `text` if it is already known, without
    /// interning it.
    ///
    /// This is what a selection uses: a name that no structure contains should
    /// resolve to nothing rather than growing the dictionary.
    ///
    /// Expected running time is `O(L)`, where `L` is `text.len()`. The operation
    /// performs no allocation and does not mutate the dictionary.
    #[must_use]
    pub fn get(&self, text: &str) -> Option<SymbolId> {
        if let Some(ordinal) = canonical::ordinal_of(text) {
            return Some(SymbolId(ordinal));
        }

        let hash = self.hasher.hash_one(text);
        find_local(
            &self.storage.table,
            hash,
            &self.storage.arena,
            &self.storage.extents,
            text,
        )
        .and_then(local_to_id)
    }

    /// Returns the string an identifier names.
    ///
    /// Runs in `O(1)` time and allocates no memory. The returned slice borrows
    /// either static canonical storage or this interner's arena.
    #[must_use]
    pub fn resolve(&self, id: SymbolId) -> Option<&str> {
        if let Some(text) = canonical::text_of(id.0) {
            return Some(text);
        }

        let local = id_to_local(id)?;
        extent_text(&self.storage.arena, &self.storage.extents, local)
    }

    /// The number of local identifiers issued.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.storage.extents.len()
    }

    /// Returns true when no local identifier has been issued.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.storage.extents.is_empty()
    }

    /// The bytes the local strings occupy.
    ///
    /// Runs in `O(1)` time and allocates no memory.
    #[must_use]
    #[inline]
    pub fn arena_len(&self) -> usize {
        self.storage.arena.len()
    }

    /// Walks canonical and structure-local identifiers in stable id order.
    ///
    /// This is the dictionary-side half of compiling glob predicates: pattern
    /// work is proportional to distinct identifiers and the atom loop remains
    /// integer membership. Iteration allocates nothing.
    pub fn iter(&self) -> impl Iterator<Item = (SymbolId, &str)> {
        let canonical = canonical::CANONICAL
            .iter()
            .enumerate()
            .filter_map(|(position, text)| {
                u32::try_from(position)
                    .ok()
                    .map(|position| (SymbolId(position), *text))
            });
        let local = self
            .storage
            .extents
            .iter()
            .enumerate()
            .filter_map(|(position, _extent)| {
                let local = u32::try_from(position).ok().map(LocalOrdinal)?;
                let symbol = local_to_id(local)?;
                extent_text(&self.storage.arena, &self.storage.extents, local)
                    .map(|text| (symbol, text))
            });
        canonical.chain(local)
    }
}

impl Default for Interner {
    /// Creates an empty interner with [`Interner::DEFAULT_LIMIT`].
    ///
    /// Runs in `O(1)` time and creates no local-string allocation.
    fn default() -> Self {
        Self::new()
    }
}

impl fmt::Debug for Interner {
    /// Formats the interner's local-count and arena-size metadata.
    ///
    /// Runs in `O(1)` time and does not traverse or copy interned strings.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Interner")
            .field("local", &self.storage.extents.len())
            .field("arena_bytes", &self.storage.arena.len())
            .finish_non_exhaustive()
    }
}

/// Prepares the metadata required to insert one new local string.
///
/// Returns the next local ordinal, its global [`SymbolId`], and its arena
/// extent. No state is mutated, so failure cannot leave a partial insertion.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn pending_local(
    local_count: usize,
    limit: u32,
    arena_len: usize,
    text_len: usize,
) -> Result<(LocalOrdinal, SymbolId, Extent), DictionaryFull> {
    let local = match u32::try_from(local_count) {
        Ok(local) if local < limit => LocalOrdinal(local),
        Ok(_) | Err(_) => return Err(DictionaryFull),
    };

    let id = local_to_id(local).ok_or(DictionaryFull)?;
    let extent = extent_for_append(arena_len, text_len).ok_or(DictionaryFull)?;

    Ok((local, id, extent))
}

/// Computes the arena extent for appending a string.
///
/// Returns `None` if either offset cannot be represented by the compact `u32`
/// extent format or if the byte-length addition overflows.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn extent_for_append(arena_len: usize, text_len: usize) -> Option<Extent> {
    let end = arena_len.checked_add(text_len)?;

    Some(Extent {
        start: u32::try_from(arena_len).ok()?,
        end: u32::try_from(end).ok()?,
    })
}

/// Finds a local ordinal matching `text`.
///
/// Expected running time is `O(L)`, where `L` is `text.len()`. The operation
/// performs no allocation.
#[inline]
fn find_local(
    table: &HashTable<LocalOrdinal>,
    hash: u64,
    arena: &str,
    extents: &[Extent],
    text: &str,
) -> Option<LocalOrdinal> {
    table
        .find(hash, |&local| {
            extent_text(arena, extents, local) == Some(text)
        })
        .copied()
}

/// Converts a local ordinal into its globally stored identifier.
///
/// Returns `None` if the canonical offset would overflow or produce the
/// reserved `u32::MAX` value.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn local_to_id(local: LocalOrdinal) -> Option<SymbolId> {
    let canonical_count = u32::try_from(canonical::CANONICAL.len()).ok()?;
    let raw = canonical_count.checked_add(local.0)?;

    if raw > SymbolId::MAX_VALID_RAW {
        None
    } else {
        Some(SymbolId(raw))
    }
}

/// Converts a non-canonical identifier into its local ordinal.
///
/// Returns `None` for canonical identifiers, the reserved maximum value, or
/// identifiers that cannot be represented by the current dictionary layout.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn id_to_local(id: SymbolId) -> Option<LocalOrdinal> {
    if id.0 > SymbolId::MAX_VALID_RAW {
        return None;
    }

    let canonical_count = u32::try_from(canonical::CANONICAL.len()).ok()?;
    id.0.checked_sub(canonical_count).map(LocalOrdinal)
}

/// Returns the arena slice represented by a local ordinal.
///
/// Invalid ordinals, invalid UTF-8 boundaries, and malformed extents produce
/// `None` rather than panicking.
///
/// Runs in `O(1)` time and allocates no memory.
#[inline]
fn extent_text<'a>(arena: &'a str, extents: &[Extent], local: LocalOrdinal) -> Option<&'a str> {
    let local_index = usize::try_from(local.0).ok()?;
    let extent = extents.get(local_index)?;
    let start = usize::try_from(extent.start).ok()?;
    let end = usize::try_from(extent.end).ok()?;

    arena.get(start..end)
}

#[cfg(test)]
#[path = "interner_tests.rs"]
mod tests;
