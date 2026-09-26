//! Producer-side memory of compiled queries and of one pending result.
//!
//! A consumer extension cannot hold a Rust `Query` — the capsule carries no
//! Rust values across the boundary — so every selection arrives as text. An
//! interactive consumer re-sends the same handful of queries over and over; a
//! small cache keyed by that text means each distinct query is lexed, parsed
//! and planned once per structure snapshot rather than once per call. The
//! snapshot a capsule retains is immutable, and every capsule selection uses
//! the default policy, so the query text alone is a complete key.
//!
//! The cache is bounded: at most [`QueryCache::CAPACITY`] queries are kept, and
//! the least recently used one is evicted first. Lookup is a linear scan over
//! that small fixed bound.

use std::sync::Arc;

/// Compiled queries keyed by their exact source text.
#[derive(Debug, Default)]
pub(crate) struct QueryCache {
    /// Most recently used last.
    entries: Vec<(Box<str>, Arc<molframe::Query>)>,
}

impl QueryCache {
    /// How many distinct queries one structure snapshot remembers.
    pub(crate) const CAPACITY: usize = 64;

    /// The compiled query for `source`, compiling it on first use.
    pub(crate) fn get_or_compile(
        &mut self,
        source: &str,
    ) -> Result<Arc<molframe::Query>, molframe::Findings> {
        if let Some(index) = self
            .entries
            .iter()
            .position(|(known, _)| known.as_ref() == source)
        {
            let entry = self.entries.remove(index);
            let query = Arc::clone(&entry.1);
            self.entries.push(entry);
            return Ok(query);
        }
        let query = Arc::new(molframe::Query::compile(source)?);
        if self.entries.len() >= Self::CAPACITY {
            let _ = self.entries.remove(0);
        }
        self.entries.push((source.into(), Arc::clone(&query)));
        Ok(query)
    }

    /// Number of remembered queries.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }
}

/// One evaluated selection waiting to be copied out.
///
/// The capsule's selection callback answers a size request by evaluating the
/// query and keeping the rows here; the immediately following copy request for
/// the same text takes them instead of evaluating again. This lets a consumer
/// allocate exactly as many rows as were selected rather than one per atom.
#[derive(Debug, Default)]
pub(crate) struct PendingRows {
    rows: Option<(Box<str>, molframe::query::Evaluation)>,
}

impl PendingRows {
    /// Keeps `selection` for the next copy of `source`.
    pub(crate) fn keep(&mut self, source: &str, selection: molframe::query::Evaluation) {
        self.rows = Some((source.into(), selection));
    }

    /// Takes the kept selection when it answers `source`.
    pub(crate) fn take(&mut self, source: &str) -> Option<molframe::query::Evaluation> {
        match self.rows.take() {
            Some((known, selection)) if known.as_ref() == source => Some(selection),
            _ => None,
        }
    }
}

#[cfg(test)]
#[path = "query_cache_tests.rs"]
mod tests;
