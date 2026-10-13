//! Strings shared between the columns of one category.
//!
//! Columns of a category repeat each other's text: the atom names of a
//! component appear in its atom, alias and bond columns alike. Routing every
//! column through one interner stores each distinct text once, so a lookup is
//! `O(1)` expected and a repeated text costs a reference count, not a copy.

use std::collections::HashSet;
use std::sync::Arc;

/// Distinct shared strings seen so far.
#[derive(Debug, Default)]
pub(crate) struct Interner {
    seen: HashSet<Arc<str>>,
}

impl Interner {
    /// The shared copy of `value`, created on first sight.
    pub(crate) fn share(&mut self, value: &str) -> Arc<str> {
        if let Some(shared) = self.seen.get(value) {
            return Arc::clone(shared);
        }
        let shared: Arc<str> = Arc::from(value);
        self.seen.insert(Arc::clone(&shared));
        shared
    }
}
