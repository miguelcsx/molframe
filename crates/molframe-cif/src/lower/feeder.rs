//! The destination a reader pushes borrowed `atom_site` rows into.
//!
//! Readers feed millions of rows, so the row type stays a generic parameter of
//! [`RowFeeder::feed`]: the lowering code is compiled once per reader rather
//! than reached through a trait object on every field access. The feeder
//! itself picks the destination with one branch per row.

use super::atoms::{AtomBuilder, AtomSiteRow};
use super::ragged::RaggedBuilder;
use std::fmt;

enum Target<'a, 'b> {
    Atom(&'a mut AtomBuilder<'b>),
    Ragged(&'a mut RaggedBuilder<'b>),
}

/// Synchronous destination for borrowed `atom_site` rows.
///
/// Rows must be fed in deposition order. The lowerer consumes the row and every
/// borrow it returns before [`RowFeeder::feed`] returns.
#[doc(hidden)]
pub struct RowFeeder<'a, 'b> {
    target: Target<'a, 'b>,
}

impl<'a, 'b> RowFeeder<'a, 'b> {
    pub(super) fn atom(builder: &'a mut AtomBuilder<'b>) -> Self {
        Self {
            target: Target::Atom(builder),
        }
    }

    pub(super) fn ragged(builder: &'a mut RaggedBuilder<'b>) -> Self {
        Self {
            target: Target::Ragged(builder),
        }
    }

    /// Lowers one complete row.
    pub fn feed<R: AtomSiteRow + ?Sized>(&mut self, row: &R) {
        match &mut self.target {
            Target::Atom(builder) => builder.feed(row),
            Target::Ragged(builder) => builder.feed(row),
        }
    }
}

impl fmt::Debug for RowFeeder<'_, '_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.debug_struct("RowFeeder").finish_non_exhaustive()
    }
}
