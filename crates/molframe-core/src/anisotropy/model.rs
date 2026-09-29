//! Anisotropic displacement ellipsoids stored as a sparse per-atom table.
//!
//! Most structures carry no anisotropic ellipsoids, so the table is optional
//! in the same sense connectivity is: `available` records whether the source
//! resolved the category at all, distinct from an empty but known set. Rows
//! are kept sorted by atom position, so building costs O(rows log rows) once
//! and every per-atom lookup costs O(log rows); a structure without the
//! category pays nothing beyond a shared empty table.

use crate::index::{AnisotropyIndex, AtomIndex};
use std::sync::Arc;

/// One anisotropic displacement ellipsoid before it is packed into columns.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnisotropicDisplacement {
    /// The atom row this ellipsoid belongs to.
    pub atom: AtomIndex,
    /// Displacement tensor in ångström squared, stored as
    /// `[U11, U22, U33, U12, U13, U23]`.
    pub u: [f32; 6],
}

/// A compact, immutable anisotropic displacement table sorted by atom.
#[derive(Clone, Debug)]
pub struct AnisotropyTable {
    atom: Arc<Vec<AtomIndex>>,
    u: Arc<Vec<[f32; 6]>>,
    available: bool,
}

impl Default for AnisotropyTable {
    fn default() -> Self {
        Self {
            atom: Arc::new(Vec::new()),
            u: Arc::new(Vec::new()),
            available: false,
        }
    }
}

impl AnisotropyTable {
    /// Number of ellipsoid rows.
    #[must_use]
    pub fn len(&self) -> usize {
        self.atom.len()
    }

    /// Whether no ellipsoid is recorded.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atom.is_empty()
    }

    /// Whether anisotropic displacement was resolved, including a known empty set.
    #[must_use]
    pub const fn is_available(&self) -> bool {
        self.available
    }

    /// Returns one ellipsoid by position.
    #[must_use]
    pub fn get(&self, row: AnisotropyIndex) -> Option<AnisotropicDisplacement> {
        let position = row.as_usize();
        Some(AnisotropicDisplacement {
            atom: *self.atom.get(position)?,
            u: *self.u.get(position)?,
        })
    }

    /// Returns the tensor of one atom, or `None` when it is isotropic.
    #[must_use]
    pub fn for_atom(&self, atom: AtomIndex) -> Option<[f32; 6]> {
        let position = self.atom.binary_search(&atom).ok()?;
        Some(self.u[position])
    }

    /// Walks records in stable atom order.
    #[must_use]
    pub fn iter(&self) -> impl ExactSizeIterator<Item = AnisotropicDisplacement> + '_ {
        self.atom
            .iter()
            .zip(self.u.iter())
            .map(|(atom, u)| AnisotropicDisplacement { atom: *atom, u: *u })
    }
}

/// Deterministic anisotropy-table construction with atom deduplication.
///
/// Records accumulate in insertion order and are sorted and deduplicated once
/// at `finish`. The first record pushed for an atom wins, so a later duplicate
/// never replaces an earlier assignment.
#[derive(Debug, Default)]
pub struct AnisotropyTableBuilder {
    records: Vec<AnisotropicDisplacement>,
}

impl AnisotropyTableBuilder {
    /// Creates an empty builder.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            records: Vec::new(),
        }
    }

    /// Adds one ellipsoid. The first assignment wins on duplicates.
    pub fn push(&mut self, record: AnisotropicDisplacement) {
        self.records.push(record);
    }

    /// Packs records into immutable columns.
    #[must_use]
    pub fn finish(self) -> AnisotropyTable {
        self.finish_with_availability(true)
    }

    /// Packs records while preserving whether the category was resolved.
    ///
    /// A merge can retain explicit ellipsoids from one input while marking the
    /// combined table unavailable because another input carried no data.
    #[must_use]
    pub fn finish_with_availability(mut self, available: bool) -> AnisotropyTable {
        // The sort is stable and `dedup_by_key` keeps the first of each run, so
        // the earliest record for an atom wins exactly as an insert-if-absent
        // would, and the output stays in ascending atom order.
        self.records.sort_by_key(|record| record.atom.get());
        self.records.dedup_by_key(|record| record.atom.get());

        let mut atom = Vec::with_capacity(self.records.len());
        let mut u = Vec::with_capacity(self.records.len());
        for record in self.records {
            atom.push(record.atom);
            u.push(record.u);
        }
        AnisotropyTable {
            atom: Arc::new(atom),
            u: Arc::new(u),
            available,
        }
    }
}

#[cfg(test)]
#[path = "model_tests.rs"]
mod tests;
