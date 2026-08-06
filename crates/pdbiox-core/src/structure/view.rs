//! A subset of a structure, costing two words and a selection.
//!
//! A view copies nothing. It holds a reference to the structure, the set of
//! atoms it covers, and the coordinate generation it was made against. Composing
//! views is set algebra on the selections, so narrowing a view never touches the
//! atoms it discards.
//!
//! The generation is what makes a stale view detectable. A view made before an
//! edit still refers to the snapshot it was made from, and asking it about a
//! newer structure is refused rather than answered with positions that have
//! since moved.

use super::data::{Structure, StructureData};
use crate::coords::{Aabb, CoordinateGeneration};
use crate::index::ModelIndex;
use crate::selection::AtomSelection;
use std::sync::Arc;

/// A zero-copy subset of a structure.
///
/// # Examples
///
/// ```
/// use pdbiox_core::structure::{Structure, StructureData};
/// use pdbiox_core::AtomSelection;
///
/// let structure = Structure::new(StructureData::empty());
/// let everything = structure.view();
/// let nothing = everything.narrow(&AtomSelection::Empty);
///
/// assert!(nothing.is_empty());
/// assert_eq!(everything.len(), structure.atom_count());
/// ```
#[derive(Clone, Debug)]
pub struct StructureView {
    structure: Arc<StructureData>,
    atoms: AtomSelection,
    generation: CoordinateGeneration,
}

impl StructureView {
    /// Creates a view over `atoms` of `structure`.
    #[must_use]
    pub fn new(structure: &Structure, atoms: AtomSelection) -> Self {
        Self {
            structure: structure.shared(),
            atoms,
            generation: structure.generation(),
        }
    }

    /// The atoms this view covers.
    #[must_use]
    pub const fn selection(&self) -> &AtomSelection {
        &self.atoms
    }

    /// The structure this view is over.
    #[must_use]
    pub fn data(&self) -> &StructureData {
        &self.structure
    }

    /// The number of atoms covered.
    #[must_use]
    pub fn len(&self) -> u32 {
        self.atoms.len()
    }

    /// Returns true when the view covers no atoms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// The generation this view was made against.
    #[must_use]
    pub const fn generation(&self) -> CoordinateGeneration {
        self.generation
    }

    /// Returns true when `structure` has moved on since this view was made.
    ///
    /// A view of an earlier snapshot is still perfectly usable *against that
    /// snapshot* — it holds its own reference to it. This asks a different
    /// question: whether the view may be used to talk about the newer one.
    #[must_use]
    pub fn is_stale_for(&self, structure: &Structure) -> bool {
        self.generation != structure.generation()
    }

    /// A view of the atoms in both this view and `other`.
    #[must_use]
    pub fn narrow(&self, other: &AtomSelection) -> Self {
        self.with_selection(self.atoms.intersect(other))
    }

    /// A view of the atoms in either this view or `other`.
    #[must_use]
    pub fn union(&self, other: &Self) -> Self {
        self.with_selection(self.atoms.union(&other.atoms))
    }

    /// A view of the atoms in both views.
    #[must_use]
    pub fn intersect(&self, other: &Self) -> Self {
        self.with_selection(self.atoms.intersect(&other.atoms))
    }

    /// A view of the atoms in this view but not the other.
    #[must_use]
    pub fn difference(&self, other: &Self) -> Self {
        self.with_selection(self.atoms.difference(&other.atoms))
    }

    fn with_selection(&self, atoms: AtomSelection) -> Self {
        Self {
            structure: Arc::clone(&self.structure),
            atoms,
            generation: self.generation,
        }
    }

    /// The positions of the covered atoms in one model.
    ///
    /// Yields values rather than a slice: the atoms of a view are not in general
    /// contiguous, so a slice would either be a lie or a copy.
    pub fn positions(&self, model: ModelIndex) -> impl Iterator<Item = [f32; 3]> + '_ {
        let block = self.structure.coords.block(model);
        self.atoms
            .iter()
            .filter_map(move |atom| block?.as_slice().get(atom as usize).copied())
    }

    /// The box enclosing the covered atoms in one model.
    #[must_use]
    pub fn bounds(&self, model: ModelIndex) -> Aabb {
        let mut bounds = Aabb::EMPTY;
        for position in self.positions(model) {
            bounds.extend(position);
        }
        bounds
    }
}

impl Structure {
    /// A view of every atom.
    #[must_use]
    pub fn view(&self) -> StructureView {
        StructureView::new(self, AtomSelection::All(self.atom_count()))
    }

    /// A view of the atoms in `selection`.
    #[must_use]
    pub fn view_of(&self, selection: AtomSelection) -> StructureView {
        StructureView::new(self, selection)
    }
}

#[cfg(test)]
#[path = "view_tests.rs"]
mod tests;
