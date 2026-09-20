//! The curated name for a zero-copy subset of a structure.

use std::ops::{BitAnd, BitOr, Sub};

use molframe_core::index::{AtomIndex, ModelIndex};
use molframe_core::selection::AtomSelection;
use molframe_core::structure::{AtomRef, StructureView};

use super::handle::Structure;

/// A zero-copy subset of a [`Structure`]'s atoms.
///
/// This is the public name for [`molframe_core::structure::StructureView`]:
/// it costs a shared structure reference and a selection, never a copy.
/// Materializing coordinates for a (possibly discontiguous) selection is
/// [`Selection::to_coordinates`], never a plain field, so the name says when a
/// gather happens.
///
/// # Examples
///
/// ```
/// # #[cfg(feature = "query")]
/// # {
/// use molframe::prelude::*;
///
/// # const PDB: &str = "\
/// # ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N
/// # ATOM      2  CA  ALA A   1      12.560   6.195  -6.504  1.00  0.00           C
/// # END
/// # ";
/// let (structure, _) = read_bytes(PDB.into(), Some("1abc.pdb"), &ReadOptions::new())?;
/// let alphas = structure.select("name CA", &AnalysisPolicy::default())?;
/// assert_eq!(alphas.len(), 1);
/// # }
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, Debug)]
pub struct Selection(StructureView);

impl Selection {
    /// The number of atoms covered.
    #[must_use]
    pub fn len(&self) -> u64 {
        self.0.len()
    }

    /// Returns true when the selection covers no atoms.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns true when `structure` has moved on since this selection was
    /// made, so its atom positions may no longer describe `structure`.
    #[must_use]
    pub fn is_stale_for(&self, structure: &Structure) -> bool {
        self.0.is_stale_for(structure.engine())
    }

    /// The covered atoms, in structure order.
    pub fn atoms(&self) -> impl Iterator<Item = AtomRef<'_>> + '_ {
        let data = self.0.data();
        self.0
            .selection()
            .iter()
            .filter_map(move |atom| data.atom(AtomIndex::new(atom)))
    }

    /// Gathers the positions of the covered atoms in one model.
    ///
    /// This is the one place a selection copies: its atoms are not in general
    /// contiguous, so a borrowed slice would either be a lie or a copy.
    #[must_use]
    pub fn to_coordinates(&self, model: ModelIndex) -> Vec<[f32; 3]> {
        self.0.positions(model).collect()
    }

    /// The selection itself, for the kernels that take one directly.
    pub(crate) fn atom_selection(&self) -> &AtomSelection {
        self.0.selection()
    }
}

impl From<StructureView> for Selection {
    fn from(view: StructureView) -> Self {
        Self(view)
    }
}

impl BitOr for &Selection {
    type Output = Selection;

    fn bitor(self, rhs: Self) -> Selection {
        Selection(self.0.union(&rhs.0))
    }
}

impl BitAnd for &Selection {
    type Output = Selection;

    fn bitand(self, rhs: Self) -> Selection {
        Selection(self.0.intersect(&rhs.0))
    }
}

impl Sub for &Selection {
    type Output = Selection;

    fn sub(self, rhs: Self) -> Selection {
        Selection(self.0.difference(&rhs.0))
    }
}
