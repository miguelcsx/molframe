//! Thin, indexable views over a structure's chains, residues, atoms and models.
//!
//! Each collection is two words — a borrow of the snapshot and nothing else —
//! and each one hands out `Copy` handles rather than objects. Iterating a
//! million atoms allocates nothing, and `get` is a bounds check over a table
//! rather than a walk, because the tables are already position-indexed.
//!
//! The iterators are concrete types rather than `Box<dyn Iterator<…>>`: a
//! collection that is entered in a loop should not allocate on the way in, and
//! a named type is what lets `len`, `Clone` and the size hint stay exact.

use std::iter::FusedIterator;

use molframe_core::index::{AtomIndex, ChainIndex, ModelIndex, ResidueIndex};
use molframe_core::structure::{
    AtomRef, ChainRef, ModelRef, ResidueRef, Structure as CoreStructure,
};

/// Defines the counter iterator one collection hands out.
///
/// Every hierarchy level walks the same shape — a position, an end, and one
/// table lookup per step — so the four differ only in the handle they build.
macro_rules! hierarchy_iter {
    ($(#[$meta:meta])* $name:ident, $item:ty, $index:ty, $handle:ident) => {
        $(#[$meta])*
        #[derive(Clone, Debug)]
        pub struct $name<'a> {
            structure: &'a CoreStructure,
            position: usize,
            end: usize,
        }

        impl<'a> Iterator for $name<'a> {
            type Item = $item;

            fn next(&mut self) -> Option<Self::Item> {
                if self.position >= self.end {
                    return None;
                }
                let ordinal = u32::try_from(self.position).ok()?;
                self.position += 1;
                self.structure.$handle(<$index>::new(ordinal))
            }

            fn size_hint(&self) -> (usize, Option<usize>) {
                let remaining = self.end - self.position;
                (remaining, Some(remaining))
            }
        }

        impl ExactSizeIterator for $name<'_> {
            fn len(&self) -> usize {
                self.end - self.position
            }
        }

        impl FusedIterator for $name<'_> {}
    };
}

/// Defines one hierarchy collection: the range view and its iterator.
macro_rules! hierarchy_collection {
    (
        $(#[$meta:meta])* $name:ident, $iter:ident, $item:ty, $index:ty, $handle:ident
    ) => {
        $(#[$meta])*
        #[derive(Clone, Copy, Debug)]
        pub struct $name<'a> {
            structure: &'a CoreStructure,
        }

        impl<'a> $name<'a> {
            pub(super) fn new(structure: &'a CoreStructure) -> Self {
                Self { structure }
            }

            /// The number of entries.
            #[must_use]
            pub fn len(&self) -> usize {
                self.iteration_end()
            }

            /// Returns true when there are none.
            #[must_use]
            pub fn is_empty(&self) -> bool {
                self.iteration_end() == 0
            }

            /// The entry at `position`, or `None` past the end.
            ///
            /// A bounds check over the table, so this costs the same at the
            /// first position and the last.
            #[must_use]
            pub fn get(&self, position: usize) -> Option<$item> {
                let ordinal = u32::try_from(position).ok()?;
                self.structure.$handle(<$index>::new(ordinal))
            }

            /// Every entry, in structure order.
            #[must_use]
            pub fn iter(&self) -> $iter<'a> {
                self.structure_iter()
            }

            fn structure_iter(&self) -> $iter<'a> {
                $iter {
                    structure: self.structure,
                    position: 0,
                    end: self.iteration_end(),
                }
            }
        }

        impl<'a> IntoIterator for $name<'a> {
            type Item = $item;
            type IntoIter = $iter<'a>;

            fn into_iter(self) -> Self::IntoIter {
                self.structure_iter()
            }
        }

        impl<'a> IntoIterator for &$name<'a> {
            type Item = $item;
            type IntoIter = $iter<'a>;

            fn into_iter(self) -> Self::IntoIter {
                self.structure_iter()
            }
        }
    };
}

hierarchy_iter!(
    /// Every chain of a structure, in structure order.
    ChainsIter,
    ChainRef<'a>,
    ChainIndex,
    chain
);

hierarchy_collection!(
    /// Every chain of a structure, in structure order.
    Chains,
    ChainsIter,
    ChainRef<'a>,
    ChainIndex,
    chain
);

impl Chains<'_> {
    fn iteration_end(self) -> usize {
        self.structure.chain_count()
    }
}

hierarchy_iter!(
    /// Every residue of a structure, in structure order.
    ResiduesIter,
    ResidueRef<'a>,
    ResidueIndex,
    residue
);

hierarchy_collection!(
    /// Every residue of a structure, in structure order.
    Residues,
    ResiduesIter,
    ResidueRef<'a>,
    ResidueIndex,
    residue
);

impl Residues<'_> {
    fn iteration_end(self) -> usize {
        self.structure.residue_count()
    }
}

hierarchy_iter!(
    /// Every atom of a structure, in structure order.
    AtomsIter,
    AtomRef<'a>,
    AtomIndex,
    atom
);

hierarchy_collection!(
    /// Every atom of a structure, in structure order.
    Atoms,
    AtomsIter,
    AtomRef<'a>,
    AtomIndex,
    atom
);

impl Atoms<'_> {
    fn iteration_end(self) -> usize {
        self.structure.atom_count() as usize
    }
}

hierarchy_iter!(
    /// Every model of a structure, in structure order.
    ModelsIter,
    ModelRef<'a>,
    ModelIndex,
    model
);

hierarchy_collection!(
    /// Every model of a structure.
    Models,
    ModelsIter,
    ModelRef<'a>,
    ModelIndex,
    model
);

impl Models<'_> {
    fn iteration_end(self) -> usize {
        self.structure.model_count()
    }
}
