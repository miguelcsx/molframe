//! Typed positions into the structure's tables.
//!
//! Every level of the hierarchy gets its own index type. They share a
//! representation but never a type, so handing a residue position to something
//! expecting an atom position is a compile error rather than a wrong number
//! discovered downstream.
//!
//! The representation is `u32`. Four billion atoms is beyond any structure that
//! will ever be deposited, and the narrower width halves the memory of every
//! column that stores a position — of which there are several, one per atom.

use std::fmt;

macro_rules! index_newtype {
    ($(#[$meta:meta])* $name:ident, $label:literal) => {
        $(#[$meta])*
        #[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
        #[repr(transparent)]
        pub struct $name(u32);

        impl $name {
            #[doc = concat!("Creates a ", $label, " position from its ordinal.")]
            ///
            /// # Examples
            ///
            /// ```
            #[doc = concat!("use pdbiox_core::", stringify!($name), ";")]
            #[doc = concat!("assert_eq!(", stringify!($name), "::new(7).get(), 7);")]
            /// ```
            #[must_use]
            pub const fn new(ordinal: u32) -> Self {
                Self(ordinal)
            }

            #[doc = concat!("Returns the ordinal of this ", $label, " position.")]
            #[must_use]
            pub const fn get(self) -> u32 {
                self.0
            }

            #[doc = concat!("Returns the ordinal widened for slice indexing.")]
            #[must_use]
            pub const fn as_usize(self) -> usize {
                self.0 as usize
            }

            #[doc = concat!("Returns the next ", $label, " position, or `None` at the representation limit.")]
            #[must_use]
            pub const fn next(self) -> Option<Self> {
                match self.0.checked_add(1) {
                    Some(next) => Some(Self(next)),
                    None => None,
                }
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}({})", $label, self.0)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{}", self.0)
            }
        }
    };
}

index_newtype!(
    /// Position of an atom row in the structure's flat atom order.
    ///
    /// Atom order is the order the source file presented, never a reordering
    /// pdbiox chose, so this position is stable for the lifetime of a snapshot.
    AtomIndex,
    "atom"
);
index_newtype!(
    /// Position of a residue in the residue table.
    ResidueIndex,
    "residue"
);
index_newtype!(
    /// Position of a chain in the chain table.
    ChainIndex,
    "chain"
);
index_newtype!(
    /// Position of an entity in the entity table.
    ///
    /// An entity is a distinct chemical species, and several chains may be
    /// copies of one. It is therefore not a level of the hierarchy: chains point
    /// at entities, and nothing is contained by an entity.
    EntityIndex,
    "entity"
);
index_newtype!(
    /// Position of a model in the model table.
    ///
    /// This is the model's position, not the number it was deposited under.
    /// The deposited number is stored separately and is never renumbered.
    ModelIndex,
    "model"
);
index_newtype!(
    /// Position of a bond in the bond table.
    BondIndex,
    "bond"
);
index_newtype!(
    /// Identifier of one generated copy of an assembly's contents.
    ///
    /// Assemblies name their own instances; materialising one keeps that name
    /// rather than assigning a fresh sequence.
    InstanceId,
    "instance"
);
#[cfg(test)]
#[path = "typed_tests.rs"]
mod tests;
