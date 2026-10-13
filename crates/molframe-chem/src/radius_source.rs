//! Where the radii of a calculation come from: a built-in set or the caller's.
//!
//! Algorithms that need radii take a [`RadiusSource`] and ask it for an
//! element's radius; none of them names a table. The defaults are the built-in
//! sets of [`radii_defaults`](crate::radii_defaults), and a caller's own table
//! ([`CustomRadii`]) replaces them without touching any algorithm. Either way
//! the source can state its [`TableIdentity`], which is what a result records.

use crate::{CustomRadii, RadiusSet, TableIdentity, vdw_radius};
use molframe_core::Element;
use std::sync::Arc;

/// A source of van der Waals radii.
#[derive(Clone, Debug, PartialEq)]
pub enum RadiusSource {
    /// One of the built-in sets.
    Set(RadiusSet),
    /// A table the caller supplied, possibly layered over a built-in set.
    Custom(Arc<CustomRadii>),
}

impl Default for RadiusSource {
    fn default() -> Self {
        Self::Set(RadiusSet::Bondi)
    }
}

impl From<RadiusSet> for RadiusSource {
    fn from(set: RadiusSet) -> Self {
        Self::Set(set)
    }
}

impl From<CustomRadii> for RadiusSource {
    fn from(table: CustomRadii) -> Self {
        Self::Custom(Arc::new(table))
    }
}

impl RadiusSource {
    /// The radius of `element` in ångström, or `None` when the source has none.
    #[must_use]
    pub fn radius(&self, element: Element) -> Option<f32> {
        match self {
            Self::Set(set) => vdw_radius(element, *set),
            Self::Custom(table) => table.radius(element),
        }
    }

    /// Name, version, source and licence of the table in force.
    #[must_use]
    pub fn identity(&self) -> TableIdentity {
        match self {
            Self::Set(set) => set.identity(),
            Self::Custom(table) => table.identity().clone(),
        }
    }

    /// The built-in set behind this source, when it is one.
    #[must_use]
    pub const fn set(&self) -> Option<RadiusSet> {
        match self {
            Self::Set(set) => Some(*set),
            Self::Custom(_) => None,
        }
    }
}

#[cfg(test)]
#[path = "radius_source_tests.rs"]
mod tests;
