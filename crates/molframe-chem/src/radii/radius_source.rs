//! Where the radii of a calculation come from: the table that was selected.
//!
//! Algorithms that need radii hold a [`RadiusSource`] and ask it for an element's
//! radius; none of them names a table or a value. The source wraps the decoded
//! table the store selected (by the profile's `vdw-radii:` line, by name, or one the
//! caller built), so a caller's own table replaces any other without touching an
//! algorithm. Either way the source states its [`TableIdentity`], which is what a
//! result records.
//!
//! Cost: cloning is an `Arc` clone; a lookup is one array index.

use crate::{Provided, TableError, TableIdentity, TableRequest, TableStore, VdwRadii};
use molframe_core::Element;
use std::sync::Arc;

/// A source of van der Waals radii.
#[derive(Clone, Debug, PartialEq)]
pub struct RadiusSource {
    table: Arc<Provided<VdwRadii>>,
}

impl From<Arc<Provided<VdwRadii>>> for RadiusSource {
    fn from(table: Arc<Provided<VdwRadii>>) -> Self {
        Self { table }
    }
}

impl From<Provided<VdwRadii>> for RadiusSource {
    fn from(table: Provided<VdwRadii>) -> Self {
        Self::from(Arc::new(table))
    }
}

impl RadiusSource {
    /// Resolves `request` against `store`.
    ///
    /// # Errors
    ///
    /// [`TableError`] when nothing selects a table or the selected file is refused.
    pub fn resolve(
        request: &TableRequest<VdwRadii>,
        store: &TableStore,
    ) -> Result<Self, TableError> {
        request.resolve(store, "--radii").map(Self::from)
    }

    /// The request that selects exactly this table again.
    #[must_use]
    pub fn request(&self) -> TableRequest<VdwRadii> {
        TableRequest::Given(Arc::clone(&self.table))
    }

    /// The radius of `element` in ångström, or `None` when the table lists none.
    #[must_use]
    pub fn radius(&self, element: Element) -> Option<f32> {
        self.table.table.radius(element)
    }

    /// Name, version, source and licence of the table in force.
    #[must_use]
    pub fn identity(&self) -> TableIdentity {
        self.table.identity.clone()
    }
}

#[cfg(test)]
#[path = "radius_source_tests.rs"]
mod tests;
