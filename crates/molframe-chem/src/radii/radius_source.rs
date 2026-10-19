//! Where the radii of a calculation come from: the table that was selected.
//!
//! Algorithms that need radii hold a [`RadiusSource`] and ask it for an element's
//! radius; none of them names a table or a value. The source wraps the decoded
//! table the store selected (by the profile's `vdw-radii:` line, by name, or one the
//! caller built), so a caller's own table replaces any other without touching an
//! algorithm. Either way the source states its [`TableIdentity`], which is what a
//! result records.
//!
//! The radius role is a type parameter. [`RadiusSource`] is the van der Waals source
//! the analyses measure atoms with; [`ElementRadiusSource<Covalent>`] selects the
//! covalent table the same way. Both resolve and look up identically, so the role
//! changes the table kind and nothing else.
//!
//! Cost: cloning is an `Arc` clone; a lookup is one array index.

use crate::{
    ElementRadii, Provided, RadiusRole, TableError, TableIdentity, TableRequest, TableStore, Vdw,
};
use molframe_core::Element;
use std::sync::Arc;

/// A source of the radii of one role, held as the one table that was selected.
#[derive(Clone, Debug, PartialEq)]
pub struct ElementRadiusSource<R> {
    table: Arc<Provided<ElementRadii<R>>>,
}

/// A source of van der Waals radii, the radii the analyses measure atoms with.
pub type RadiusSource = ElementRadiusSource<Vdw>;

impl<R> From<Arc<Provided<ElementRadii<R>>>> for ElementRadiusSource<R> {
    fn from(table: Arc<Provided<ElementRadii<R>>>) -> Self {
        Self { table }
    }
}

impl<R> From<Provided<ElementRadii<R>>> for ElementRadiusSource<R> {
    fn from(table: Provided<ElementRadii<R>>) -> Self {
        Self::from(Arc::new(table))
    }
}

impl<R: RadiusRole + 'static> ElementRadiusSource<R> {
    /// Resolves `request` against `store`.
    ///
    /// # Errors
    ///
    /// [`TableError`] when nothing selects a table or the selected file is refused.
    pub fn resolve(
        request: &TableRequest<ElementRadii<R>>,
        store: &TableStore,
    ) -> Result<Self, TableError> {
        request.resolve(store, "--radii").map(Self::from)
    }

    /// The request that selects exactly this table again.
    #[must_use]
    pub fn request(&self) -> TableRequest<ElementRadii<R>> {
        TableRequest::Given(Arc::clone(&self.table))
    }
}

impl<R> ElementRadiusSource<R> {
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
