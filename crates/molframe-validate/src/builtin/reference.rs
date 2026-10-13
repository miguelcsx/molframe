//! The reference tables a built-in validation runs against.
//!
//! The scoring code is generic: it asks a [`ReferenceSet`] for each table it
//! needs and reads every parameter from it. The set holds one
//! [`TableRequest`] per kind (the profile's choice, a path or cache entry name, or
//! a table already decoded) and a [`TableStore`] that resolves them, so the
//! library compiles in no table: a check whose table is not selected fails with
//! [`TableError::Missing`](molframe_chem::TableError::Missing) naming the option,
//! cache entry or fetch command that supplies it. Every resolved table carries
//! its [`TableIdentity`](molframe_chem::TableIdentity), which results record.
//!
//! Cost: a request resolves in an `Arc` clone once its table has been decoded
//! (see [`TableStore`]); nothing is copied per call.

use std::sync::Arc;

use molframe_chem::{Provided, TableRequest, TableStore};

use super::BuiltinError;
use super::nucleic_policy::NucleicPolicy;
use super::reference_data::{RamaGrids, RotamerTables};
use super::reference_tables::{ChiDefinitions, GeometryTargets};

/// A decoded table with its identity, shared.
pub type Part<T> = Arc<Provided<T>>;

/// Which table each built-in check reads, and the store that resolves them.
#[derive(Clone, Debug)]
pub struct ReferenceSet {
    store: Arc<TableStore>,
    /// Ramachandran density grids and class cut-offs.
    pub rama: TableRequest<RamaGrids>,
    /// Side-chain rotamer classification tables.
    pub rotamers: TableRequest<RotamerTables>,
    /// Bond-length and bond-angle targets.
    pub geometry: TableRequest<GeometryTargets>,
    /// Chi-angle atom definitions and the rotamer table each uses.
    pub chi: TableRequest<ChiDefinitions>,
    /// Nucleic pucker and suite windows.
    pub nucleic: TableRequest<NucleicPolicy>,
}

impl Default for ReferenceSet {
    fn default() -> Self {
        Self::with_store(TableStore::global())
    }
}

impl ReferenceSet {
    /// A set whose every request is the profile's, resolved through `store`.
    #[must_use]
    pub fn with_store(store: Arc<TableStore>) -> Self {
        Self {
            store,
            rama: TableRequest::profile(),
            rotamers: TableRequest::profile(),
            geometry: TableRequest::profile(),
            chi: TableRequest::profile(),
            nucleic: TableRequest::profile(),
        }
    }

    /// The Ramachandran grids.
    ///
    /// # Errors
    ///
    /// [`BuiltinError::Table`] when no table is selected or the selected one is bad.
    pub fn rama(&self) -> Result<Part<RamaGrids>, BuiltinError> {
        Ok(self.rama.resolve(&self.store, "--rama-table")?)
    }

    /// The rotamer tables.
    ///
    /// # Errors
    ///
    /// As [`ReferenceSet::rama`].
    pub fn rotamers(&self) -> Result<Part<RotamerTables>, BuiltinError> {
        Ok(self.rotamers.resolve(&self.store, "--rotamer-table")?)
    }

    /// The bond and angle targets.
    ///
    /// # Errors
    ///
    /// As [`ReferenceSet::rama`].
    pub fn geometry(&self) -> Result<Part<GeometryTargets>, BuiltinError> {
        Ok(self.geometry.resolve(&self.store, "--geometry-table")?)
    }

    /// The chi definitions.
    ///
    /// # Errors
    ///
    /// As [`ReferenceSet::rama`].
    pub fn chi(&self) -> Result<Part<ChiDefinitions>, BuiltinError> {
        Ok(self.chi.resolve(&self.store, "--chi-table")?)
    }

    /// The nucleic policy windows.
    ///
    /// # Errors
    ///
    /// As [`ReferenceSet::rama`].
    pub fn nucleic(&self) -> Result<Part<NucleicPolicy>, BuiltinError> {
        Ok(self.nucleic.resolve(&self.store, "--nucleic-table")?)
    }
}

#[cfg(test)]
#[path = "reference_tests.rs"]
mod tests;
