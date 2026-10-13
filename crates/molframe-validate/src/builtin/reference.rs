//! The set of reference tables a built-in validation runs against.
//!
//! The scoring code is generic: it receives a [`ReferenceSet`] and reads
//! every parameter from it. Each table travels with a [`TableInfo`] naming its
//! source, version and licence, so a result can report exactly what was
//! applied. The bundled defaults live in the `bundled` module; a caller
//! replaces any subset with its own tables through [`ReferenceOverrides`],
//! and the remaining tables stay the defaults.

use std::ops::Deref;
use std::sync::Arc;

use super::reference_data::{RamaGrids, RotamerTables};
use super::reference_tables::{ChiDefinitions, GeometryTargets};

/// Name, version, source and licence of one reference table.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TableInfo {
    /// Short stable name of the table.
    pub name: String,
    /// Version label, for example a commit or release.
    pub version: String,
    /// Where the data comes from.
    pub source: String,
    /// Licence the data is distributed under.
    pub licence: String,
}

impl TableInfo {
    /// Describes a table.
    #[must_use]
    pub fn new(
        name: impl Into<String>,
        version: impl Into<String>,
        source: impl Into<String>,
        licence: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            version: version.into(),
            source: source.into(),
            licence: licence.into(),
        }
    }
}

impl std::fmt::Display for TableInfo {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            formatter,
            "{} {} ({}; {})",
            self.name, self.version, self.source, self.licence
        )
    }
}

/// A table together with the description of where it came from.
#[derive(Debug)]
pub struct Part<T> {
    info: TableInfo,
    value: Arc<T>,
}

impl<T> Clone for Part<T> {
    fn clone(&self) -> Self {
        Self {
            info: self.info.clone(),
            value: Arc::clone(&self.value),
        }
    }
}

impl<T> Part<T> {
    /// Pairs a table with its description.
    #[must_use]
    pub fn new(info: TableInfo, value: T) -> Self {
        Self {
            info,
            value: Arc::new(value),
        }
    }

    /// Pairs an already shared table with its description.
    #[must_use]
    pub const fn shared(info: TableInfo, value: Arc<T>) -> Self {
        Self { info, value }
    }

    /// The description of this table.
    #[must_use]
    pub const fn info(&self) -> &TableInfo {
        &self.info
    }
}

impl<T> Deref for Part<T> {
    type Target = T;

    fn deref(&self) -> &T {
        &self.value
    }
}

/// Every table the built-in protein checks read.
#[derive(Clone, Debug)]
pub struct ReferenceSet {
    /// Ramachandran density grids and class cut-offs.
    pub rama: Part<RamaGrids>,
    /// Side-chain rotamer classification tables.
    pub rotamers: Part<RotamerTables>,
    /// Bond-length and bond-angle targets.
    pub geometry: Part<GeometryTargets>,
    /// Chi-angle atom definitions and the rotamer table each uses.
    pub chi: Part<ChiDefinitions>,
}

impl ReferenceSet {
    /// The description of every table, in a fixed order.
    #[must_use]
    pub fn infos(&self) -> Vec<TableInfo> {
        vec![
            self.rama.info().clone(),
            self.rotamers.info().clone(),
            self.geometry.info().clone(),
            self.chi.info().clone(),
        ]
    }

    /// Replaces the tables present in `overrides`; the rest stay as they are.
    #[must_use]
    pub fn overlaid(self, overrides: ReferenceOverrides) -> Self {
        Self {
            rama: pick(overrides.rama, self.rama),
            rotamers: pick(overrides.rotamers, self.rotamers),
            geometry: pick(overrides.geometry, self.geometry),
            chi: pick(overrides.chi, self.chi),
        }
    }
}

fn pick<T>(replacement: Option<Part<T>>, default: Part<T>) -> Part<T> {
    match replacement {
        Some(part) => part,
        None => default,
    }
}

/// Tables a caller supplies in place of the defaults.
#[derive(Clone, Debug, Default)]
pub struct ReferenceOverrides {
    /// Replacement Ramachandran grids.
    pub rama: Option<Part<RamaGrids>>,
    /// Replacement rotamer tables.
    pub rotamers: Option<Part<RotamerTables>>,
    /// Replacement bond and angle targets.
    pub geometry: Option<Part<GeometryTargets>>,
    /// Replacement chi definitions.
    pub chi: Option<Part<ChiDefinitions>>,
}

#[cfg(test)]
#[path = "reference_tests.rs"]
mod tests;
