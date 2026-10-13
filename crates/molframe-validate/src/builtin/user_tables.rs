//! Loading user-supplied reference tables from bytes or files.
//!
//! The parsing is the same code that decodes the bundled tables, so a user
//! table is held to the same format rules. Files are described by their path
//! unless the caller names the table. Nothing here knows about Python or the
//! command line; those layers only pass a path or bytes through.

use std::path::Path;

use super::reference_data::{RamaGrids, RotamerTables};
use super::reference_tables::{ChiDefinitions, GeometryTargets};
use super::{BuiltinError, Part, ReferenceOverrides, TableInfo};

/// Description used for a table read from `path`, unless `info` is given.
fn describe(kind: &str, path: &Path, info: Option<TableInfo>) -> TableInfo {
    match info {
        Some(info) => info,
        None => TableInfo::new(
            format!("user-{kind}"),
            "user-supplied",
            path.display().to_string(),
            "unspecified",
        ),
    }
}

fn read(path: &Path) -> Result<Vec<u8>, BuiltinError> {
    std::fs::read(path).map_err(|error| BuiltinError::Load {
        path: path.display().to_string(),
        reason: error.to_string(),
    })
}

fn text(path: &Path) -> Result<String, BuiltinError> {
    String::from_utf8(read(path)?).map_err(|error| BuiltinError::Load {
        path: path.display().to_string(),
        reason: error.to_string(),
    })
}

fn nonempty<T>(
    path: &Path,
    empty: bool,
    kind: &str,
    part: Part<T>,
) -> Result<Part<T>, BuiltinError> {
    if empty {
        Err(BuiltinError::Load {
            path: path.display().to_string(),
            reason: format!("no {kind} rows could be parsed"),
        })
    } else {
        Ok(part)
    }
}

impl ReferenceOverrides {
    /// Reads whichever tables have a path; every other table stays default.
    ///
    /// # Errors
    ///
    /// Returns [`BuiltinError`] if a file cannot be read or is malformed.
    pub fn from_paths(
        rama: Option<&Path>,
        rotamers: Option<&Path>,
        geometry: Option<&Path>,
        chi: Option<&Path>,
    ) -> Result<Self, BuiltinError> {
        let mut out = Self::default();
        if let Some(path) = rama {
            let grids = RamaGrids::from_bytes(&read(path)?)?;
            out.rama = Some(Part::new(describe("ramachandran", path, None), grids));
        }
        if let Some(path) = rotamers {
            let tables = RotamerTables::from_bytes(&read(path)?)?;
            out.rotamers = Some(Part::new(describe("rotamer", path, None), tables));
        }
        if let Some(path) = geometry {
            let targets = GeometryTargets::parse(&text(path)?);
            let part = Part::new(describe("geometry", path, None), targets);
            out.geometry = Some(nonempty(path, part.is_empty(), "geometry", part)?);
        }
        if let Some(path) = chi {
            let definitions = ChiDefinitions::parse(&text(path)?);
            let part = Part::new(describe("chi", path, None), definitions);
            out.chi = Some(nonempty(path, part.is_empty(), "chi", part)?);
        }
        Ok(out)
    }

    /// Replaces the Ramachandran grids with a decoded table and its description.
    ///
    /// # Errors
    ///
    /// Returns [`BuiltinError`] if the bytes are not a well-formed table.
    pub fn with_rama_bytes(mut self, bytes: &[u8], info: TableInfo) -> Result<Self, BuiltinError> {
        self.rama = Some(Part::new(info, RamaGrids::from_bytes(bytes)?));
        Ok(self)
    }

    /// Replaces the rotamer tables with a decoded table and its description.
    ///
    /// # Errors
    ///
    /// Returns [`BuiltinError`] if the bytes are not a well-formed table.
    pub fn with_rotamer_bytes(
        mut self,
        bytes: &[u8],
        info: TableInfo,
    ) -> Result<Self, BuiltinError> {
        self.rotamers = Some(Part::new(info, RotamerTables::from_bytes(bytes)?));
        Ok(self)
    }

    /// Replaces the bond and angle targets with a parsed table.
    #[must_use]
    pub fn with_geometry_text(mut self, text: &str, info: TableInfo) -> Self {
        self.geometry = Some(Part::new(info, GeometryTargets::parse(text)));
        self
    }

    /// Replaces the chi definitions with a parsed table.
    #[must_use]
    pub fn with_chi_text(mut self, text: &str, info: TableInfo) -> Self {
        self.chi = Some(Part::new(info, ChiDefinitions::parse(text)));
        self
    }
}

#[cfg(test)]
#[path = "user_tables_tests.rs"]
mod tests;
