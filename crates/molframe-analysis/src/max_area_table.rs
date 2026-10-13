//! A table of maximum residue areas, built-in or the caller's own.
//!
//! Relative accessibility is an area divided by a maximum, and which maximum is
//! a choice. [`MaxAreaTable`] holds one set of maxima together with the
//! [`TableIdentity`] that says where it came from; the built-in Tien et al.
//! table is [`MaxAreaTable::builtin`], and [`MaxAreaTable::parse`] reads a
//! caller's own from text:
//!
//! ```text
//! name: lab-max-area
//! version: 1
//! source: Gly-X-Gly ensemble, our force field
//! licence: CC0-1.0
//!
//! # residue maximum-area-in-square-angstrom
//! ALA 130.5
//! GLY 105
//! ```
//!
//! `name`, `version`, `source` and `licence` are required and every row is a
//! residue name and a positive finite area. A residue the table lacks has no
//! relative accessibility. Protonation variants (`HID`, `CYX`, ...) take their
//! parent's value. Every problem is reported with its line number and nothing
//! is accepted partially.

use crate::max_area_defaults::{TIEN_2013, TIEN_2013_AREAS};
use crate::protonation::parent_residue;
use crate::table_text::{IdentityFields, TableTextError};
use molframe_chem::TableIdentity;
use std::collections::BTreeMap;
use std::sync::OnceLock;

/// The largest area a table may declare, in square ångström; anything above is
/// a unit mistake.
const MAXIMUM_AREA: f64 = 2000.0;

/// Why a maximum-area table was refused.
pub type MaxAreaTableError = TableTextError;

const WHAT: &str = "maximum-area table";

/// Maximum accessible area per residue name, with the identity of its source.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxAreaTable {
    identity: TableIdentity,
    areas: BTreeMap<String, f64>,
}

impl Default for MaxAreaTable {
    /// A copy of the Tien et al. (2013) theoretical maxima; prefer
    /// [`MaxAreaTable::builtin`], which borrows them.
    fn default() -> Self {
        Self::builtin().clone()
    }
}

impl MaxAreaTable {
    /// The Tien et al. (2013) theoretical maxima, built once and borrowed for the
    /// life of the process.
    #[must_use]
    pub fn builtin() -> &'static Self {
        static BUILTIN: OnceLock<MaxAreaTable> = OnceLock::new();
        BUILTIN.get_or_init(|| Self {
            identity: TIEN_2013,
            areas: TIEN_2013_AREAS
                .iter()
                .map(|&(name, area)| (name.to_owned(), area))
                .collect(),
        })
    }

    /// Parses and validates a table.
    ///
    /// # Errors
    ///
    /// [`MaxAreaTableError`] for a missing, unknown or repeated header, a
    /// repeated residue, or an area that is not a positive finite number.
    pub fn parse(text: &str) -> Result<Self, MaxAreaTableError> {
        let mut fields = IdentityFields::default();
        let mut areas = BTreeMap::new();
        for (index, raw) in text.lines().enumerate() {
            let number = index + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let Some((key, value)) = line.split_once(':') {
                if !fields.set(WHAT, number, key.trim(), value.trim())? {
                    return Err(MaxAreaTableError::at(
                        WHAT,
                        number,
                        format!(
                            "unknown header `{}`; expected name, version, source or licence",
                            key.trim()
                        ),
                    ));
                }
            } else {
                read_row(number, line, &mut areas)?;
            }
        }
        if areas.is_empty() {
            return Err(MaxAreaTableError::at(WHAT, 0, "the table lists no residue"));
        }
        Ok(Self {
            identity: fields.finish(WHAT)?,
            areas,
        })
    }

    /// The declared identity of this table.
    #[must_use]
    pub const fn identity(&self) -> &TableIdentity {
        &self.identity
    }

    /// The maximum area of residue `name`, Å², or `None` when the table has
    /// no value for it.
    #[must_use]
    pub fn max_area(&self, name: &str) -> Option<f64> {
        self.areas
            .get(name)
            .or_else(|| self.areas.get(parent_residue(name)))
            .copied()
    }

    /// `area` as a fraction of the residue's maximum, or `None` without one.
    #[must_use]
    pub fn relative(&self, name: &str, area: f64) -> Option<f64> {
        self.max_area(name).map(|maximum| area / maximum)
    }
}

fn read_row(
    line: usize,
    text: &str,
    areas: &mut BTreeMap<String, f64>,
) -> Result<(), MaxAreaTableError> {
    let mut words = text.split_whitespace();
    let (Some(name), Some(value), None) = (words.next(), words.next(), words.next()) else {
        return Err(MaxAreaTableError::at(
            WHAT,
            line,
            "expected `residue area`, two words",
        ));
    };
    let area: f64 = value
        .parse()
        .map_err(|_| MaxAreaTableError::at(WHAT, line, format!("`{value}` is not a number")))?;
    if !area.is_finite() || area <= 0.0 || area > MAXIMUM_AREA {
        return Err(MaxAreaTableError::at(
            WHAT,
            line,
            format!("area {value} must be positive and at most {MAXIMUM_AREA}"),
        ));
    }
    if areas.insert(name.to_ascii_uppercase(), area).is_some() {
        return Err(MaxAreaTableError::at(
            WHAT,
            line,
            format!("residue `{name}` is listed twice"),
        ));
    }
    Ok(())
}

#[cfg(test)]
#[path = "max_area_table_tests.rs"]
mod tests;
