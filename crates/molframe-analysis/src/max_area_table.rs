//! A table of maximum residue areas, built-in or the caller's own.
//!
//! Relative accessibility is an area divided by a maximum, and which maximum is
//! a choice. [`MaxAreaTable`] holds one set of maxima together with the
//! [`TableIdentity`] that says where it came from; the built-in Tien et al.
//! table is [`MaxAreaTable::default`], and [`MaxAreaTable::parse`] reads a
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
use core::fmt;
use molframe_chem::TableIdentity;
use std::collections::BTreeMap;

/// The largest area a table may declare, in square ångström; anything above is
/// a unit mistake.
const MAXIMUM_AREA: f64 = 2000.0;

/// Why a maximum-area table was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MaxAreaTableError {
    line: usize,
    message: String,
}

impl MaxAreaTableError {
    fn at(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }

    /// The one-based line the problem is on, or zero for a whole-table problem.
    #[must_use]
    pub const fn line(&self) -> usize {
        self.line
    }
}

impl fmt::Display for MaxAreaTableError {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(out, "maximum-area table: {}", self.message)
        } else {
            write!(
                out,
                "maximum-area table, line {}: {}",
                self.line, self.message
            )
        }
    }
}

impl core::error::Error for MaxAreaTableError {}

/// Maximum accessible area per residue name, with the identity of its source.
#[derive(Clone, Debug, PartialEq)]
pub struct MaxAreaTable {
    identity: TableIdentity,
    areas: BTreeMap<String, f64>,
}

impl Default for MaxAreaTable {
    /// The Tien et al. (2013) theoretical maxima.
    fn default() -> Self {
        Self {
            identity: TIEN_2013,
            areas: TIEN_2013_AREAS
                .iter()
                .map(|&(name, area)| (name.to_owned(), area))
                .collect(),
        }
    }
}

impl MaxAreaTable {
    /// Parses and validates a table.
    ///
    /// # Errors
    ///
    /// [`MaxAreaTableError`] for a missing, unknown or repeated header, a
    /// repeated residue, or an area that is not a positive finite number.
    pub fn parse(text: &str) -> Result<Self, MaxAreaTableError> {
        let mut fields: [Option<String>; 4] = Default::default();
        let mut areas = BTreeMap::new();
        for (index, raw) in text.lines().enumerate() {
            let number = index + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let Some((key, value)) = line.split_once(':') {
                set_header(number, key.trim(), value.trim(), &mut fields)?;
            } else {
                read_row(number, line, &mut areas)?;
            }
        }
        if areas.is_empty() {
            return Err(MaxAreaTableError::at(0, "the table lists no residue"));
        }
        let [name, version, source, licence] = fields;
        let need = |value: Option<String>, key: &str| {
            value.ok_or_else(|| MaxAreaTableError::at(0, format!("missing `{key}` header")))
        };
        Ok(Self {
            identity: TableIdentity::supplied(
                &need(name, "name")?,
                &need(version, "version")?,
                &need(source, "source")?,
                &need(licence, "licence")?,
            ),
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

fn set_header(
    line: usize,
    key: &str,
    value: &str,
    fields: &mut [Option<String>; 4],
) -> Result<(), MaxAreaTableError> {
    if value.is_empty() {
        return Err(MaxAreaTableError::at(line, format!("`{key}` has no value")));
    }
    let slot = match key {
        "name" => 0,
        "version" => 1,
        "source" => 2,
        "licence" | "license" => 3,
        other => {
            return Err(MaxAreaTableError::at(
                line,
                format!("unknown header `{other}`; expected name, version, source or licence"),
            ));
        }
    };
    match fields.get_mut(slot) {
        Some(field) if field.is_none() => *field = Some(value.to_owned()),
        _ => {
            return Err(MaxAreaTableError::at(
                line,
                format!("`{key}` is declared twice"),
            ));
        }
    }
    Ok(())
}

fn read_row(
    line: usize,
    text: &str,
    areas: &mut BTreeMap<String, f64>,
) -> Result<(), MaxAreaTableError> {
    let mut words = text.split_whitespace();
    let (Some(name), Some(value), None) = (words.next(), words.next(), words.next()) else {
        return Err(MaxAreaTableError::at(
            line,
            "expected `residue area`, two words",
        ));
    };
    let area: f64 = value
        .parse()
        .map_err(|_| MaxAreaTableError::at(line, format!("`{value}` is not a number")))?;
    if !area.is_finite() || area <= 0.0 || area > MAXIMUM_AREA {
        return Err(MaxAreaTableError::at(
            line,
            format!("area {value} must be positive and at most {MAXIMUM_AREA}"),
        ));
    }
    if areas.insert(name.to_ascii_uppercase(), area).is_some() {
        return Err(MaxAreaTableError::at(
            line,
            format!("residue `{name}` is listed twice"),
        ));
    }
    Ok(())
}

#[cfg(test)]
#[path = "max_area_table_tests.rs"]
mod tests;
