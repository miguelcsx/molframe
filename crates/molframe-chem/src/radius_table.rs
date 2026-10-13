//! A caller's own radius table, parsed and validated.
//!
//! The table is plain text so it can be written by hand, kept next to the data
//! it was used for and passed unchanged from Rust, the command line or Python.
//! It declares its own identity, so a result computed with it says which table
//! produced the numbers:
//!
//! ```text
//! name: lab-radii
//! version: 2
//! source: Smith lab, internal survey of the CSD
//! licence: CC0-1.0
//! base: bondi
//!
//! # element radius-in-angstrom
//! C 1.75
//! N 1.60
//! ```
//!
//! `name`, `version`, `source` and `licence` are required. `base` is optional
//! and layers the table over one of the built-in sets: an element the caller
//! lists takes the caller's value, any other element falls through to the base,
//! and without a base it has no radius; a table must list an element or name a
//! base. Blank lines and lines starting with `#`
//! are ignored. Every problem is reported with its line number and nothing is
//! accepted partially.

use crate::{RadiusSet, TableIdentity, vdw_radius};
use core::fmt;
use molframe_core::Element;

/// The largest radius a table may declare, in ångström; anything above is a
/// unit mistake rather than an atom.
const MAXIMUM_RADIUS: f32 = 20.0;

/// Why a radius table was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RadiusTableError {
    line: usize,
    message: String,
}

impl RadiusTableError {
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

impl fmt::Display for RadiusTableError {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(out, "radius table: {}", self.message)
        } else {
            write!(out, "radius table, line {}: {}", self.line, self.message)
        }
    }
}

impl core::error::Error for RadiusTableError {}

/// A validated, caller-supplied radius table.
#[derive(Clone, Debug, PartialEq)]
pub struct CustomRadii {
    identity: TableIdentity,
    base: Option<RadiusSet>,
    values: [Option<f32>; 119],
}

#[derive(Default)]
struct Header {
    name: Option<String>,
    version: Option<String>,
    source: Option<String>,
    licence: Option<String>,
    base: Option<RadiusSet>,
}

impl CustomRadii {
    /// Parses and validates a table.
    ///
    /// # Errors
    ///
    /// [`RadiusTableError`] for a missing or repeated header, an unknown
    /// element, a radius that is not a positive finite number, or a repeated
    /// element.
    pub fn parse(text: &str) -> Result<Self, RadiusTableError> {
        let mut header = Header::default();
        let mut values = [None; 119];
        let mut listed = 0_usize;
        for (index, raw) in text.lines().enumerate() {
            let number = index + 1;
            let line = raw.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let Some((key, value)) = line.split_once(':') {
                header.set(number, key.trim(), value.trim())?;
                continue;
            }
            read_entry(number, line, &mut values)?;
            listed += 1;
        }
        let identity = header.identity()?;
        if listed == 0 && header.base.is_none() {
            return Err(RadiusTableError::at(
                0,
                "the table lists no element and has no base",
            ));
        }
        Ok(Self {
            identity,
            base: header.base,
            values,
        })
    }

    /// The declared identity of this table.
    #[must_use]
    pub const fn identity(&self) -> &TableIdentity {
        &self.identity
    }

    /// The built-in set this table is layered over, if any.
    #[must_use]
    pub const fn base(&self) -> Option<RadiusSet> {
        self.base
    }

    /// The radius of `element`: the caller's value, else the base's, else none.
    #[must_use]
    pub fn radius(&self, element: Element) -> Option<f32> {
        let own = self
            .values
            .get(element.atomic_number() as usize)
            .copied()
            .flatten();
        match (own, self.base) {
            (Some(radius), _) => Some(radius),
            (None, Some(base)) => vdw_radius(element, base),
            (None, None) => None,
        }
    }
}

impl Header {
    fn set(&mut self, line: usize, key: &str, value: &str) -> Result<(), RadiusTableError> {
        if value.is_empty() {
            return Err(RadiusTableError::at(line, format!("`{key}` has no value")));
        }
        let slot = match key {
            "name" => &mut self.name,
            "version" => &mut self.version,
            "source" => &mut self.source,
            "licence" | "license" => &mut self.licence,
            "base" => {
                if self.base.is_some() {
                    return Err(RadiusTableError::at(line, "`base` is declared twice"));
                }
                self.base = Some(value.parse().map_err(|_| {
                    RadiusTableError::at(
                        line,
                        format!(
                            "unknown base `{value}`; expected one of {}",
                            RadiusSet::NAMES.join(", ")
                        ),
                    )
                })?);
                return Ok(());
            }
            other => {
                return Err(RadiusTableError::at(
                    line,
                    format!(
                        "unknown header `{other}`; expected name, version, source, licence or base"
                    ),
                ));
            }
        };
        if slot.is_some() {
            return Err(RadiusTableError::at(
                line,
                format!("`{key}` is declared twice"),
            ));
        }
        *slot = Some(value.to_owned());
        Ok(())
    }

    fn identity(&self) -> Result<TableIdentity, RadiusTableError> {
        let need = |field: &Option<String>, key: &str| match field {
            Some(value) => Ok(value.clone()),
            None => Err(RadiusTableError::at(
                0,
                format!("the header `{key}` is required"),
            )),
        };
        Ok(TableIdentity::supplied(
            &need(&self.name, "name")?,
            &need(&self.version, "version")?,
            &need(&self.source, "source")?,
            &need(&self.licence, "licence")?,
        ))
    }
}

fn read_entry(
    line: usize,
    text: &str,
    values: &mut [Option<f32>; 119],
) -> Result<(), RadiusTableError> {
    let mut words = text.split_whitespace();
    let (Some(symbol), Some(radius), None) = (words.next(), words.next(), words.next()) else {
        return Err(RadiusTableError::at(
            line,
            "expected `<element> <radius>` or a `key: value` header",
        ));
    };
    let Some(element) = Element::from_symbol(symbol) else {
        return Err(RadiusTableError::at(
            line,
            format!("`{symbol}` is not an element symbol"),
        ));
    };
    let radius: f32 = radius.parse().map_err(|_| {
        RadiusTableError::at(line, format!("`{radius}` is not a number of ångström"))
    })?;
    if !(radius.is_finite() && radius > 0.0 && radius <= MAXIMUM_RADIUS) {
        return Err(RadiusTableError::at(
            line,
            format!("the radius of {symbol} must be positive and at most {MAXIMUM_RADIUS} Å"),
        ));
    }
    match values.get_mut(element.atomic_number() as usize) {
        Some(slot @ None) => {
            *slot = Some(radius);
            Ok(())
        }
        Some(Some(_)) => Err(RadiusTableError::at(
            line,
            format!("{symbol} is listed twice"),
        )),
        None => Err(RadiusTableError::at(
            line,
            format!("`{symbol}` is not an element symbol"),
        )),
    }
}

#[cfg(test)]
#[path = "radius_table_tests.rs"]
mod tests;
