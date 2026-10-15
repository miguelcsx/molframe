//! A table of maximum residue areas, supplied by the caller.
//!
//! Relative accessibility is an area divided by a maximum, and which maximum is
//! a choice. [`MaxAreaTable`] holds one set of maxima; the identity that says
//! where they came from travels with it in [`molframe_chem::Provided`], and the
//! table is selected by name like every other table. Its text is the shared
//! table format with one row per residue:
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
//! Every row is a residue name and a positive finite area. A residue the table
//! lacks has no relative accessibility. Protonation variants (`HID`, `CYX`, ...)
//! take their parent's value. Every problem is reported with its line number and
//! nothing is accepted partially.
//!
//! Complexity: decoding is O(rows); a lookup is O(log rows).

use crate::protonation::parent_residue;
use molframe_chem::{
    Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};
use std::collections::BTreeMap;

/// The largest area a table may declare, in square ångström; anything above is
/// a unit mistake.
const MAXIMUM_AREA: f64 = 2000.0;

const WHAT: &str = "maximum-area table";

/// Maximum accessible area per residue name.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct MaxAreaTable {
    areas: BTreeMap<String, f64>,
}

impl TableKind for MaxAreaTable {
    const KIND: &'static str = "max-area";
    const EXTENSION: &'static str = "maxarea";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(WHAT, text_of(WHAT, bytes)?, &mut table)?;
        if table.areas.is_empty() {
            return Err(TableTextError::at(WHAT, 0, "the table lists no residue"));
        }
        Ok(Provided { identity, table })
    }
}

impl TableLines for MaxAreaTable {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let row = Row::new(WHAT, line);
        let [name, value] = row.words(text, "`residue area`, two words")?;
        let area = row.bounded::<f64>(value, MAXIMUM_AREA, &format!("area {value}"))?;
        if self.areas.insert(name.to_ascii_uppercase(), area).is_some() {
            return Err(row.fail(format!("residue `{name}` is listed twice")));
        }
        Ok(())
    }
}

impl MaxAreaTable {
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

#[cfg(test)]
#[path = "max_area_table_tests.rs"]
mod tests;
