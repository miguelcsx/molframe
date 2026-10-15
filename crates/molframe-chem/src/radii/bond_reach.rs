//! How far each element reaches when bonds are inferred from coordinates.
//!
//! Two atoms are perceived as bonded when their distance lies below the sum of
//! their reaches divided by 1.95 (and above a closest plausible contact). The
//! reach of an element, and whether it is a metal that coordinates rather than
//! bonds covalently, is reference data selected by name (profile line
//! `bond-reach`, files `<name>.reach`), never compiled in:
//!
//! ```text
//! # element reach-in-angstrom [metal]
//! C 1.75
//! Zn 2.85 metal
//! ```
//!
//! An element the table does not list takes part in no inferred bond; perception
//! leaves it to deposited connectivity rather than guessing a reach.
//!
//! Cost: decoding is O(rows); a lookup is one array index.

use crate::{
    Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};
use molframe_core::Element;

const WHAT: &str = "bond reach table";
const MAXIMUM_REACH: f32 = 20.0;

/// A validated bonding reach per element.
#[derive(Clone, Debug, PartialEq)]
pub struct BondReach {
    reach: [Option<(f32, bool)>; 119],
}

impl Default for BondReach {
    fn default() -> Self {
        Self { reach: [None; 119] }
    }
}

impl TableKind for BondReach {
    const KIND: &'static str = "bond-reach";
    const EXTENSION: &'static str = "reach";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(WHAT, text_of(WHAT, bytes)?, &mut table)?;
        if table.reach.iter().all(Option::is_none) {
            return Err(TableTextError::at(WHAT, 0, "the table lists no element"));
        }
        Ok(Provided { identity, table })
    }
}

impl TableLines for BondReach {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let row = Row::new(WHAT, line);
        let words: Vec<&str> = text.split_whitespace().collect();
        let (symbol, value, metal) = match words[..] {
            [symbol, value] => (symbol, value, false),
            [symbol, value, "metal"] => (symbol, value, true),
            _ => return Err(row.fail("expected `element reach [metal]`")),
        };
        let element = row.element(symbol)?;
        let reach = row.bounded::<f32>(
            value,
            f64::from(MAXIMUM_REACH),
            &format!("the reach of {symbol}"),
        )?;
        row.claim(
            self.reach.get_mut(element.atomic_number() as usize),
            symbol,
            (reach, metal),
        )
    }
}

impl BondReach {
    /// The reach of `element` in ångström, or `None` when the table lists none.
    #[must_use]
    pub fn reach(&self, element: Element) -> Option<f32> {
        self.entry(element).map(|(reach, _)| reach)
    }

    /// Whether the table marks `element` as a metal.
    #[must_use]
    pub fn is_metal(&self, element: Element) -> bool {
        self.entry(element).is_some_and(|(_, metal)| metal)
    }

    fn entry(&self, element: Element) -> Option<(f32, bool)> {
        self.reach
            .get(element.atomic_number() as usize)
            .copied()
            .flatten()
    }
}

#[cfg(test)]
#[path = "bond_reach_tests.rs"]
mod tests;
