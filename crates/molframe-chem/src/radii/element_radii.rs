//! One radius per element, read from a table file: van der Waals or covalent.
//!
//! Which radii a calculation uses is a scientific choice (Bondi, Alvarez, a force
//! field's, a lab's own), so no radius lives in code. A table is plain text in the
//! shared format, selected by name from the cache or the profile like every other
//! table, and it declares its own identity so a result says which numbers produced it:
//!
//! ```text
//! name: lab-radii
//! version: 2
//! source: Smith lab, internal survey of the CSD
//! licence: CC0-1.0
//!
//! # element radius-in-angstrom
//! C 1.75
//! N 1.60
//! ```
//!
//! The same reader serves two kinds through the [`RadiusRole`] marker: [`VdwRadii`]
//! (profile line `vdw-radii`, files `<name>.radii`) and [`CovalentRadii`] (profile line
//! `covalent-radii`, files `<name>.covalent`). An element a table does not list has no
//! radius; the library never fills one table from another. Every row is an element
//! symbol and a positive finite radius at most 20 Å, every problem is reported with
//! its line number, and a table is accepted whole or not at all.
//!
//! Cost: decoding is O(rows); a lookup is one array index and allocates nothing.

use crate::{
    Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};
use core::marker::PhantomData;
use molframe_core::Element;

/// The largest radius a table may declare, in ångström; anything above is a
/// unit mistake rather than an atom.
const MAXIMUM_RADIUS: f32 = 20.0;

/// What a radius table measures, and how its kind is named.
pub trait RadiusRole {
    /// Kind name, also the profile line that selects a table of this role.
    const KIND: &'static str;
    /// File extension of a cache entry of this role.
    const EXTENSION: &'static str;
    /// How a table of this role is called in messages.
    const WHAT: &'static str;
}

/// The van der Waals role: the size of an atom as an obstacle.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Vdw {}

/// The single-bond covalent role: the size of an atom as a bonding partner.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Covalent {}

impl RadiusRole for Vdw {
    const KIND: &'static str = "vdw-radii";
    const EXTENSION: &'static str = "radii";
    const WHAT: &'static str = "van der Waals radius table";
}

impl RadiusRole for Covalent {
    const KIND: &'static str = "covalent-radii";
    const EXTENSION: &'static str = "covalent";
    const WHAT: &'static str = "covalent radius table";
}

/// A validated radius per element for the role `R`.
#[derive(Clone, Debug, PartialEq)]
pub struct ElementRadii<R> {
    values: [Option<f32>; 119],
    role: PhantomData<fn() -> R>,
}

/// Van der Waals radii.
pub type VdwRadii = ElementRadii<Vdw>;

/// Single-bond covalent radii.
pub type CovalentRadii = ElementRadii<Covalent>;

impl<R> Default for ElementRadii<R> {
    fn default() -> Self {
        Self {
            values: [None; 119],
            role: PhantomData,
        }
    }
}

impl<R: RadiusRole> TableKind for ElementRadii<R> {
    const KIND: &'static str = R::KIND;
    const EXTENSION: &'static str = R::EXTENSION;

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(R::WHAT, text_of(R::WHAT, bytes)?, &mut table)?;
        if table.values.iter().all(Option::is_none) {
            return Err(TableTextError::at(R::WHAT, 0, "the table lists no element"));
        }
        Ok(Provided { identity, table })
    }
}

impl<R: RadiusRole> TableLines for ElementRadii<R> {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let row = Row::new(R::WHAT, line);
        let [symbol, radius] = row.words(text, "`<element> <radius>` or a `key: value` header")?;
        let element = row.element(symbol)?;
        let radius = row.bounded::<f32>(
            radius,
            f64::from(MAXIMUM_RADIUS),
            &format!("the radius of {symbol}"),
        )?;
        row.claim(
            self.values.get_mut(element.atomic_number() as usize),
            symbol,
            radius,
        )
    }
}

impl<R> ElementRadii<R> {
    /// The radius of `element` in ångström, or `None` when the table lists none.
    #[must_use]
    pub fn radius(&self, element: Element) -> Option<f32> {
        self.values
            .get(element.atomic_number() as usize)
            .copied()
            .flatten()
    }
}

#[cfg(test)]
#[path = "element_radii_tests.rs"]
mod tests;
