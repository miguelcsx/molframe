//! Ionic radii indexed by element, charge and coordination, read from a table file.
//!
//! The catalogue is the Shannon-style list of effective and crystal radii per
//! element, charge, coordination and spin state. It is a table like any other:
//! selected by name (profile line `ionic-radii`, files `<name>.ionic`), with its
//! identity recorded in results, never compiled in. One row per entry:
//!
//! ```text
//! # element charge coordination spin(-|HS|LS) effective crystal reliable(yes|no|-)
//! Fe 2 VI HS 0.78 0.92 yes
//! ```
//!
//! Rows must be grouped by element in atomic-number order. A missing element has
//! no entries and the library never substitutes an estimate.
//!
//! Cost: decoding is O(rows); a lookup is one array index and a slice.

use crate::{
    Provided, Row, TableIdentity, TableKind, TableLines, TableTextError, read_table, text_of,
};
use molframe_core::Element;
use std::ops::Range;

const WHAT: &str = "ionic radius table";

/// Spin state attached to an ionic-radius entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IonicSpin {
    /// No spin distinction applies or was supplied.
    Unspecified,
    /// High-spin configuration.
    High,
    /// Low-spin configuration.
    Low,
}

impl IonicSpin {
    /// Reads the word a table spells a spin with: `-`, `HS` or `LS`.
    ///
    /// # Errors
    ///
    /// A message naming the word when it is none of those.
    pub fn from_word(word: &str) -> Result<Self, String> {
        match word {
            "-" => Ok(Self::Unspecified),
            "HS" => Ok(Self::High),
            "LS" => Ok(Self::Low),
            other => Err(format!("spin `{other}` must be `-`, `HS` or `LS`")),
        }
    }
}

/// One charge- and coordination-specific ionic radius.
#[derive(Clone, Debug, PartialEq)]
pub struct IonicRadius {
    /// Formal ionic charge.
    pub charge: i8,
    /// Coordination geometry/count in Shannon's notation (`VI`, `IVSQ`, ...).
    pub coordination: Box<str>,
    /// Spin state where the table distinguishes it.
    pub spin: IonicSpin,
    /// Effective ionic radius in ångström.
    pub ionic_radius: f32,
    /// Crystal radius in ångström.
    pub crystal_radius: f32,
    /// Source marks this entry as among the most reliable, when specified.
    pub most_reliable: Option<bool>,
}

/// A validated ionic radius catalogue.
#[derive(Clone, Debug, PartialEq)]
pub struct IonicRadii {
    records: Vec<IonicRadius>,
    ranges: Vec<Range<usize>>,
    last: usize,
}

impl Default for IonicRadii {
    fn default() -> Self {
        Self {
            records: Vec::new(),
            ranges: vec![0..0; 119],
            last: 0,
        }
    }
}

impl TableKind for IonicRadii {
    const KIND: &'static str = "ionic-radii";
    const EXTENSION: &'static str = "ionic";

    fn decode(bytes: &[u8], _: TableIdentity) -> Result<Provided<Self>, TableTextError> {
        let mut table = Self::default();
        let identity = read_table(WHAT, text_of(WHAT, bytes)?, &mut table)?;
        if table.records.is_empty() {
            return Err(TableTextError::at(WHAT, 0, "the table lists no entry"));
        }
        Ok(Provided { identity, table })
    }
}

impl TableLines for IonicRadii {
    fn header(&mut self, _: usize, _: &str, _: &str) -> Result<bool, TableTextError> {
        Ok(false)
    }

    fn row(&mut self, line: usize, text: &str) -> Result<(), TableTextError> {
        let row = Row::new(WHAT, line);
        let [
            symbol,
            charge,
            coordination,
            spin,
            effective,
            crystal,
            reliable,
        ] = row.words(
            text,
            "`element charge coordination spin effective crystal reliable`",
        )?;
        let element = row.element(symbol)?;
        let z = element.atomic_number() as usize;
        if z < self.last {
            return Err(row.fail(format!(
                "`{symbol}` is out of atomic-number order; group rows by element"
            )));
        }
        self.last = z;
        let fail = |message: String| row.fail(message);
        let number = |word: &str| row.finite::<f32>(word);
        let record = IonicRadius {
            charge: charge
                .parse()
                .map_err(|_| fail(format!("`{charge}` is not a charge")))?,
            coordination: coordination.into(),
            spin: IonicSpin::from_word(spin).map_err(fail)?,
            ionic_radius: number(effective)?,
            crystal_radius: number(crystal)?,
            most_reliable: match reliable {
                "yes" => Some(true),
                "no" => Some(false),
                "-" => None,
                other => return Err(fail(format!("`{other}` must be `yes`, `no` or `-`"))),
            },
        };
        let end = self.records.len() + 1;
        let range = self
            .ranges
            .get_mut(z)
            .ok_or_else(|| fail(format!("`{symbol}` is not an element symbol")))?;
        if *range == (0..0) {
            range.start = end - 1;
        }
        range.end = end;
        self.records.push(record);
        Ok(())
    }
}

impl IonicRadii {
    /// Every entry for `element`, in table order; empty when it has none.
    #[must_use]
    pub fn entries(&self, element: Element) -> &[IonicRadius] {
        match self
            .ranges
            .get(element.atomic_number() as usize)
            .and_then(|range| self.records.get(range.clone()))
        {
            Some(found) => found,
            None => &[],
        }
    }
}

#[cfg(test)]
#[path = "ionic_tests.rs"]
mod tests;
