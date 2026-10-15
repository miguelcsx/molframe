//! The one reader for the rows of a reference table.
//!
//! Every table kind reads its rows the same way: split the words, parse a finite
//! number, check a bound, resolve an element symbol, refuse a key given twice.
//! [`Row`] carries the table's name and the one-based line, so each refusal is a
//! [`TableTextError`] that names both; a kind keeps only what is its own, the
//! layout of its columns and what it stores.
//!
//! Cost: constant per row; no allocation unless a refusal is built.

use crate::TableTextError;
use core::str::FromStr;
use molframe_core::Element;

/// One row of a table, for reading its words and refusing it with its line.
#[derive(Clone, Copy, Debug)]
pub struct Row {
    what: &'static str,
    line: usize,
}

impl Row {
    /// The row on one-based `line` of a table called `what`.
    #[must_use]
    pub const fn new(what: &'static str, line: usize) -> Self {
        Self { what, line }
    }

    /// A refusal of this row.
    #[must_use]
    pub fn fail(self, message: impl Into<String>) -> TableTextError {
        TableTextError::at(self.what, self.line, message)
    }

    /// The row's words when there are exactly `N`, else a refusal that says `layout`.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when the row does not have `N` words.
    pub fn words<'a, const N: usize>(
        self,
        text: &'a str,
        layout: &str,
    ) -> Result<[&'a str; N], TableTextError> {
        let mut words = text.split_whitespace();
        let mut found = [""; N];
        for slot in &mut found {
            *slot = words
                .next()
                .ok_or_else(|| self.fail(format!("expected {layout}")))?;
        }
        if words.next().is_some() {
            return Err(self.fail(format!("expected {layout}")));
        }
        Ok(found)
    }

    /// The element a symbol names.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when `symbol` is not an element symbol.
    pub fn element(self, symbol: &str) -> Result<Element, TableTextError> {
        Element::from_symbol(symbol)
            .ok_or_else(|| self.fail(format!("`{symbol}` is not an element symbol")))
    }

    /// A finite number.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when `word` is not a finite number of type `T`.
    pub fn finite<T: FromStr + Into<f64> + Copy>(self, word: &str) -> Result<T, TableTextError> {
        word.parse::<T>()
            .ok()
            .filter(|value| (*value).into().is_finite())
            .ok_or_else(|| self.fail(format!("`{word}` is not a number")))
    }

    /// A number above zero and at most `limit`; `noun` names it in the refusal.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when `word` is not a finite number in `(0, limit]`.
    pub fn bounded<T: FromStr + Into<f64> + Copy>(
        self,
        word: &str,
        limit: f64,
        noun: &str,
    ) -> Result<T, TableTextError> {
        let value = self.finite::<T>(word)?;
        let number: f64 = value.into();
        if number > 0.0 && number <= limit {
            Ok(value)
        } else {
            Err(self.fail(format!("{noun} must be positive and at most {limit}")))
        }
    }

    /// Stores `value` in `slot`, refusing a `key` that was already given.
    ///
    /// # Errors
    ///
    /// [`TableTextError`] when `slot` is absent or already holds a value.
    pub fn claim<T>(
        self,
        slot: Option<&mut Option<T>>,
        key: &str,
        value: T,
    ) -> Result<(), TableTextError> {
        match slot {
            Some(slot @ None) => {
                *slot = Some(value);
                Ok(())
            }
            _ => Err(self.fail(format!("`{key}` is listed twice"))),
        }
    }
}

#[cfg(test)]
#[path = "table_row_tests.rs"]
mod tests;
