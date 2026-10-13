//! What every caller-supplied text table shares: the error and the identity header.
//!
//! A table file opens with `name:`, `version:`, `source:` and `licence:` lines
//! that become its [`TableIdentity`]; the rows after them belong to the table's
//! own format. Every problem is reported with its one-based line number, and a
//! table is accepted whole or not at all.

use core::fmt;
use molframe_chem::TableIdentity;

/// Why a text table was refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TableTextError {
    what: &'static str,
    line: usize,
    message: String,
}

impl TableTextError {
    pub(crate) fn at(what: &'static str, line: usize, message: impl Into<String>) -> Self {
        Self {
            what,
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

impl fmt::Display for TableTextError {
    fn fmt(&self, out: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.line == 0 {
            write!(out, "{}: {}", self.what, self.message)
        } else {
            write!(out, "{}, line {}: {}", self.what, self.line, self.message)
        }
    }
}

impl core::error::Error for TableTextError {}

/// The four identity headers of a table, as read.
#[derive(Default)]
pub(crate) struct IdentityFields {
    fields: [Option<String>; 4],
}

impl IdentityFields {
    /// Reads an identity header; `Ok(false)` when `key` is not one of them.
    pub(crate) fn set(
        &mut self,
        what: &'static str,
        line: usize,
        key: &str,
        value: &str,
    ) -> Result<bool, TableTextError> {
        let slot = match key {
            "name" => 0,
            "version" => 1,
            "source" => 2,
            "licence" | "license" => 3,
            _ => return Ok(false),
        };
        if value.is_empty() {
            return Err(TableTextError::at(
                what,
                line,
                format!("`{key}` has no value"),
            ));
        }
        match self.fields.get_mut(slot) {
            Some(field) if field.is_none() => *field = Some(value.to_owned()),
            _ => {
                return Err(TableTextError::at(
                    what,
                    line,
                    format!("`{key}` is declared twice"),
                ));
            }
        }
        Ok(true)
    }

    /// The declared identity.
    pub(crate) fn finish(self, what: &'static str) -> Result<TableIdentity, TableTextError> {
        let [name, version, source, licence] = self.fields;
        let need = |value: Option<String>, key: &str| {
            value.ok_or_else(|| TableTextError::at(what, 0, format!("missing `{key}` header")))
        };
        Ok(TableIdentity::supplied(
            &need(name, "name")?,
            &need(version, "version")?,
            &need(source, "source")?,
            &need(licence, "licence")?,
        ))
    }
}
