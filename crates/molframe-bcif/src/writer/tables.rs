//! Encodes caller-owned tables column by column into `BinaryCIF` categories.
//!
//! A table is read through [`ColumnTable`], one cell at a time. No text form
//! and no parsed document is built: each column is scanned once to choose its
//! codec and once to fill a single typed buffer, so peak workspace is one
//! column of the largest table rather than the whole table.

use super::structure::write_with_extra;
use crate::container::{EncodedCategory, EncodedColumn};
use crate::encode::StringColumnEncoder;
use crate::{EncodedData, encode_floats, encode_integers};
use molframe_cif::CifWriteOptions;
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::io::DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES;
use molframe_core::structure::Structure;
use num_traits::ToPrimitive;
use std::fmt::Write as _;

/// One borrowed cell of a [`ColumnTable`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Cell<'a> {
    /// The CIF `.` value.
    Inapplicable,
    /// The CIF `?` value.
    Unknown,
    /// Text.
    Text(&'a str),
    /// Integer.
    Integer(i64),
    /// Floating point.
    Float(f64),
}

/// A category exposed cell by cell, without materialised rows.
pub trait ColumnTable {
    /// Category name without the leading underscore.
    fn name(&self) -> &str;
    /// Number of rows.
    fn row_count(&self) -> usize;
    /// Number of columns.
    fn column_count(&self) -> usize;
    /// Item name of one column.
    fn item(&self, column: usize) -> &str;
    /// One cell; out-of-range cells read as unknown.
    fn cell(&self, column: usize, row: usize) -> Cell<'_>;
}

/// Writes a structure and appends `tables` after its coordinates.
///
/// # Errors
///
/// Returns canonical projection or `BinaryCIF` encoding diagnostics.
pub fn write_structure_with_tables(
    structure: &Structure,
    options: &CifWriteOptions,
    tables: &[&dyn ColumnTable],
) -> Result<Vec<u8>, Vec<Diagnostic>> {
    let mut extra = Vec::with_capacity(tables.len());
    for table in tables {
        extra.push(encode_table(*table).map_err(|error| vec![error])?);
    }
    let mut output = Vec::new();
    write_with_extra(
        structure,
        options,
        DEFAULT_OUTPUT_MEMORY_LIMIT_BYTES,
        extra,
        &mut output,
    )?;
    Ok(output)
}

pub(super) fn encode_table(table: &dyn ColumnTable) -> Result<EncodedCategory, Diagnostic> {
    let rows = table.row_count();
    let mut columns = Vec::with_capacity(table.column_count());
    for index in 0..table.column_count() {
        let (data, mask) = encode_column(table, index, rows)?;
        columns.push(EncodedColumn {
            name: table.item(index).to_owned(),
            data,
            mask,
        });
    }
    Ok(EncodedCategory {
        name: format!("_{}", table.name()),
        row_count: rows,
        columns,
    })
}

#[derive(Clone, Copy)]
enum Kind {
    Integer,
    Float,
    String,
}

/// What one pass over a column learns before any value is stored.
struct Scan {
    kind: Option<Kind>,
    masked: bool,
}

fn scan(table: &dyn ColumnTable, column: usize, rows: usize) -> Scan {
    let mut kind = None;
    let mut masked = false;
    for row in 0..rows {
        let candidate = match table.cell(column, row) {
            Cell::Integer(_) => Kind::Integer,
            Cell::Float(_) => Kind::Float,
            Cell::Text(_) => Kind::String,
            Cell::Inapplicable | Cell::Unknown => {
                masked = true;
                continue;
            }
        };
        kind = Some(match (kind, candidate) {
            (None, candidate) => candidate,
            (Some(Kind::Integer), Kind::Integer) => Kind::Integer,
            (Some(Kind::Integer | Kind::Float), Kind::Integer | Kind::Float) => Kind::Float,
            // A column of numbers beside text is a string column.
            _ => Kind::String,
        });
    }
    Scan { kind, masked }
}

fn mask_value(cell: Cell<'_>) -> i64 {
    match cell {
        Cell::Inapplicable => 1,
        Cell::Unknown => 2,
        _ => 0,
    }
}

fn encode_column(
    table: &dyn ColumnTable,
    column: usize,
    rows: usize,
) -> Result<(EncodedData, Option<EncodedData>), Diagnostic> {
    let scanned = scan(table, column, rows);
    let kind = match scanned.kind {
        Some(kind) => kind,
        None => Kind::String,
    };
    let data = match kind {
        Kind::Integer => {
            let mut values = Vec::with_capacity(rows);
            for row in 0..rows {
                values.push(match table.cell(column, row) {
                    Cell::Integer(value) => value,
                    _ => 0,
                });
            }
            encode_integers(&values)?
        }
        Kind::Float => {
            let mut values = Vec::with_capacity(rows);
            for row in 0..rows {
                values.push(match table.cell(column, row) {
                    Cell::Integer(value) => value.to_f64().ok_or_else(type_error)?,
                    Cell::Float(value) => value,
                    _ => 0.0,
                });
            }
            encode_floats(&values)?
        }
        Kind::String => {
            let mut encoder = StringColumnEncoder::with_capacity(rows);
            let mut digits = String::new();
            for row in 0..rows {
                match table.cell(column, row) {
                    Cell::Text(value) => encoder.push(value),
                    Cell::Integer(value) => {
                        digits.clear();
                        let _ = write!(digits, "{value}");
                        encoder.push(&digits);
                    }
                    Cell::Float(value) => {
                        digits.clear();
                        let _ = write!(digits, "{value}");
                        encoder.push(&digits);
                    }
                    Cell::Inapplicable | Cell::Unknown => encoder.push(""),
                }
            }
            encoder.finish()?
        }
    };
    let mask = if scanned.masked {
        let mut values = Vec::with_capacity(rows);
        for row in 0..rows {
            values.push(mask_value(table.cell(column, row)));
        }
        Some(encode_integers(&values)?)
    } else {
        None
    };
    Ok((data, mask))
}

fn type_error() -> Diagnostic {
    Diagnostic::new(Code::E1403).with_message("a BinaryCIF column mixes incompatible value types")
}

#[cfg(test)]
#[path = "tables_tests.rs"]
mod tests;
