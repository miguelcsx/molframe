//! Canonical `ModelCIF` projection over the shared streaming mmCIF renderer.

use crate::{ModelCategory, ModelCif, ValueRef};
use molframe_cif::{CifWriteOptions, CifWriteToError, write_quoted};
use molframe_core::Structure;
use molframe_core::io::TextOutput;
use std::fmt;
use std::io::{self, Write as IoWrite};

/// Writes canonical coordinates and retained `ModelCIF` extension categories.
///
/// # Errors
///
/// Returns the mmCIF preflight error without returning partial output.
pub fn write_canonical(
    structure: &Structure,
    model_cif: &ModelCif,
) -> Result<String, CifWriteToError> {
    write_canonical_with_options(structure, model_cif, &CifWriteOptions::new())
}

/// Writes canonical coordinates and retained extension categories with explicit options.
///
/// # Errors
///
/// Returns the mmCIF preflight error without returning partial output.
pub fn write_canonical_with_options(
    structure: &Structure,
    model_cif: &ModelCif,
    options: &CifWriteOptions,
) -> Result<String, CifWriteToError> {
    let mut output = Vec::with_capacity(structure.atom_count() as usize * 100);
    write_canonical_to(structure, model_cif, options, &mut output)?;
    String::from_utf8(output)
        .map_err(|error| CifWriteToError::Output(io::Error::new(io::ErrorKind::InvalidData, error)))
}

/// Streams canonical coordinates and compact `ModelCIF` categories.
///
/// Structural preflight completes before the first byte is written. Extension
/// values are visited directly from compact columns without materialising rows
/// or a full output buffer.
///
/// # Errors
///
/// Returns the shared mmCIF preflight refusal or the destination I/O error.
pub fn write_canonical_to<W: IoWrite>(
    structure: &Structure,
    model_cif: &ModelCif,
    options: &CifWriteOptions,
    output: &mut W,
) -> Result<(), CifWriteToError> {
    molframe_cif::write_canonical_to(structure, options, output)?;
    let mut output = TextOutput::new(output);
    for category in &model_cif.categories {
        write_category(&mut output, category);
    }
    output.finish().map_err(CifWriteToError::Output)
}

fn write_category(output: &mut impl fmt::Write, category: &ModelCategory) {
    if category.row_count() == 0 {
        return;
    }
    if category.row_count() == 1 {
        for (column, item) in category.items().iter().enumerate() {
            let _ = write!(output, "_{}.{item} ", category.name());
            write_value(output, category.value_ref(column, 0));
            let _ = output.write_char('\n');
        }
        let _ = output.write_str("#\n");
        return;
    }
    let _ = output.write_str("loop_\n");
    for item in category.items() {
        let _ = writeln!(output, "_{}.{item}", category.name());
    }
    for row in category.rows() {
        for index in 0..category.items().len() {
            if index > 0 {
                let _ = output.write_char(' ');
            }
            write_value(output, row.value_ref(index));
        }
        let _ = output.write_char('\n');
    }
    let _ = output.write_str("#\n");
}

fn write_value(output: &mut impl fmt::Write, value: Option<ValueRef<'_>>) {
    match value {
        Some(ValueRef::Inapplicable) => {
            let _ = output.write_char('.');
        }
        Some(ValueRef::Unknown) | None => {
            let _ = output.write_char('?');
        }
        Some(ValueRef::Text(value)) => write_quoted(output, value),
        Some(ValueRef::Integer(value)) => {
            let _ = write!(output, "{value}");
        }
        Some(ValueRef::Float(value)) => {
            let _ = write!(output, "{value}");
        }
    }
}

#[cfg(test)]
#[path = "canonical_tests.rs"]
mod tests;
