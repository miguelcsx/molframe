//! Canonical `ModelCIF` projection over the shared mmCIF renderer.

use crate::{ModelCategory, ModelCif};
use pdbiox_cif::{CifWriteError, CifWriteOptions, render_value};
use pdbiox_core::Structure;
use std::fmt::Write as _;

/// Writes canonical coordinates and retained `ModelCIF` extension categories.
///
/// # Errors
///
/// Returns the mmCIF preflight error without returning partial output.
pub fn write_canonical(
    structure: &Structure,
    model_cif: &ModelCif,
) -> Result<String, CifWriteError> {
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
) -> Result<String, CifWriteError> {
    let mut output = pdbiox_cif::write_canonical_with_options(structure, options)?;
    for category in &model_cif.categories {
        write_category(&mut output, category);
    }
    Ok(output)
}

fn write_category(output: &mut String, category: &ModelCategory) {
    if category.rows().is_empty() {
        return;
    }
    if category.rows().len() == 1 {
        for (item, value) in category.items().iter().zip(category.rows()[0].values()) {
            let _ = writeln!(
                output,
                "_{}.{item} {}",
                category.name(),
                render_value(value)
            );
        }
        output.push_str("#\n");
        return;
    }
    output.push_str("loop_\n");
    for item in category.items() {
        let _ = writeln!(output, "_{}.{item}", category.name());
    }
    for row in category.rows() {
        for (index, value) in row.values().iter().enumerate() {
            if index > 0 {
                output.push(' ');
            }
            output.push_str(&render_value(value));
        }
        output.push('\n');
    }
    output.push_str("#\n");
}

#[cfg(test)]
#[path = "canonical_tests.rs"]
mod tests;
