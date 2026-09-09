//! Incremental PDB output over the shared bounded text adapter.

use super::pdb::{PdbOptions, render_selected};
use pdbiox_core::diagnostic::{Code, Diagnostic};
use pdbiox_core::io::{Select, SelectAll, TextOutput};
use pdbiox_core::structure::Structure;
use std::io::Write;

/// Streams a complete PDB directly to a byte destination.
///
/// # Errors
///
/// Returns every fixed-column refusal before the first byte, or an output
/// diagnostic if the destination stops accepting bytes.
pub fn write_to<W: Write>(
    structure: &Structure,
    options: &PdbOptions,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    write_selected_to(structure, options, &SelectAll, output)
}

/// Streams selected PDB records without materialising the complete text.
///
/// # Errors
///
/// Returns capacity, required-field, or destination diagnostics.
pub fn write_selected_to<W: Write>(
    structure: &Structure,
    options: &PdbOptions,
    select: &impl Select,
    output: &mut W,
) -> Result<(), Vec<Diagnostic>> {
    let mut output = TextOutput::new(output);
    render_selected(&mut output, structure, options, select)?;
    output
        .finish()
        .map_err(|error| vec![output_diagnostic(error)])
}

fn output_diagnostic(error: impl std::fmt::Display) -> Diagnostic {
    Diagnostic::new(Code::E7901).with_context("reason", error.to_string())
}
