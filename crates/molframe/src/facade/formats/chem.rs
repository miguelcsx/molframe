//! Chemical Component Dictionary reading.

use molframe_core::contract::DictionaryVersion;
use molframe_core::diagnostic::{Diagnostic, Findings};
use molframe_core::io::{InputBuffer, Limits};
use std::path::Path;

/// Reads a Chemical Component Dictionary file as a versioned provider.
///
/// This explicit entry point avoids guessing whether a `.cif` file is a
/// structure or a component dictionary. Parsing and component lowering remain
/// in their owning Rust crates.
///
/// # Errors
///
/// Returns file, CIF syntax, or component-definition diagnostics.
pub fn read_component_dictionary(
    path: impl AsRef<Path>,
    version: DictionaryVersion,
) -> Result<(molframe_chem::CifProvider, Vec<Diagnostic>), Findings> {
    let input = InputBuffer::open(path, Limits::default()).map_err(Findings::from)?;
    molframe_chem::read_ccd(&input, version).map_err(Findings::from)
}
