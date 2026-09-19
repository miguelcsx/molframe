//! Arrow IPC serialisation for the atom table.
//!
//! Each atom chunk is written before the next one is materialised, so peak
//! memory is one record batch rather than the whole table.

use super::AtomTable;
use super::table_file::{TableFileError, schema_with_metadata, visit_batches};
use arrow::ipc::writer::FileWriter;
use molframe_core::Structure;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

/// Writes a structure atom table as an Arrow IPC file.
///
/// # Errors
///
/// Returns [`TableFileError`] when the destination, atom table or IPC stream
/// cannot be created or completed.
pub fn write_atom_ipc(path: impl AsRef<Path>, structure: &Structure) -> Result<(), TableFileError> {
    write_atom_ipc_with_metadata(path, structure, BTreeMap::new())
}

/// Writes an Arrow IPC atom table with caller-supplied schema metadata.
///
/// # Errors
///
/// Returns [`TableFileError`] for destination, schema, batch or stream errors.
pub fn write_atom_ipc_with_metadata(
    path: impl AsRef<Path>,
    structure: &Structure,
    metadata: BTreeMap<String, String>,
) -> Result<(), TableFileError> {
    let table = AtomTable::new(structure);
    let schema = schema_with_metadata(&table, metadata);
    let file = File::create(path)?;
    let mut writer = FileWriter::try_new(file, &schema)?;
    visit_batches(&table, &schema, |batch| {
        writer.write(batch)?;
        Ok(())
    })?;
    writer.finish()?;
    Ok(())
}
