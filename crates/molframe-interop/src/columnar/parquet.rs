//! Parquet serialisation for the atom table.
//!
//! Each atom chunk is encoded before the next one is materialised, so peak
//! memory is one record batch rather than the whole table.

use super::AtomTable;
use super::table_file::{TableFileError, schema_with_metadata, visit_batches};
use molframe_core::Structure;
use parquet::arrow::ArrowWriter;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;

/// Writes a structure atom table as a Parquet file.
///
/// # Errors
///
/// Returns [`TableFileError`] when the destination, atom table or Parquet file
/// cannot be created or completed.
pub fn write_atom_parquet(
    path: impl AsRef<Path>,
    structure: &Structure,
) -> Result<(), TableFileError> {
    write_atom_parquet_with_metadata(path, structure, BTreeMap::new())
}

/// Writes a Parquet atom table with caller-supplied Arrow schema metadata.
///
/// # Errors
///
/// Returns [`TableFileError`] for destination, schema, batch or encoding errors.
pub fn write_atom_parquet_with_metadata(
    path: impl AsRef<Path>,
    structure: &Structure,
    metadata: BTreeMap<String, String>,
) -> Result<(), TableFileError> {
    let table = AtomTable::new(structure);
    let schema = schema_with_metadata(&table, metadata);
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema.clone(), None)?;
    visit_batches(&table, &schema, |batch| {
        writer.write(batch)?;
        Ok(())
    })?;
    writer.close()?;
    Ok(())
}
