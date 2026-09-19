//! Standard file serialisation for the Arrow atom table.

use crate::AtomTable;
use crate::stream::ArrowTableExport;
use arrow::datatypes::{Schema, SchemaRef};
use arrow::error::ArrowError;
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;
use molframe_core::Structure;
use parquet::arrow::ArrowWriter;
use std::collections::BTreeMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

/// Failure while writing an Arrow-backed table file.
#[derive(Debug, thiserror::Error)]
pub enum TableFileError {
    /// Destination could not be created or written.
    #[error("table file I/O failed: {0}")]
    Io(#[from] std::io::Error),
    /// Arrow schema or batch construction failed.
    #[error("Arrow export failed: {0}")]
    Arrow(#[from] ArrowError),
    /// Parquet encoding failed.
    #[error("Parquet export failed: {0}")]
    Parquet(#[from] parquet::errors::ParquetError),
}

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
    let _ = writer.close()?;
    Ok(())
}

fn schema_with_metadata(table: &AtomTable, metadata: BTreeMap<String, String>) -> SchemaRef {
    let base = table.schema();
    Arc::new(Schema::new_with_metadata(
        base.fields().clone(),
        metadata.into_iter().collect(),
    ))
}

/// Writes each table chunk before materialising the next one.
fn visit_batches(
    table: &AtomTable,
    schema: &SchemaRef,
    mut visit: impl FnMut(&RecordBatch) -> Result<(), TableFileError>,
) -> Result<(), TableFileError> {
    for index in 0..table.batch_count() {
        let batch = table.batch(index)?.with_schema(schema.clone())?;
        visit(&batch)?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "files_tests.rs"]
mod tests;
