//! Standard file serialisation for the Arrow atom table.

use crate::AtomTable;
use arrow::error::ArrowError;
use arrow::ipc::writer::FileWriter;
use parquet::arrow::ArrowWriter;
use pdbiox_core::Structure;
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
    let (schema, batches) = batches_with_metadata(&table, metadata)?;
    let file = File::create(path)?;
    let mut writer = FileWriter::try_new(file, &schema)?;
    for batch in &batches {
        writer.write(batch)?;
    }
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
    let (schema, batches) = batches_with_metadata(&table, metadata)?;
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, None)?;
    for batch in &batches {
        writer.write(batch)?;
    }
    let _ = writer.close()?;
    Ok(())
}

fn batches_with_metadata(
    table: &AtomTable,
    metadata: BTreeMap<String, String>,
) -> Result<
    (
        arrow::datatypes::SchemaRef,
        Vec<arrow::record_batch::RecordBatch>,
    ),
    ArrowError,
> {
    let base = table.schema();
    let schema = Arc::new(arrow::datatypes::Schema::new_with_metadata(
        base.fields().clone(),
        metadata.into_iter().collect(),
    ));
    let batches = table
        .record_batches()?
        .into_iter()
        .map(|batch| {
            arrow::record_batch::RecordBatch::try_new(schema.clone(), batch.columns().to_vec())
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((schema, batches))
}

#[cfg(test)]
#[path = "files_tests.rs"]
mod tests;
