//! Shared plumbing for Arrow-backed atom-table file writers.

use super::AtomTable;
use super::stream::ArrowTableExport;
use arrow::datatypes::{Schema, SchemaRef};
use arrow::error::ArrowError;
use arrow::record_batch::RecordBatch;
use std::collections::BTreeMap;
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

molframe_core::diagnostic_from!(TableFileError, |_error| molframe_core::Code::E7901);

pub(crate) fn schema_with_metadata(
    table: &AtomTable,
    metadata: BTreeMap<String, String>,
) -> SchemaRef {
    let base = table.schema();
    // Table provenance, units and schema version first; caller keys win.
    let mut merged = base.metadata().clone();
    merged.extend(metadata);
    Arc::new(Schema::new_with_metadata(base.fields().clone(), merged))
}

/// Writes each table chunk before materialising the next one.
pub(crate) fn visit_batches(
    table: &AtomTable,
    schema: &SchemaRef,
    mut visit: impl FnMut(&RecordBatch) -> Result<(), TableFileError>,
) -> Result<(), TableFileError> {
    for index in 0..table.batch_count() {
        // Rebind rather than `with_schema`: the file schema may override a
        // table metadata key, which `with_schema` rejects.
        let batch = table.batch(index)?;
        let batch = RecordBatch::try_new(schema.clone(), batch.columns().to_vec())?;
        visit(&batch)?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "table_file_tests.rs"]
mod tests;
