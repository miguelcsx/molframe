//! Arrow export for native contact tables.

use super::extension::{ExportCost, field};
use super::stream::{ArrowStream, ArrowTableExport};
use arrow::array::{ArrayRef, Float32Array, UInt32Array};
use arrow::datatypes::{DataType, Schema, SchemaRef};
use arrow::error::{ArrowError, Result};
use arrow::record_batch::RecordBatch;
use molframe_analysis::ContactTable;
use std::sync::Arc;

/// Lazy Arrow view over an owner-backed native contact table.
#[derive(Clone, Debug)]
pub struct ContactArrowTable {
    contacts: Arc<ContactTable>,
    schema: SchemaRef,
}

impl ContactArrowTable {
    /// Retains an existing shared table without copying its native columns.
    #[must_use]
    pub fn new(contacts: Arc<ContactTable>) -> Self {
        let fields = vec![
            field("first", DataType::UInt32, false, None, ExportCost::Copy),
            field("second", DataType::UInt32, false, None, ExportCost::Copy),
            field("distance", DataType::Float32, false, None, ExportCost::Copy),
        ];
        Self {
            contacts,
            schema: Arc::new(Schema::new(fields)),
        }
    }

    /// Arrow schema, including the explicit per-column copy cost.
    #[must_use]
    pub fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    /// Materializes the aligned columns as one Arrow record batch.
    ///
    /// # Errors
    ///
    /// Returns an Arrow error if the native column lengths are inconsistent.
    pub fn record_batches(&self) -> Result<Vec<RecordBatch>> {
        <Self as ArrowTableExport>::record_batches(self)
    }

    /// Creates a lazy Arrow C stream retaining the native table owner.
    ///
    /// # Errors
    ///
    /// Reserved for failures constructing the stream adapter.
    pub fn arrow_stream(&self) -> Result<ArrowStream> {
        <Self as ArrowTableExport>::arrow_stream(self)
    }
}

impl ArrowTableExport for ContactArrowTable {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn batch_count(&self) -> usize {
        usize::from(!self.contacts.is_empty())
    }

    fn batch(&self, index: usize) -> Result<RecordBatch> {
        if index != 0 || self.contacts.is_empty() {
            return Err(ArrowError::InvalidArgumentError(
                "contact table batch index is out of bounds".to_owned(),
            ));
        }
        let first =
            UInt32Array::from_iter_values(self.contacts.first().iter().map(|value| value.get()));
        let second =
            UInt32Array::from_iter_values(self.contacts.second().iter().map(|value| value.get()));
        let distance = Float32Array::from(self.contacts.distances().to_vec());
        let columns: Vec<ArrayRef> = vec![Arc::new(first), Arc::new(second), Arc::new(distance)];
        RecordBatch::try_new(self.schema.clone(), columns)
    }
}

#[cfg(test)]
#[path = "contact_tests.rs"]
mod tests;
