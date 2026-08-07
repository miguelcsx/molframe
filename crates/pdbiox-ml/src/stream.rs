//! Owned Arrow C Stream values.

use arrow::datatypes::SchemaRef;
use arrow::error::Result;
use arrow::ffi_stream::FFI_ArrowArrayStream;
use arrow::record_batch::{RecordBatch, RecordBatchIterator};

/// One consumable Arrow C Stream Interface value.
#[derive(Debug)]
pub struct ArrowStream(FFI_ArrowArrayStream);

impl ArrowStream {
    pub(crate) fn from_batches(schema: SchemaRef, batches: Vec<RecordBatch>) -> Self {
        let reader = RecordBatchIterator::new(batches.into_iter().map(Ok), schema);
        Self(FFI_ArrowArrayStream::new(Box::new(reader)))
    }

    /// Borrows the ABI-compatible stream value.
    #[must_use]
    pub const fn as_ffi(&self) -> &FFI_ArrowArrayStream {
        &self.0
    }

    /// Consumes the wrapper and returns the ABI-compatible stream value.
    #[must_use]
    pub fn into_ffi(self) -> FFI_ArrowArrayStream {
        self.0
    }
}

pub(crate) trait ArrowTableExport {
    fn schema(&self) -> SchemaRef;
    fn record_batches(&self) -> Result<Vec<RecordBatch>>;

    fn arrow_stream(&self) -> Result<ArrowStream> {
        Ok(ArrowStream::from_batches(
            self.schema(),
            self.record_batches()?,
        ))
    }
}
