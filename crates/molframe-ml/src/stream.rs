//! Owned Arrow C Stream values.

use arrow::datatypes::SchemaRef;
use arrow::error::Result;
use arrow::ffi_stream::FFI_ArrowArrayStream;
use arrow::record_batch::{RecordBatch, RecordBatchReader};

/// One consumable Arrow C Stream Interface value.
///
/// The stream owns its pull-based producer. Each array handed to a consumer
/// owns its backing buffers independently and may outlive this value.
#[derive(Debug)]
pub struct ArrowStream(FFI_ArrowArrayStream);

impl ArrowStream {
    fn from_reader(reader: impl RecordBatchReader + Send + 'static) -> Self {
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

pub(crate) trait ArrowTableExport: Clone + Send + 'static {
    fn schema(&self) -> SchemaRef;
    fn batch_count(&self) -> usize;
    fn batch(&self, index: usize) -> Result<RecordBatch>;

    fn record_batches(&self) -> Result<Vec<RecordBatch>> {
        (0..self.batch_count())
            .map(|index| self.batch(index))
            .collect()
    }

    fn arrow_stream(&self) -> Result<ArrowStream> {
        Ok(ArrowStream::from_reader(LazyBatchReader::new(self.clone())))
    }
}

/// Pull-based reader that owns its source and materialises at most one batch
/// per call. Creating the reader is `O(1)` in both time and auxiliary memory.
struct LazyBatchReader<S> {
    source: S,
    next: usize,
    end: usize,
}

impl<S: ArrowTableExport> LazyBatchReader<S> {
    fn new(source: S) -> Self {
        let end = source.batch_count();
        Self {
            source,
            next: 0,
            end,
        }
    }
}

impl<S: ArrowTableExport> Iterator for LazyBatchReader<S> {
    type Item = Result<RecordBatch>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.next >= self.end {
            return None;
        }
        let index = self.next;
        self.next += 1;
        let batch = self.source.batch(index);
        if batch.is_err() {
            self.next = self.end;
        }
        Some(batch)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.end.saturating_sub(self.next);
        (remaining, Some(remaining))
    }
}

impl<S: ArrowTableExport> RecordBatchReader for LazyBatchReader<S> {
    fn schema(&self) -> SchemaRef {
        self.source.schema()
    }
}

#[cfg(test)]
#[path = "stream_tests.rs"]
mod tests;
