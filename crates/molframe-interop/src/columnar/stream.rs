//! Owned Arrow C Stream values.

use arrow::datatypes::SchemaRef;
use arrow::error::{ArrowError, Result};
use arrow::ffi_stream::FFI_ArrowArrayStream;
use arrow::record_batch::{RecordBatch, RecordBatchReader};
use std::ops::Range;

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
    /// Rows in the whole table.
    fn row_count(&self) -> usize;
    /// The rows batch `index` holds, known without materialising it.
    fn batch_rows(&self, index: usize) -> Range<usize>;

    fn record_batches(&self) -> Result<Vec<RecordBatch>> {
        (0..self.batch_count())
            .map(|index| self.batch(index))
            .collect()
    }

    fn arrow_stream(&self) -> Result<ArrowStream> {
        Ok(ArrowStream::from_reader(LazyBatchReader::new(self.clone())))
    }

    /// A lazy stream over only the rows in `rows`.
    ///
    /// Batches wholly outside the window are never materialised; the first and
    /// last overlapping batches are sliced, which shares their buffers.
    fn arrow_stream_rows(&self, rows: Range<usize>) -> Result<ArrowStream> {
        if rows.start > rows.end || rows.end > self.row_count() {
            return Err(ArrowError::InvalidArgumentError(
                "row window lies outside the table".to_owned(),
            ));
        }
        Ok(ArrowStream::from_reader(LazyBatchReader::new(
            RowWindow::new(self.clone(), rows),
        )))
    }
}

/// A table restricted to a contiguous window of rows.
#[derive(Clone)]
struct RowWindow<S> {
    source: S,
    rows: Range<usize>,
    first_batch: usize,
    batches: usize,
}

impl<S: ArrowTableExport> RowWindow<S> {
    fn new(source: S, rows: Range<usize>) -> Self {
        let total = source.batch_count();
        let mut overlapping = (0..total).filter(|index| {
            let held = source.batch_rows(*index);
            held.start < rows.end && held.end > rows.start
        });
        let first_batch = match overlapping.next() {
            Some(index) => index,
            None => 0,
        };
        let batches = usize::from(rows.start < rows.end) * (1 + overlapping.count());
        Self {
            source,
            rows,
            first_batch,
            batches,
        }
    }
}

impl<S: ArrowTableExport> ArrowTableExport for RowWindow<S> {
    fn schema(&self) -> SchemaRef {
        self.source.schema()
    }

    fn batch_count(&self) -> usize {
        self.batches
    }

    fn batch(&self, index: usize) -> Result<RecordBatch> {
        let absolute = self.first_batch + index;
        let held = self.source.batch_rows(absolute);
        let start = held.start.max(self.rows.start);
        let end = held.end.min(self.rows.end);
        let batch = self.source.batch(absolute)?;
        Ok(batch.slice(start - held.start, end - start))
    }

    fn row_count(&self) -> usize {
        self.rows.len()
    }

    fn batch_rows(&self, index: usize) -> Range<usize> {
        let held = self.source.batch_rows(self.first_batch + index);
        let start = held.start.max(self.rows.start) - self.rows.start;
        let end = held.end.min(self.rows.end) - self.rows.start;
        start..end
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
