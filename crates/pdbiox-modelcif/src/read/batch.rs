//! `ModelCIF` structure batches over the shared incremental CIF tokenizer.

use pdbiox_cif::MmcifBatchSource;
use pdbiox_core::{
    Backpressure, BatchDemand, BatchLease, BatchSource, ChunkId, DatasetId, ExecutionContext,
    LogicalRow, ReadOptions, SourceBytes, StructureBatch, StructureBatchError,
};

/// Incremental `ModelCIF` coordinate rows.
///
/// Structure categories use the exact mmCIF tokenizer and lowering path. Model
/// quality categories remain separately accessible through the compact
/// `ModelCIF` projection rather than being retained by every atom batch.
#[derive(Debug)]
pub struct ModelCifBatchSource<S: SourceBytes> {
    inner: MmcifBatchSource<S>,
}

impl<S: SourceBytes> ModelCifBatchSource<S> {
    /// Creates a bounded `ModelCIF` structure source.
    ///
    /// # Errors
    ///
    /// Returns a typed budget or source error from the shared CIF reader.
    pub fn new(
        source: S,
        options: ReadOptions,
        dataset: DatasetId,
        first_chunk: ChunkId,
        first_row: LogicalRow,
        window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, StructureBatchError> {
        Ok(Self {
            inner: MmcifBatchSource::new(
                source,
                options,
                dataset,
                first_chunk,
                first_row,
                window_bytes,
                context,
            )?,
        })
    }
}

impl<S: SourceBytes> BatchSource for ModelCifBatchSource<S> {
    type Batch = StructureBatch;
    type Error = StructureBatchError;

    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<Self::Batch>>, Self::Error> {
        self.inner.next_batch(demand, context)
    }
}
