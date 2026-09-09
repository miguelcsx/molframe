//! One transactional BCIF pull, separated from the public source state machine.

use super::{
    BcifBatchSource, ContinuityLevel, StructureBatchError, address_overflow, append_rows,
    decode_columns, reserve_bytes,
};
use pdbiox_core::{
    Backpressure, Batch, BatchContinuity, BatchDemand, BatchLease, ExecutionContext, LogicalRow,
    SourceBytes, StructureBatch, StructureBatchBuffer, StructureBatchBuilder,
};

struct CandidateBatch {
    batch: StructureBatchBuffer,
    emitted: u32,
    after: ContinuityLevel,
}

struct DecodedCandidate {
    batch: Option<CandidateBatch>,
    remaining: u64,
}

type PullResult =
    Result<(Option<Backpressure<BatchLease<StructureBatch>>>, bool), StructureBatchError>;

impl<S: SourceBytes> BcifBatchSource<S> {
    pub(super) fn pull_once(
        &mut self,
        demand: BatchDemand,
        demanded_capacity: u32,
        context: &ExecutionContext,
    ) -> PullResult {
        let source_rows = u64::from(demanded_capacity).min(self.rows_remaining);
        let source_rows = usize::try_from(source_rows).map_err(|_| address_overflow())?;
        let capacity = u32::try_from(source_rows).map_err(|_| address_overflow())?;
        if !self.ensure_decode_capacity(source_rows)? {
            return Ok((None, true));
        }
        let reserve = reserve_bytes(capacity, self.max_text_bytes_per_row);
        let descriptor = pdbiox_core::ChunkDescriptor::new(
            self.dataset,
            self.chunk,
            self.logical_row,
            capacity,
        )?;
        let Some((builder, reservation)) =
            self.pool
                .acquire(descriptor, reserve, demand.max_bytes, context)?
        else {
            return Ok((None, true));
        };
        let candidate = self.decode_candidate(source_rows, builder)?;
        let Some(candidate_batch) = candidate.batch else {
            self.rows_remaining = candidate.remaining;
            self.finished = candidate.remaining == 0;
            return Ok((None, false));
        };
        if !demand.accepts(candidate_batch.batch.batch()) {
            self.checkpoints.restore(&mut self.columns);
            let required = candidate_batch.batch.batch().retained_bytes();
            self.pool.recycle(candidate_batch.batch, reservation);
            return Err(StructureBatchError::DemandTooSmall {
                required,
                available: demand.max_bytes,
            });
        }
        let Some(next_row) = self
            .logical_row
            .get()
            .checked_add(u64::from(candidate_batch.emitted))
            .map(LogicalRow::new)
        else {
            self.checkpoints.restore(&mut self.columns);
            self.pool.recycle(candidate_batch.batch, reservation);
            return Err(address_overflow());
        };
        let Some(next_chunk) = self.chunk.next() else {
            self.checkpoints.restore(&mut self.columns);
            self.pool.recycle(candidate_batch.batch, reservation);
            return Err(address_overflow());
        };
        let lease = self
            .pool
            .lease(candidate_batch.batch, reservation)
            .inspect_err(|_error| {
                self.checkpoints.restore(&mut self.columns);
            })?;
        self.rows_remaining = candidate.remaining;
        self.logical_row = next_row;
        self.chunk = next_chunk;
        self.continuation = candidate_batch.after;
        self.finished = candidate.remaining == 0;
        Ok((Some(Backpressure::Ready(lease)), false))
    }

    fn decode_candidate(
        &mut self,
        source_rows: usize,
        mut builder: StructureBatchBuilder,
    ) -> Result<DecodedCandidate, StructureBatchError> {
        self.checkpoints.capture(&self.columns);
        let decoded = (|| {
            decode_columns(
                &mut self.source,
                &mut self.columns,
                &mut self.scratch,
                source_rows,
            )?;
            let emitted = append_rows(
                &mut builder,
                &self.columns,
                &self.fields,
                &self.scratch,
                &self.options,
                source_rows,
            )?;
            let consumed = u64::try_from(source_rows).map_err(|_| address_overflow())?;
            let remaining = self
                .rows_remaining
                .checked_sub(consumed)
                .ok_or_else(address_overflow)?;
            if emitted == 0 {
                return Ok(DecodedCandidate {
                    batch: None,
                    remaining,
                });
            }
            let after = if remaining == 0 {
                ContinuityLevel::None
            } else {
                ContinuityLevel::Residue
            };
            builder.set_continuity(BatchContinuity {
                before: self.continuation,
                after,
            });
            Ok(DecodedCandidate {
                batch: Some(CandidateBatch {
                    batch: builder.finish()?,
                    emitted,
                    after,
                }),
                remaining,
            })
        })();
        if decoded.is_err() {
            self.checkpoints.restore(&mut self.columns);
        }
        decoded
    }
}
