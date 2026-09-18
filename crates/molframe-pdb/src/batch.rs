//! Pull-based fixed-column reading over bounded source windows.

use crate::{fixed, hybrid36};
use molframe_core::{
    Backpressure, Batch, BatchContinuity, BatchDemand, BatchLease, BatchSource, ChunkId, Code,
    ContinuityLevel, DatasetId, Diagnostic, Element, ExecutionContext, LogicalRow,
    MissingElementPolicy, Presence, ReadOptions, SourceBytes, StructureAtomRecord, StructureBatch,
    StructureBatchBuffer, StructureBatchBuilder, StructureBatchError, StructureBatchPool,
};
use num_traits::ToPrimitive;

const RETAINED_BYTES_PER_ROW: usize = 128;
const DICTIONARY_HEADROOM: usize = 64 * 1024;

/// Incremental PDB rows over a seekable or resident byte source.
#[derive(Debug)]
pub struct PdbBatchSource<S> {
    source: S,
    options: ReadOptions,
    dataset: DatasetId,
    next_chunk: ChunkId,
    next_row: LogicalRow,
    cursor: u64,
    window_bytes: usize,
    model: i32,
    continuation: ContinuityLevel,
    finished: bool,
    pool: StructureBatchPool,
}

impl<S: SourceBytes> PdbBatchSource<S> {
    /// Creates a bounded reader. `window_bytes` is an upper bound and never
    /// grows in response to source length.
    ///
    /// # Errors
    ///
    /// Returns [`StructureBatchError::DemandTooSmall`] for a zero-sized window.
    pub fn new(
        source: S,
        options: ReadOptions,
        dataset: DatasetId,
        first_chunk: ChunkId,
        first_row: LogicalRow,
        window_bytes: usize,
    ) -> Result<Self, StructureBatchError> {
        if window_bytes == 0 {
            return Err(StructureBatchError::DemandTooSmall {
                required: 1,
                available: 0,
            });
        }
        Ok(Self {
            source,
            options,
            dataset,
            next_chunk: first_chunk,
            next_row: first_row,
            cursor: 0,
            window_bytes,
            model: 1,
            continuation: ContinuityLevel::None,
            finished: false,
            pool: StructureBatchPool::new(),
        })
    }

    fn row_capacity(demand: BatchDemand) -> Result<u32, StructureBatchError> {
        let byte_rows = match demand
            .max_bytes
            .saturating_sub(DICTIONARY_HEADROOM)
            .checked_div(RETAINED_BYTES_PER_ROW)
        {
            Some(rows) => rows,
            None => 0,
        };
        let rows = demand.max_rows.min(byte_rows).min(u32::MAX as usize);
        u32::try_from(rows).map_err(|_| StructureBatchError::DemandTooSmall {
            required: RETAINED_BYTES_PER_ROW,
            available: demand.max_bytes,
        })
    }

    fn reserve_bytes(capacity: u32, demand: BatchDemand) -> usize {
        (capacity as usize)
            .saturating_mul(RETAINED_BYTES_PER_ROW)
            .saturating_add(DICTIONARY_HEADROOM)
            .min(demand.max_bytes)
    }

    fn pull_batches(
        &mut self,
        demand: BatchDemand,
        capacity: u32,
        reserved: usize,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<StructureBatch>>, StructureBatchError> {
        loop {
            let descriptor = molframe_core::ChunkDescriptor::new(
                self.dataset,
                self.next_chunk,
                self.next_row,
                capacity,
            )?;
            let Some((builder, reservation)) =
                self.pool
                    .acquire(descriptor, reserved, demand.max_bytes, context)?
            else {
                return Ok(Backpressure::Pending);
            };
            let requested = self.window_bytes.min(demand.max_bytes.max(1));
            let length = self.source.len_hint();
            let window = self.source.window(self.cursor, requested)?;
            let bytes = window.bytes();
            if bytes.is_empty() {
                self.finished = true;
                return Ok(Backpressure::Finished);
            }
            let source_finished = length
                .is_some_and(|length| window.end().is_ok_and(|window_end| window_end >= length));
            let complete = complete_prefix(bytes, source_finished);
            if complete == 0 {
                let required = bytes.len().saturating_add(1);
                return Err(StructureBatchError::RecordExceedsBudget {
                    required,
                    available: requested,
                });
            }
            let model_checkpoint = self.model;
            let parsed = match parse_window(
                &bytes[..complete],
                builder,
                capacity,
                self.continuation,
                &self.options,
                &mut self.model,
            ) {
                Ok(parsed) => parsed,
                Err(error) => {
                    self.model = model_checkpoint;
                    return Err(error);
                }
            };
            if parsed.rows == 0 {
                self.pool.recycle(parsed.batch, reservation);
                self.cursor = advance(self.cursor, parsed.consumed)?;
                if source_finished && parsed.consumed == complete {
                    self.finished = true;
                    return Ok(Backpressure::Finished);
                }
                continue;
            }
            let batch = parsed.batch;
            if !demand.accepts(batch.batch()) {
                self.model = model_checkpoint;
                let required = batch.batch().retained_bytes();
                self.pool.recycle(batch, reservation);
                return Err(StructureBatchError::DemandTooSmall {
                    required,
                    available: demand.max_bytes,
                });
            }
            let next_cursor = advance(self.cursor, parsed.consumed)?;
            let (next_row, next_chunk) =
                match next_identity(self.next_row, self.next_chunk, parsed.rows) {
                    Ok(identity) => identity,
                    Err(error) => {
                        self.model = model_checkpoint;
                        self.pool.recycle(batch, reservation);
                        return Err(error);
                    }
                };
            let continuation = if parsed.ended_for_demand {
                ContinuityLevel::Residue
            } else {
                ContinuityLevel::None
            };
            let lease = match self.pool.lease(batch, reservation) {
                Ok(lease) => lease,
                Err(error) => {
                    self.model = model_checkpoint;
                    return Err(error);
                }
            };
            self.cursor = next_cursor;
            self.next_row = next_row;
            self.next_chunk = next_chunk;
            self.continuation = continuation;
            return Ok(Backpressure::Ready(lease));
        }
    }
}

impl<S: SourceBytes> BatchSource for PdbBatchSource<S> {
    type Batch = StructureBatch;
    type Error = StructureBatchError;

    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<Self::Batch>>, Self::Error> {
        if self.finished {
            return Ok(Backpressure::Finished);
        }
        if !demand.can_accept_work() {
            return Ok(Backpressure::Pending);
        }
        let capacity = Self::row_capacity(demand)?;
        if capacity == 0 {
            return Err(StructureBatchError::DemandTooSmall {
                required: RETAINED_BYTES_PER_ROW + DICTIONARY_HEADROOM,
                available: demand.max_bytes,
            });
        }
        self.pull_batches(
            demand,
            capacity,
            Self::reserve_bytes(capacity, demand),
            context,
        )
    }
}

fn next_identity(
    row: LogicalRow,
    chunk: ChunkId,
    rows: u32,
) -> Result<(LogicalRow, ChunkId), StructureBatchError> {
    let row = row
        .get()
        .checked_add(u64::from(rows))
        .map(LogicalRow::new)
        .ok_or_else(identity_overflow)?;
    let chunk = chunk.next().ok_or_else(identity_overflow)?;
    Ok((row, chunk))
}

struct ParsedWindow {
    batch: StructureBatchBuffer,
    consumed: usize,
    rows: u32,
    ended_for_demand: bool,
}

fn parse_window(
    bytes: &[u8],
    mut builder: StructureBatchBuilder,
    capacity: u32,
    continuity: ContinuityLevel,
    options: &ReadOptions,
    model: &mut i32,
) -> Result<ParsedWindow, StructureBatchError> {
    let mut consumed = 0usize;
    let mut rows = 0u32;
    for raw in bytes.split_inclusive(|byte| *byte == b'\n') {
        let line_bytes = match raw.strip_suffix(b"\n") {
            Some(value) => value,
            None => raw,
        };
        let line_bytes = match line_bytes.strip_suffix(b"\r") {
            Some(value) => value,
            None => line_bytes,
        };
        let line = std::str::from_utf8(line_bytes)
            .map_err(|_| Diagnostic::new(Code::E1201).with_message("input is not valid text"))?;
        let record = fixed::record(line);
        if record == "MODEL" {
            *model = match fixed::integer(line, 11, 14).and_then(|value| i32::try_from(value).ok())
            {
                Some(value) => value,
                None => 1,
            };
        } else if matches!(record, "ATOM" | "HETATM") {
            if rows == capacity {
                break;
            }
            builder.push(atom_record(line, *model, options))?;
            rows += 1;
        }
        consumed = consumed.saturating_add(raw.len());
    }
    let ended_for_demand = rows == capacity && consumed < bytes.len();
    builder.set_continuity(BatchContinuity {
        before: continuity,
        after: if ended_for_demand {
            ContinuityLevel::Residue
        } else {
            ContinuityLevel::None
        },
    });
    Ok(ParsedWindow {
        batch: builder.finish()?,
        consumed,
        rows,
        ended_for_demand,
    })
}

fn atom_record<'a>(line: &'a str, model: i32, options: &ReadOptions) -> StructureAtomRecord<'a> {
    let raw_atom = fixed::raw(line, 13, 16);
    let declared = fixed::text(line, 77, 78);
    let element = match Element::from_symbol(declared) {
        Some(element) => element,
        None => match options.missing_element_policy {
            MissingElementPolicy::PreserveUnknown => Element::UNKNOWN,
            MissingElementPolicy::InferFromAtomName => Element::infer_from_pdb_atom_name(raw_atom),
        },
    };
    let position = fixed::real(line, 31, 38)
        .and_then(|value| value.to_f32())
        .zip(fixed::real(line, 39, 46).and_then(|value| value.to_f32()))
        .zip(fixed::real(line, 47, 54).and_then(|value| value.to_f32()))
        .map(|((x, y), z)| [x, y, z]);
    StructureAtomRecord {
        model,
        chain: fixed::text(line, 22, 22),
        component: fixed::text(line, 18, 20),
        sequence: match hybrid36::decode(fixed::raw(line, 23, 26), 4)
            .and_then(|value| i32::try_from(value).ok())
        {
            Some(sequence) => sequence,
            None => i32::MIN,
        },
        insertion: fixed::text(line, 27, 27),
        atom: raw_atom.trim(),
        alternate: fixed::text(line, 17, 17),
        element,
        position,
        occupancy: optional_real(line, 55, 60, 1.0),
        b_factor: optional_real(line, 61, 66, 0.0),
        formal_charge: formal_charge(fixed::text(line, 79, 80)),
        atom_site_id: match hybrid36::decode(fixed::raw(line, 7, 11), 5)
            .and_then(|value| u32::try_from(value).ok())
        {
            Some(identifier) => identifier,
            None => 0,
        },
        heterogen: fixed::record(line) == "HETATM",
    }
}

fn optional_real(line: &str, start: usize, end: usize, absent: f32) -> (f32, Presence) {
    match fixed::real(line, start, end).and_then(|value| value.to_f32()) {
        Some(value) => (value, Presence::Present),
        None => (absent, Presence::Unknown),
    }
}

fn formal_charge(text: &str) -> (i8, Presence) {
    let bytes = text.as_bytes();
    if bytes.len() != 2 || !bytes[0].is_ascii_digit() {
        return (0, Presence::Inapplicable);
    }
    let magnitude = (bytes[0] - b'0').cast_signed();
    match bytes[1] {
        b'+' => (magnitude, Presence::Present),
        b'-' => (-magnitude, Presence::Present),
        _ => (0, Presence::Unknown),
    }
}

fn complete_prefix(bytes: &[u8], source_finished: bool) -> usize {
    if source_finished {
        return bytes.len();
    }
    bytes
        .iter()
        .rposition(|byte| *byte == b'\n')
        .map_or(0, |position| position + 1)
}

fn advance(start: u64, count: usize) -> Result<u64, StructureBatchError> {
    let count = u64::try_from(count).map_err(|_| identity_overflow())?;
    start.checked_add(count).ok_or_else(identity_overflow)
}

fn identity_overflow() -> StructureBatchError {
    StructureBatchError::Diagnostic(
        Diagnostic::new(Code::E1903).with_message("PDB batch identity exceeds u64"),
    )
}

#[cfg(test)]
#[path = "batch_tests.rs"]
mod tests;
