//! Pull-based `BinaryCIF` coordinate decoding with memory independent of file size.

mod append;
mod budget;
mod decoded;
mod pull;
mod stream;
mod stream_primitive;

use super::index::BcifOffsetIndex;
use append::append_rows;
use budget::{address_overflow, range_bytes, reserve_bytes};
use decoded::{ColumnChunk, DecodedChunk};
use molframe_core::{
    Backpressure, BatchDemand, BatchLease, BatchSource, ChunkId, Code, ContinuityLevel, DatasetId,
    Diagnostic, Element, ExecutionContext, LogicalRow, MemoryReservation, Presence, ReadOptions,
    SourceBytes, StructureAtomRecord, StructureBatch, StructureBatchError, StructureBatchPool,
};
use num_traits::ToPrimitive;
use stream::{ColumnCheckpoints, ColumnStream, PAYLOAD_BUFFER_BYTES};

const RETAINED_BYTES_PER_ROW: usize = 128;
const DICTIONARY_HEADROOM: usize = 64 * 1024;
const FIELD_COUNT: usize = 16;
const INDEX_WINDOW_BYTES: usize = 64 * 1024;

/// Incremental `BinaryCIF` coordinate rows decoded directly from indexed payloads.
#[derive(Debug)]
pub struct BcifBatchSource<S: SourceBytes> {
    source: S,
    columns: Vec<ColumnStream>,
    fields: [Option<usize>; FIELD_COUNT],
    options: ReadOptions,
    dataset: DatasetId,
    chunk: ChunkId,
    logical_row: LogicalRow,
    rows_remaining: u64,
    max_text_bytes_per_row: usize,
    continuation: ContinuityLevel,
    finished: bool,
    checkpoints: ColumnCheckpoints,
    scratch: Vec<DecodedChunk>,
    _index: BcifOffsetIndex,
    workspace: MemoryReservation,
    pool: StructureBatchPool,
}

impl<S: SourceBytes> BcifBatchSource<S> {
    /// Indexes container metadata without reading large binary payloads.
    ///
    /// # Errors
    ///
    /// Returns a syntax, addressability, codec or memory-budget diagnostic.
    pub fn new(
        mut source: S,
        options: ReadOptions,
        dataset: DatasetId,
        first_chunk: ChunkId,
        first_row: LogicalRow,
        window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, StructureBatchError> {
        let index =
            BcifOffsetIndex::build(&mut source, window_bytes.min(INDEX_WINDOW_BYTES), context)?;
        let Some(atom_site) = index
            .categories
            .iter()
            .find(|category| category.name.trim_start_matches('_') == "atom_site")
        else {
            return Err(Diagnostic::new(Code::E2001)
                .with_message("BinaryCIF block contains no atom_site category")
                .into());
        };
        let rows_remaining = atom_site.rows;
        let selected = atom_site
            .columns
            .iter()
            .filter(|column| field_index(&column.name).is_some());
        let workspace_bytes = selected.clone().try_fold(0usize, |total, column| {
            let data = range_bytes(&column.data.encoding)?;
            let mask = match &column.mask {
                Some(mask) => range_bytes(&mask.encoding)?,
                None => 0,
            };
            total
                .checked_add(data)
                .and_then(|value| value.checked_add(mask))
                .and_then(|value| value.checked_add(PAYLOAD_BUFFER_BYTES * 2))
                .ok_or_else(address_overflow)
        })?;
        let mut workspace = context.try_reserve(workspace_bytes)?;
        let mut columns = Vec::new();
        let mut fields = [None; FIELD_COUNT];
        for indexed in atom_site
            .columns
            .iter()
            .filter_map(|column| field_index(&column.name).map(|field| (field, column.clone())))
        {
            let (field, column) = indexed;
            if fields[field].is_some() {
                return Err(Diagnostic::new(Code::E1401)
                    .with_message("BinaryCIF atom_site contains a duplicate column")
                    .with_context("item", column.name)
                    .into());
            }
            let item = column.name.clone();
            let stream = ColumnStream::new(&mut source, column, window_bytes).map_err(|error| {
                Diagnostic::new(error.code())
                    .with_message(format!("BinaryCIF column {item}: {error}"))
                    .with_context("item", item)
            })?;
            fields[field] = Some(columns.len());
            columns.push(stream);
        }
        let max_text_bytes_per_row = columns
            .iter()
            .map(ColumnStream::max_text_width)
            .sum::<usize>();
        let checkpoints = ColumnCheckpoints::new(&columns);
        let mut scratch = Vec::with_capacity(columns.len());
        for column in &columns {
            scratch.push(column.empty_decoded());
        }
        let retained = columns
            .iter()
            .map(ColumnStream::retained_bytes)
            .sum::<usize>()
            .saturating_add(checkpoints.retained_bytes())
            .saturating_add(
                scratch
                    .capacity()
                    .saturating_mul(std::mem::size_of::<DecodedChunk>()),
            );
        if retained > workspace.bytes() {
            workspace.try_grow(retained - workspace.bytes())?;
        }
        workspace.shrink_to(retained);
        Ok(Self {
            source,
            columns,
            fields,
            options,
            dataset,
            chunk: first_chunk,
            logical_row: first_row,
            rows_remaining,
            max_text_bytes_per_row,
            continuation: ContinuityLevel::None,
            finished: false,
            checkpoints,
            scratch,
            _index: index,
            workspace,
            pool: StructureBatchPool::new(),
        })
    }

    fn ensure_decode_capacity(&mut self, rows: usize) -> Result<bool, Diagnostic> {
        let before = self
            .scratch
            .iter()
            .map(DecodedChunk::capacity_bytes)
            .sum::<usize>();
        let additional = self
            .scratch
            .iter()
            .map(|chunk| chunk.additional_capacity_bytes(rows))
            .sum::<usize>();
        if additional > 0 && self.workspace.try_grow(additional).is_err() {
            return Ok(false);
        }
        for chunk in &mut self.scratch {
            chunk.ensure_capacity(rows)?;
        }
        let after = self
            .scratch
            .iter()
            .map(DecodedChunk::capacity_bytes)
            .sum::<usize>();
        let actual_growth = after.saturating_sub(before);
        if actual_growth > additional {
            self.workspace
                .try_grow(actual_growth - additional)
                .map_err(|error| {
                    Diagnostic::new(Code::E1902)
                        .with_message("BinaryCIF decode workspace exceeds its reservation")
                        .with_context("reason", error.to_string())
                })?;
        }
        Ok(true)
    }
}

impl<S: SourceBytes> BatchSource for BcifBatchSource<S> {
    type Batch = StructureBatch;
    type Error = StructureBatchError;

    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<Self::Batch>>, Self::Error> {
        if !demand.can_accept_work() {
            return Ok(Backpressure::Pending);
        }
        let mut capacity = row_capacity(demand)?;
        let minimum = reserve_bytes(1, self.max_text_bytes_per_row);
        if minimum > context.memory_budget().bytes() {
            return Err(StructureBatchError::RecordExceedsBudget {
                required: minimum,
                available: context.memory_budget().bytes(),
            });
        }
        loop {
            if self.finished || self.rows_remaining == 0 {
                return Ok(Backpressure::Finished);
            }
            if context.cancellation().is_cancelled() {
                return Err(Diagnostic::new(Code::E1901)
                    .with_message("execution cancelled")
                    .into());
            }
            match self.pull_once(demand, capacity, context)? {
                (Some(outcome), _) => return Ok(outcome),
                (None, false) => {}
                (None, true) if capacity > 1 => {
                    capacity = capacity.div_ceil(2);
                }
                (None, true) => return Ok(Backpressure::Pending),
            }
        }
    }
}

fn decode_columns<S: SourceBytes>(
    source: &mut S,
    columns: &mut [ColumnStream],
    chunks: &mut [DecodedChunk],
    rows: usize,
) -> Result<(), Diagnostic> {
    for (column, chunk) in columns.iter_mut().zip(chunks) {
        if let Err(error) = column.decode_into(source, rows, chunk) {
            return Err(Diagnostic::new(error.code())
                .with_message(format!("BinaryCIF column {}: {error}", column.name))
                .with_context("item", column.name.clone()));
        }
    }
    Ok(())
}

struct RowView<'a> {
    columns: &'a [ColumnStream],
    fields: &'a [Option<usize>; FIELD_COUNT],
    chunks: &'a [DecodedChunk],
    row: usize,
}

impl RowView<'_> {
    fn accepted(&self, options: &ReadOptions) -> bool {
        let element = Element::from_symbol(self.text(Field::Element));
        let hydrogen = element.is_some_and(Element::is_hydrogen);
        !(options.discard_hydrogens && hydrogen
            || options.only_first_model && self.integer_i32(Field::Model, 1) != 1)
    }

    fn record(&self) -> StructureAtomRecord<'_> {
        let position = self
            .real(Field::X)
            .zip(self.real(Field::Y))
            .zip(self.real(Field::Z))
            .and_then(|((x, y), z)| {
                x.to_f32()
                    .zip(y.to_f32())
                    .zip(z.to_f32())
                    .map(|((x, y), z)| [x, y, z])
            });
        StructureAtomRecord {
            model: self.integer_i32(Field::Model, 1),
            chain: self.text(Field::Chain),
            component: self.text(Field::Component),
            sequence: self.integer_i32(Field::Sequence, i32::MIN),
            insertion: self.text(Field::Insertion),
            atom: self.text(Field::Atom),
            alternate: self.text(Field::Alternate),
            element: match Element::from_symbol(self.text(Field::Element)) {
                Some(element) => element,
                None => Element::UNKNOWN,
            },
            position,
            occupancy: self.real_presence(Field::Occupancy, 1.0),
            b_factor: self.real_presence(Field::BFactor, 0.0),
            formal_charge: self.integer_presence(Field::Charge),
            atom_site_id: match self
                .integer(Field::Id)
                .and_then(|value| u32::try_from(value).ok())
            {
                Some(value) => value,
                None => 0,
            },
            heterogen: self.text(Field::Group).eq_ignore_ascii_case("HETATM"),
        }
    }

    fn text(&self, field: Field) -> &str {
        let Some((column, chunk)) = self.column(field) else {
            return "";
        };
        if presence(chunk, self.row) != Presence::Present {
            return "";
        }
        let ColumnChunk::Text(values) = &chunk.values else {
            return "";
        };
        let Some(index) = values.get(self.row) else {
            return "";
        };
        column.text(*index).map_or("", |text| text)
    }

    fn integer(&self, field: Field) -> Option<i64> {
        let (_, chunk) = self.column(field)?;
        if presence(chunk, self.row) != Presence::Present {
            return None;
        }
        match &chunk.values {
            ColumnChunk::Integer(values) => values.get(self.row).copied(),
            ColumnChunk::Float(_) | ColumnChunk::Text(_) => None,
        }
    }

    fn real(&self, field: Field) -> Option<f64> {
        let (_, chunk) = self.column(field)?;
        if presence(chunk, self.row) != Presence::Present {
            return None;
        }
        match &chunk.values {
            ColumnChunk::Integer(values) => values.get(self.row).and_then(ToPrimitive::to_f64),
            ColumnChunk::Float(values) => values.get(self.row).map(|value| f64::from(*value)),
            ColumnChunk::Text(_) => None,
        }
    }

    fn integer_i32(&self, field: Field, absent: i32) -> i32 {
        match self
            .integer(field)
            .and_then(|value| i32::try_from(value).ok())
        {
            Some(value) => value,
            None => absent,
        }
    }

    fn real_presence(&self, field: Field, absent: f32) -> (f32, Presence) {
        let Some((_, chunk)) = self.column(field) else {
            return (absent, Presence::Unknown);
        };
        let validity = presence(chunk, self.row);
        let value = match self.real(field).and_then(|value| value.to_f32()) {
            Some(value) => value,
            None => absent,
        };
        (value, validity)
    }

    fn integer_presence(&self, field: Field) -> (i8, Presence) {
        let Some((_, chunk)) = self.column(field) else {
            return (0, Presence::Unknown);
        };
        let validity = presence(chunk, self.row);
        let value = match self
            .integer(field)
            .and_then(|value| i8::try_from(value).ok())
        {
            Some(value) => value,
            None => 0,
        };
        (value, validity)
    }

    fn column(&self, field: Field) -> Option<(&ColumnStream, &DecodedChunk)> {
        let index = (*self.fields.get(field as usize)?)?;
        Some((self.columns.get(index)?, self.chunks.get(index)?))
    }
}

fn presence(chunk: &DecodedChunk, row: usize) -> Presence {
    match chunk.mask.as_ref().and_then(|mask| mask.get(row)).copied() {
        Some(1) => Presence::Inapplicable,
        Some(0) | None => Presence::Present,
        Some(2..) => Presence::Unknown,
    }
}

#[derive(Clone, Copy)]
#[repr(usize)]
enum Field {
    Group = 0,
    Id,
    Element,
    Atom,
    Alternate,
    Component,
    Chain,
    Sequence,
    Insertion,
    X,
    Y,
    Z,
    Occupancy,
    BFactor,
    Charge,
    Model,
}

fn field_index(name: &str) -> Option<usize> {
    let field = match name.as_bytes() {
        b"group_PDB" => Field::Group,
        b"id" => Field::Id,
        b"type_symbol" => Field::Element,
        b"label_atom_id" => Field::Atom,
        b"label_alt_id" => Field::Alternate,
        b"label_comp_id" => Field::Component,
        b"label_asym_id" => Field::Chain,
        b"label_seq_id" => Field::Sequence,
        b"pdbx_PDB_ins_code" => Field::Insertion,
        b"Cartn_x" => Field::X,
        b"Cartn_y" => Field::Y,
        b"Cartn_z" => Field::Z,
        b"occupancy" => Field::Occupancy,
        b"B_iso_or_equiv" => Field::BFactor,
        b"pdbx_formal_charge" => Field::Charge,
        b"pdbx_PDB_model_num" => Field::Model,
        _ => return None,
    };
    Some(field as usize)
}

fn row_capacity(demand: BatchDemand) -> Result<u32, StructureBatchError> {
    let retained = match demand
        .max_bytes
        .saturating_sub(DICTIONARY_HEADROOM)
        .checked_div(RETAINED_BYTES_PER_ROW)
    {
        Some(rows) => rows,
        None => 0,
    };
    let rows = retained.min(demand.max_rows).min(u32::MAX as usize);
    if rows == 0 {
        return Err(StructureBatchError::DemandTooSmall {
            required: DICTIONARY_HEADROOM + RETAINED_BYTES_PER_ROW,
            available: demand.max_bytes,
        });
    }
    u32::try_from(rows).map_err(|_| address_overflow())
}

#[cfg(test)]
#[path = "batch_tests.rs"]
mod tests;
