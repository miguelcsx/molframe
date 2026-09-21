//! Pull-based coordinate projection over bounded CIF token windows.

mod budget;
mod token;

use budget::{identity_overflow, reserve_bytes, row_capacity};
use molframe_core::{
    Backpressure, Batch, BatchContinuity, BatchDemand, BatchLease, BatchSource, ChunkId,
    ContinuityLevel, DatasetId, Element, ExecutionContext, LogicalRow, Presence, ReadOptions,
    SourceBytes, StructureAtomRecord, StructureBatch, StructureBatchBuffer, StructureBatchBuilder,
    StructureBatchError, StructureBatchPool,
};
use num_traits::ToPrimitive;
use token::TokenCursor;

/// Incremental mmCIF coordinate rows.
#[derive(Debug)]
pub struct MmcifBatchSource<S: SourceBytes> {
    tokens: TokenCursor<S>,
    options: ReadOptions,
    dataset: DatasetId,
    chunk: ChunkId,
    logical_row: LogicalRow,
    parser: ParserState,
    continuation: ContinuityLevel,
    pool: StructureBatchPool,
}

impl<S: SourceBytes> MmcifBatchSource<S> {
    /// Creates a reader whose token and source windows never grow with file size.
    ///
    /// # Errors
    ///
    /// Returns a budget error when the bounded token carry cannot be reserved.
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
            tokens: TokenCursor::new(source, window_bytes, context)?,
            options,
            dataset,
            chunk: first_chunk,
            logical_row: first_row,
            parser: ParserState::new(),
            continuation: ContinuityLevel::None,
            pool: StructureBatchPool::new(),
        })
    }
}

impl<S: SourceBytes> BatchSource for MmcifBatchSource<S> {
    type Batch = StructureBatch;
    type Error = StructureBatchError;

    fn next_batch(
        &mut self,
        demand: BatchDemand,
        context: &ExecutionContext,
    ) -> Result<Backpressure<BatchLease<Self::Batch>>, Self::Error> {
        if self.parser.finished {
            return Ok(Backpressure::Finished);
        }
        if !demand.can_accept_work() {
            return Ok(Backpressure::Pending);
        }
        let capacity = row_capacity(demand)?;
        let reserve = reserve_bytes(capacity, demand);
        let descriptor = molframe_core::ChunkDescriptor::new(
            self.dataset,
            self.chunk,
            self.logical_row,
            capacity,
        )?;
        let Some((builder, reservation)) =
            self.pool
                .acquire(descriptor, reserve, demand.max_bytes, context)?
        else {
            return Ok(Backpressure::Pending);
        };
        let token_checkpoint = self.tokens.checkpoint();
        let parser_checkpoint = self.parser.checkpoint();
        let parsed = match parse_batch(
            &mut self.tokens,
            &mut self.parser,
            &self.options,
            builder,
            capacity,
            self.continuation,
        ) {
            Ok(parsed) => parsed,
            Err(error) => {
                self.tokens.restore(token_checkpoint);
                self.parser.restore(parser_checkpoint);
                return Err(error);
            }
        };
        let Some(batch) = parsed.batch else {
            return Ok(Backpressure::Finished);
        };
        if !demand.accepts(batch.batch()) {
            self.tokens.restore(token_checkpoint);
            self.parser.restore(parser_checkpoint);
            let required = batch.batch().retained_bytes();
            self.pool.recycle(batch, reservation);
            return Err(StructureBatchError::DemandTooSmall {
                required,
                available: demand.max_bytes,
            });
        }
        let Some(next_row) = self
            .logical_row
            .get()
            .checked_add(u64::from(parsed.rows))
            .map(LogicalRow::new)
        else {
            self.tokens.restore(token_checkpoint);
            self.parser.restore(parser_checkpoint);
            self.pool.recycle(batch, reservation);
            return Err(identity_overflow());
        };
        let Some(next_chunk) = self.chunk.next() else {
            self.tokens.restore(token_checkpoint);
            self.parser.restore(parser_checkpoint);
            self.pool.recycle(batch, reservation);
            return Err(identity_overflow());
        };
        let lease = match self.pool.lease(batch, reservation) {
            Ok(lease) => lease,
            Err(error) => {
                self.tokens.restore(token_checkpoint);
                self.parser.restore(parser_checkpoint);
                return Err(error);
            }
        };
        self.logical_row = next_row;
        self.chunk = next_chunk;
        self.continuation = parsed.after;
        Ok(Backpressure::Ready(lease))
    }
}

struct ParsedBatch {
    batch: Option<StructureBatchBuffer>,
    rows: u32,
    after: ContinuityLevel,
}

fn parse_batch<S: SourceBytes>(
    tokens: &mut TokenCursor<S>,
    parser: &mut ParserState,
    options: &ReadOptions,
    mut builder: StructureBatchBuilder,
    capacity: u32,
    before: ContinuityLevel,
) -> Result<ParsedBatch, StructureBatchError> {
    let mut rows = 0u32;
    while rows < capacity && !parser.finished {
        let Some(token) = tokens.next_token()? else {
            parser.finished = true;
            break;
        };
        if parser.token(token, options) && parser.row.accepted(options) {
            builder.push(parser.row.record(options))?;
            rows += 1;
        }
    }
    if rows == 0 {
        return Ok(ParsedBatch {
            batch: None,
            rows,
            after: ContinuityLevel::None,
        });
    }
    let after = if parser.finished {
        ContinuityLevel::None
    } else {
        ContinuityLevel::Residue
    };
    builder.set_continuity(BatchContinuity { before, after });
    Ok(ParsedBatch {
        batch: Some(builder.finish()?),
        rows,
        after,
    })
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    Search,
    Tags,
    Rows,
}

#[derive(Debug)]
struct ParserState {
    mode: Mode,
    columns: Vec<Option<Field>>,
    column: usize,
    row: AtomRow,
    finished: bool,
}

#[derive(Clone, Copy, Debug)]
struct ParserCheckpoint {
    mode: Mode,
    column: usize,
    finished: bool,
}

impl ParserState {
    fn new() -> Self {
        Self {
            mode: Mode::Search,
            columns: Vec::new(),
            column: 0,
            row: AtomRow::new(),
            finished: false,
        }
    }

    fn checkpoint(&self) -> ParserCheckpoint {
        ParserCheckpoint {
            mode: self.mode,
            column: self.column,
            finished: self.finished,
        }
    }

    fn restore(&mut self, checkpoint: ParserCheckpoint) {
        self.mode = checkpoint.mode;
        self.column = checkpoint.column;
        self.finished = checkpoint.finished;
    }

    fn token(&mut self, token: &str, options: &ReadOptions) -> bool {
        match self.mode {
            Mode::Search => {
                if token.eq_ignore_ascii_case("loop_") {
                    self.mode = Mode::Tags;
                    self.columns.clear();
                }
                false
            }
            Mode::Tags => self.tag_or_value(token, options),
            Mode::Rows => self.value(token, options),
        }
    }

    fn tag_or_value(&mut self, token: &str, options: &ReadOptions) -> bool {
        if token.starts_with('_') {
            self.columns.push(Field::from_tag(token));
            return false;
        }
        if !self.columns.iter().any(Option::is_some) {
            self.mode = Mode::Search;
            return false;
        }
        self.mode = Mode::Rows;
        self.value(token, options)
    }

    fn value(&mut self, token: &str, _options: &ReadOptions) -> bool {
        if self.column == 0 {
            if is_control(token) {
                self.finished = true;
                return false;
            }
            self.row.clear();
        }
        let field = self.columns.get(self.column).copied().flatten();
        if let Some(field) = field {
            self.row.set(field, token);
        }
        self.column += 1;
        if self.column < self.columns.len() {
            return false;
        }
        self.column = 0;
        true
    }
}

#[derive(Clone, Copy, Debug)]
enum Field {
    Group,
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

impl Field {
    fn from_tag(tag: &str) -> Option<Self> {
        let item = tag.strip_prefix("_atom_site.")?;
        match item {
            "group_PDB" => Some(Self::Group),
            "id" => Some(Self::Id),
            "type_symbol" => Some(Self::Element),
            "label_atom_id" => Some(Self::Atom),
            "label_alt_id" => Some(Self::Alternate),
            "label_comp_id" => Some(Self::Component),
            "label_asym_id" => Some(Self::Chain),
            "label_seq_id" => Some(Self::Sequence),
            "pdbx_PDB_ins_code" => Some(Self::Insertion),
            "Cartn_x" => Some(Self::X),
            "Cartn_y" => Some(Self::Y),
            "Cartn_z" => Some(Self::Z),
            "occupancy" => Some(Self::Occupancy),
            "B_iso_or_equiv" => Some(Self::BFactor),
            "pdbx_formal_charge" => Some(Self::Charge),
            "pdbx_PDB_model_num" => Some(Self::Model),
            _ => None,
        }
    }
}

#[derive(Debug)]
struct AtomRow {
    text: [String; 7],
    integer: [i64; 4],
    real: [f64; 5],
    present: u32,
}

impl AtomRow {
    fn new() -> Self {
        Self {
            text: core::array::from_fn(|_| String::new()),
            integer: [0; 4],
            real: [0.0; 5],
            present: 0,
        }
    }

    fn clear(&mut self) {
        self.present = 0;
    }

    fn set(&mut self, field: Field, token: &str) {
        match field {
            Field::Group => self.set_text(0, token),
            Field::Element => self.set_text(1, token),
            Field::Atom => self.set_text(2, token),
            Field::Alternate => self.set_text(3, token),
            Field::Component => self.set_text(4, token),
            Field::Chain => self.set_text(5, token),
            Field::Insertion => self.set_text(6, token),
            Field::Id => self.set_integer(0, token),
            Field::Sequence => self.set_integer(1, token),
            Field::Charge => self.set_integer(2, token),
            Field::Model => self.set_integer(3, token),
            Field::X => self.set_real(0, token),
            Field::Y => self.set_real(1, token),
            Field::Z => self.set_real(2, token),
            Field::Occupancy => self.set_real(3, token),
            Field::BFactor => self.set_real(4, token),
        }
    }

    fn set_text(&mut self, index: usize, token: &str) {
        self.text[index].clear();
        if !matches!(token, "." | "?") {
            self.text[index].push_str(token);
            self.present |= 1 << index;
        }
    }

    fn set_integer(&mut self, index: usize, token: &str) {
        if let Ok(value) = token.parse() {
            self.integer[index] = value;
            self.present |= 1 << (7 + index);
        }
    }

    fn set_real(&mut self, index: usize, token: &str) {
        if let Ok(value) = token.parse() {
            self.real[index] = value;
            self.present |= 1 << (11 + index);
        }
    }

    fn accepted(&self, options: &ReadOptions) -> bool {
        let element = Element::from_symbol(&self.text[1]);
        let hydrogen = element.is_some_and(Element::is_hydrogen);
        let model = integer_i32(&self.integer, self.present, 3, 1);
        !(options.discard_hydrogens && hydrogen || options.only_first_model && model != 1)
    }

    fn record(&self, _options: &ReadOptions) -> StructureAtomRecord<'_> {
        let position = if self.present & (0b111 << 11) == 0b111 << 11 {
            self.real[0]
                .to_f32()
                .zip(self.real[1].to_f32())
                .zip(self.real[2].to_f32())
                .map(|((x, y), z)| [x, y, z])
        } else {
            None
        };
        StructureAtomRecord {
            model: integer_i32(&self.integer, self.present, 3, 1),
            chain: &self.text[5],
            component: &self.text[4],
            sequence: integer_i32(&self.integer, self.present, 1, i32::MIN),
            insertion: &self.text[6],
            atom: &self.text[2],
            alternate: &self.text[3],
            element: match Element::from_symbol(&self.text[1]) {
                Some(element) => element,
                None => Element::UNKNOWN,
            },
            position,
            occupancy: real_presence(&self.real, self.present, 3, 1.0),
            b_factor: real_presence(&self.real, self.present, 4, 0.0),
            formal_charge: integer_presence(&self.integer, self.present, 2),
            atom_site_id: integer_u32(&self.integer, self.present, 0),
            heterogen: self.text[0].eq_ignore_ascii_case("HETATM"),
        }
    }
}

fn integer_i32(values: &[i64; 4], present: u32, index: usize, absent: i32) -> i32 {
    if present & (1 << (7 + index)) == 0 {
        return absent;
    }
    match i32::try_from(values[index]) {
        Ok(value) => value,
        Err(_) => absent,
    }
}

fn integer_u32(values: &[i64; 4], present: u32, index: usize) -> u32 {
    if present & (1 << (7 + index)) == 0 {
        return 0;
    }
    let Ok(value) = u32::try_from(values[index]) else {
        return 0;
    };
    value
}

fn integer_presence(values: &[i64; 4], present: u32, index: usize) -> (i8, Presence) {
    if present & (1 << (7 + index)) == 0 {
        return (0, Presence::Unknown);
    }
    match i8::try_from(values[index]) {
        Ok(value) => (value, Presence::Present),
        Err(_) => (0, Presence::Unknown),
    }
}

fn real_presence(values: &[f64; 5], present: u32, index: usize, absent: f32) -> (f32, Presence) {
    if present & (1 << (11 + index)) == 0 {
        return (absent, Presence::Unknown);
    }
    match values[index].to_f32() {
        Some(value) => (value, Presence::Present),
        None => (absent, Presence::Unknown),
    }
}

fn is_control(token: &str) -> bool {
    token.starts_with('_')
        || token.eq_ignore_ascii_case("loop_")
        || token.starts_with("data_")
        || token.starts_with("save_")
        || token.eq_ignore_ascii_case("stop_")
}

#[cfg(test)]
#[path = "batch_tests.rs"]
mod tests;
