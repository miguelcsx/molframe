//! Virtual tera-scale structural inputs consumed by the real batch readers.

use super::{ResourceRecord, measure_case};
use molframe::ReadOptions;
use molframe::core::{
    Backpressure, BatchDemand, BatchSource, ByteWindow, ChunkId, DatasetId, ExecutionContext,
    LogicalRow, MemoryReservation, ScratchPolicy, SourceBytes, StructureBatch, StructureBatchError,
};
use std::hint::black_box;

const WINDOW_BYTES: usize = 4 * 1024 * 1024;
const BATCH_BYTES: usize = 16 * 1024 * 1024;
const CIF_HEADER: &[u8] = b"data_tera\nloop_\n\
_atom_site.group_PDB\n\
_atom_site.id\n\
_atom_site.type_symbol\n\
_atom_site.label_atom_id\n\
_atom_site.label_alt_id\n\
_atom_site.label_comp_id\n\
_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n\
_atom_site.pdbx_PDB_ins_code\n\
_atom_site.Cartn_x\n\
_atom_site.Cartn_y\n\
_atom_site.Cartn_z\n\
_atom_site.occupancy\n\
_atom_site.B_iso_or_equiv\n\
_atom_site.pdbx_formal_charge\n\
_atom_site.pdbx_PDB_model_num\n";
const CIF_ROW: &[u8] = b"ATOM 1 C CA . GLY A 1 ? 1.0 2.0 3.0 1.0 10.0 0 1\n";
const PDB_ROW: &[u8] =
    b"ATOM      1  CA  GLY A   1       1.000   2.000   3.000  1.00 10.00           C  \n";

#[derive(Clone, Copy)]
pub(super) enum GeneratedFormat {
    Mmcif,
    ModelCif,
    Pdb,
}

impl GeneratedFormat {
    pub(super) fn parse(value: &str) -> Option<Self> {
        match value {
            "generated_mmcif_batches" => Some(Self::Mmcif),
            "generated_modelcif_batches" => Some(Self::ModelCif),
            "generated_pdb_batches" => Some(Self::Pdb),
            _ => None,
        }
    }

    const fn name(self) -> &'static str {
        match self {
            Self::Mmcif => "generated_mmcif_batches",
            Self::ModelCif => "generated_modelcif_batches",
            Self::Pdb => "generated_pdb_batches",
        }
    }

    const fn layout(self) -> (&'static [u8], &'static [u8]) {
        match self {
            Self::Mmcif | Self::ModelCif => (CIF_HEADER, CIF_ROW),
            Self::Pdb => (&[], PDB_ROW),
        }
    }
}

pub(super) fn run(format: GeneratedFormat, minimum_bytes: u64) -> Result<ResourceRecord, String> {
    let name = format.name();
    measure_case(name, move || {
        let context = ExecutionContext::builder()
            .scratch_policy(ScratchPolicy::new(0))
            .build()
            .map_err(|error| format!("{name}: context failed: {error}"))?;
        let (header, row) = format.layout();
        let source = RepeatedSource::new(header, row, minimum_bytes, WINDOW_BYTES, &context)?;
        let options = ReadOptions::new();
        let dataset = DatasetId::new(0);
        let chunk = ChunkId::new(0);
        let logical_row = LogicalRow::new(0);
        let rows = match format {
            GeneratedFormat::Mmcif => drain(
                molframe::cif::MmcifBatchSource::new(
                    source,
                    options,
                    dataset,
                    chunk,
                    logical_row,
                    WINDOW_BYTES,
                    &context,
                )
                .map_err(|error| format!("{name}: open failed: {error}"))?,
                &context,
                name,
            )?,
            GeneratedFormat::ModelCif => drain(
                molframe::modelcif::ModelCifBatchSource::new(
                    source,
                    options,
                    dataset,
                    chunk,
                    logical_row,
                    WINDOW_BYTES,
                    &context,
                )
                .map_err(|error| format!("{name}: open failed: {error}"))?,
                &context,
                name,
            )?,
            GeneratedFormat::Pdb => drain(
                molframe::pdb::PdbBatchSource::new(
                    source,
                    options,
                    dataset,
                    chunk,
                    logical_row,
                    WINDOW_BYTES,
                )
                .map_err(|error| format!("{name}: open failed: {error}"))?,
                &context,
                name,
            )?,
        };
        Ok(black_box(rows))
    })
}

pub(super) fn drain<S>(mut source: S, context: &ExecutionContext, name: &str) -> Result<u64, String>
where
    S: BatchSource<Batch = StructureBatch, Error = StructureBatchError>,
{
    let mut rows = 0_u64;
    loop {
        match source
            .next_batch(BatchDemand::new(131_072, BATCH_BYTES), context)
            .map_err(|error| format!("{name}: batch failed: {error}"))?
        {
            Backpressure::Ready(batch) => {
                let count = u64::try_from(batch.batch().models().len())
                    .map_err(|_| format!("{name}: batch length exceeds u64"))?;
                rows = rows
                    .checked_add(count)
                    .ok_or_else(|| format!("{name}: row count exceeds u64"))?;
            }
            Backpressure::Pending => return Err(format!("{name}: unexpected backpressure")),
            Backpressure::Finished => return Ok(rows),
        }
    }
}

#[derive(Debug)]
struct RepeatedSource {
    header: &'static [u8],
    row: &'static [u8],
    rows: u64,
    length: u64,
    buffer: Vec<u8>,
    _reservation: MemoryReservation,
}

impl RepeatedSource {
    fn new(
        header: &'static [u8],
        row: &'static [u8],
        minimum_bytes: u64,
        window_bytes: usize,
        context: &ExecutionContext,
    ) -> Result<Self, String> {
        let header_bytes = u64::try_from(header.len()).map_err(|_| "header exceeds u64")?;
        let row_bytes = u64::try_from(row.len()).map_err(|_| "row exceeds u64")?;
        let payload = minimum_bytes.saturating_sub(header_bytes);
        let rows = payload
            .checked_add(row_bytes.saturating_sub(1))
            .and_then(|value| value.checked_div(row_bytes))
            .map_or(1, |value| value.max(1));
        let length = rows
            .checked_mul(row_bytes)
            .and_then(|value| value.checked_add(header_bytes))
            .ok_or_else(|| "generated source length exceeds u64".to_owned())?;
        let reservation = context
            .try_reserve(window_bytes)
            .map_err(|error| format!("generated source window failed: {error}"))?;
        let mut buffer = Vec::new();
        buffer
            .try_reserve_exact(window_bytes)
            .map_err(|error| format!("generated source allocation failed: {error}"))?;
        Ok(Self {
            header,
            row,
            rows,
            length,
            buffer,
            _reservation: reservation,
        })
    }

    fn append_at(&mut self, mut offset: u64, target: usize) -> Result<(), String> {
        while self.buffer.len() < target {
            let header_bytes =
                u64::try_from(self.header.len()).map_err(|_| "header length exceeds u64")?;
            let (source, local_u64) = if offset < header_bytes {
                (self.header, offset)
            } else {
                let row_bytes =
                    u64::try_from(self.row.len()).map_err(|_| "row length exceeds u64")?;
                (self.row, (offset - header_bytes) % row_bytes)
            };
            let local = usize::try_from(local_u64).map_err(|_| "local offset exceeds usize")?;
            let remaining = target.saturating_sub(self.buffer.len());
            let count = remaining.min(source.len().saturating_sub(local));
            let end = local.saturating_add(count);
            let bytes = source
                .get(local..end)
                .ok_or_else(|| "generated source range is invalid".to_owned())?;
            self.buffer.extend_from_slice(bytes);
            offset = offset
                .checked_add(u64::try_from(count).map_err(|_| "copy exceeds u64")?)
                .ok_or_else(|| "generated source offset exceeds u64".to_owned())?;
        }
        Ok(())
    }
}

impl SourceBytes for RepeatedSource {
    fn window(&mut self, start: u64, len: usize) -> Result<ByteWindow<'_>, molframe::Diagnostic> {
        if start >= self.length || len == 0 {
            return Ok(ByteWindow::new(start, &[]));
        }
        if len > self.buffer.capacity() {
            return Err(molframe::Diagnostic::new(molframe::Code::E1902)
                .with_context("requested bytes", len.to_string()));
        }
        let remaining = match usize::try_from(self.length - start) {
            Ok(length) => length,
            Err(_) => usize::MAX,
        };
        let target = len.min(remaining);
        self.buffer.clear();
        self.append_at(start, target).map_err(|error| {
            molframe::Diagnostic::new(molframe::Code::E1903).with_context("reason", error)
        })?;
        Ok(ByteWindow::new(start, &self.buffer))
    }

    fn len_hint(&self) -> Option<u64> {
        let row_bytes = u64::try_from(self.row.len()).ok()?;
        let header_bytes = u64::try_from(self.header.len()).ok()?;
        self.rows
            .checked_mul(row_bytes)
            .and_then(|value| value.checked_add(header_bytes))
    }
}

#[cfg(test)]
#[path = "generated_structure_tests.rs"]
mod tests;
