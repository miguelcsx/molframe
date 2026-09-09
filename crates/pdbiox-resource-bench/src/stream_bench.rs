//! Reading a generated mmCIF stream, at sizes past any committed fixture.
//!
//! These cases exist to answer one question: does peak memory follow the input
//! size, or the caller's budget? The answer today is the input size, which is
//! what bounds pdbiox to files that fit in memory. The measurement is recorded
//! at several sizes so the shape of the curve is visible — a flat line across
//! two decimal orders is the claim, and a straight line at slope one is the
//! defect.
//!
//! Generated inputs are materialised before the measured region, then consumed
//! through the public bounded batch API. This isolates ingest memory from both
//! fixture generation and optional full-structure collection.

use super::{ResourceRecord, measure_case, measure_retained_case};
use pdbiox::core::{
    Backpressure, BatchDemand, BatchSource, ExecutionContext, MemoryBudget, ScratchPolicy,
    SourceBytes, WindowedFile,
};
use pdbiox::{InputBuffer, Limits, ReadOptions};
use pdbiox_bench::{Sample, Seed, SyntheticCifSource, Tile, structure};
use std::hint::black_box;
use std::io::Write;
use std::path::Path;

/// Bytes in a gibibyte.
const GIBIBYTE: u64 = 1 << 30;
const WINDOW_BYTES: usize = 8 << 20;
const BATCH_BYTES: usize = 16 << 20;

pub(super) fn run_stream_synthetic_256m() -> Result<ResourceRecord, String> {
    read_stream("stream_synthetic_256m", GIBIBYTE / 4)
}

pub(super) fn run_stream_synthetic_1g() -> Result<ResourceRecord, String> {
    read_stream("stream_synthetic_1g", GIBIBYTE)
}

/// Reads roughly `bytes` of generated mmCIF and reports what it cost.
///
/// The tile and the copy count are prepared before the measurement opens, so
/// the record covers the read and nothing else.
fn read_stream(name: &'static str, bytes: u64) -> Result<ResourceRecord, String> {
    let tile = Tile::from_structure(&structure(Sample::Large), Seed::new(1));
    let copies = SyntheticCifSource::copies_for_bytes(&tile, bytes);
    if copies == 0 {
        return Err(format!("{name}: the tile produced no bytes"));
    }

    let directory =
        tempfile::tempdir().map_err(|error| format!("temporary directory failed: {error}"))?;
    let path = directory.path().join("synthetic.cif");
    let mut source = SyntheticCifSource::new(tile, copies);
    let mut file = std::fs::File::create(&path)
        .map_err(|error| format!("{name}: fixture create failed: {error}"))?;
    std::io::copy(&mut source, &mut file)
        .map_err(|error| format!("{name}: fixture generation failed: {error}"))?;
    file.flush()
        .map_err(|error| format!("{name}: fixture flush failed: {error}"))?;

    measure_case(name, move || drain_structure_batches(name, &path))
}

pub(super) fn drain_structure_batches(name: &'static str, path: &Path) -> Result<u64, String> {
    let context = ExecutionContext::builder()
        .scratch_policy(ScratchPolicy::new(0))
        .build()
        .map_err(|error| format!("{name}: context failed: {error}"))?;
    drain_structure_batches_with_context(name, path, &context).map(|stats| stats.rows)
}

#[derive(Clone, Copy)]
struct BatchStats {
    rows: u64,
    models: u64,
}

fn drain_structure_batches_with_context(
    name: &'static str,
    path: &Path,
    context: &ExecutionContext,
) -> Result<BatchStats, String> {
    let mut source = pdbiox::open_structure_batches(path, &ReadOptions::new(), context)
        .map_err(|error| format!("{name}: bounded open failed: {error}"))?;
    let mut rows = 0_u64;
    let mut models = 0_u64;
    let mut previous_model = None;
    loop {
        match source
            .next_batch(BatchDemand::new(131_072, BATCH_BYTES), context)
            .map_err(|error| format!("{name}: batch failed: {error}"))?
        {
            Backpressure::Ready(batch) => {
                for model in batch.batch().models() {
                    if previous_model != Some(*model) {
                        models = models
                            .checked_add(1)
                            .ok_or_else(|| format!("{name}: model count exceeds u64"))?;
                        previous_model = Some(*model);
                    }
                }
                let count = u64::try_from(batch.batch().models().len())
                    .map_err(|_| format!("{name}: local batch length exceeds u64"))?;
                rows = rows
                    .checked_add(count)
                    .ok_or_else(|| format!("{name}: row count exceeds u64"))?;
            }
            Backpressure::Pending => {
                return Err(format!(
                    "{name}: source made no progress without a live lease"
                ));
            }
            Backpressure::Finished => break,
        }
    }
    Ok(black_box(BatchStats { rows, models }))
}

pub(super) fn run_structure_batch_file(
    path: &Path,
    maximum_spill_bytes: u64,
) -> Result<ResourceRecord, String> {
    measure_retained_case("structure_batch_file", || {
        let directory = tempfile::tempdir()
            .map_err(|error| format!("compressed spill directory failed: {error}"))?;
        let context = ExecutionContext::builder()
            .scratch_policy(ScratchPolicy::new(0))
            .temp_storage_policy(pdbiox::core::TempStoragePolicy::directory(
                directory.path(),
                maximum_spill_bytes,
            ))
            .build()
            .map_err(|error| format!("structure_batch_file: context failed: {error}"))?;
        let stats = drain_structure_batches_with_context("structure_batch_file", path, &context)?;
        Ok((stats.rows, Some(stats.models), ()))
    })
}

/// Writes a generated mmCIF to `path` and returns its byte length.
fn materialise(path: &Path, bytes: u64) -> Result<u64, String> {
    let tile = Tile::from_structure(&structure(Sample::Large), Seed::new(1));
    let copies = SyntheticCifSource::copies_for_bytes(&tile, bytes);
    let mut source = SyntheticCifSource::new(tile, copies);
    let mut file =
        std::fs::File::create(path).map_err(|error| format!("create failed: {error}"))?;
    let mut buffer = vec![0_u8; 1 << 20];
    let mut written = 0_u64;
    loop {
        let count = match std::io::Read::read(&mut source, &mut buffer) {
            Ok(0) => break,
            Ok(count) => count,
            Err(error) => return Err(format!("generation failed: {error}")),
        };
        let Some(chunk) = buffer.get(..count) else {
            break;
        };
        file.write_all(chunk)
            .map_err(|error| format!("write failed: {error}"))?;
        match u64::try_from(count) {
            Ok(count) => written += count,
            Err(_) => return Err("write count does not fit a u64".to_owned()),
        }
    }
    file.flush()
        .map_err(|error| format!("flush failed: {error}"))?;
    Ok(written)
}

/// Reads a materialised file through the copying open path.
pub(super) fn run_file_copied_1g() -> Result<ResourceRecord, String> {
    let directory =
        tempfile::tempdir().map_err(|error| format!("temporary directory failed: {error}"))?;
    let path = directory.path().join("synthetic.cif");
    let _bytes = materialise(&path, GIBIBYTE)?;

    measure_case("file_copied_1g", || {
        let buffer = InputBuffer::open(&path, unbounded())
            .map_err(|finding| format!("open failed: {finding:?}"))?;
        let (structure, _findings) = pdbiox::cif::read(&buffer, &ReadOptions::new())
            .map_err(|findings| format!("read failed: {findings:?}"))?;
        Ok(black_box(u64::from(structure.atom_count())))
    })
}

/// Reads the same file through the zero-copy mapped path.
pub(super) fn run_file_mapped_1g() -> Result<ResourceRecord, String> {
    let directory =
        tempfile::tempdir().map_err(|error| format!("temporary directory failed: {error}"))?;
    let path = directory.path().join("synthetic.cif");
    let _bytes = materialise(&path, GIBIBYTE)?;

    measure_case("file_mapped_1g", || {
        let file = std::fs::File::open(&path).map_err(|error| format!("open failed: {error}"))?;
        // SAFETY: the file was written by this process into a temporary
        // directory it owns, is not open for writing anywhere, and is removed
        // only after the mapping is dropped.
        let mapped = unsafe { pdbiox_mmap::MappedFile::map_file_unchecked(&file) }
            .map_err(|error| format!("map failed: {error}"))?;
        let buffer = InputBuffer::from_mapped(mapped);
        let (structure, _findings) = pdbiox::cif::read(&buffer, &ReadOptions::new())
            .map_err(|findings| format!("read failed: {findings:?}"))?;
        Ok(black_box(u64::from(structure.atom_count())))
    })
}

/// Reads every byte through one reusable budgeted window.
///
/// This is the long-gate primitive used for 1 TiB forward-read validation. It
/// deliberately performs no output materialisation, so retained memory must be
/// independent of file length.
pub(super) fn run_window_scan_file(path: &Path) -> Result<ResourceRecord, String> {
    measure_case("window_scan_file", || {
        let budget = MemoryBudget::new(WINDOW_BYTES * 2)
            .map_err(|error| format!("window budget failed: {error}"))?;
        let context = ExecutionContext::builder()
            .memory_budget(budget)
            .scratch_policy(ScratchPolicy::new(0))
            .build()
            .map_err(|error| format!("window context failed: {error}"))?;
        let mut source = WindowedFile::open(path, WINDOW_BYTES, &context)
            .map_err(|error| format!("windowed open failed: {error}"))?;
        let expected = source
            .len_hint()
            .ok_or_else(|| "windowed file did not publish its length".to_owned())?;
        let mut offset = 0_u64;
        let mut digest = 0_u64;
        while offset < expected {
            let window = source
                .window(offset, WINDOW_BYTES)
                .map_err(|error| format!("window at {offset} failed: {error}"))?;
            if window.bytes().is_empty() {
                return Err(format!("windowed scan stopped at {offset} of {expected}"));
            }
            let bytes = window.bytes();
            let first = bytes.first().copied().map_or(0, u64::from);
            let last = bytes.last().copied().map_or(0, u64::from);
            digest = digest.rotate_left(7)
                ^ first
                ^ (last << 8)
                ^ u64::try_from(bytes.len()).map_err(|_| "window length exceeds u64")?;
            offset = window
                .end()
                .map_err(|error| format!("window end at {offset} failed: {error}"))?;
        }
        black_box(digest);
        Ok(offset)
    })
}

/// Counts coordinate cells without building anything sized by the input.
///
/// Isolating the reader from the structure it produces is the only way to ask
/// whether ingest is bounded. A read that materialises a structure has a
/// footprint proportional to the atom count however the bytes arrived, so
/// measuring it answers a question about the output rather than the input.
struct CountingSink {
    cells: u64,
}

impl pdbiox::cif::CifEventSink for CountingSink {
    type Output = u64;

    fn block(&mut self, _name: &str) {}

    fn value(
        &mut self,
        _category: &str,
        _item: &str,
        _value: pdbiox::cif::CifScalar<'_>,
        _span: pdbiox::ByteSpan,
    ) {
        self.cells = self.cells.saturating_add(1);
    }

    fn finish(self) -> u64 {
        self.cells
    }
}

pub(super) fn run_lex_synthetic_1g() -> Result<ResourceRecord, String> {
    lex_stream("lex_synthetic_1g", GIBIBYTE)
}

/// Runs the lexer alone over roughly `bytes` of generated mmCIF.
///
/// Timed against the scan case, which adds value typing and event dispatch on
/// top of the same traversal, this separates what tokenising costs from what
/// interpreting the tokens costs. Static reading is a poor guide to that split.
fn lex_stream(name: &'static str, bytes: u64) -> Result<ResourceRecord, String> {
    let tile = Tile::from_structure(&structure(Sample::Large), Seed::new(1));
    let copies = SyntheticCifSource::copies_for_bytes(&tile, bytes);

    measure_case(name, move || {
        let source = SyntheticCifSource::new(tile, copies);
        let buffer = InputBuffer::from_reader(source, unbounded())
            .map_err(|finding| format!("{name}: input was refused: {finding:?}"))?;
        let mut lexer = pdbiox::cif::lexer::Lexer::new(buffer.as_bytes())
            .map_err(|error| format!("{name}: lexer refused the input: {error:?}"))?;
        let mut tokens = 0_u64;
        loop {
            match lexer.next_token() {
                Ok(Some(_)) => tokens = tokens.saturating_add(1),
                Ok(None) => break,
                Err(error) => return Err(format!("{name}: lexing failed: {error:?}")),
            }
        }
        Ok(black_box(tokens))
    })
}

pub(super) fn run_scan_synthetic_256m() -> Result<ResourceRecord, String> {
    scan_stream("scan_synthetic_256m", GIBIBYTE / 4)
}

pub(super) fn run_scan_synthetic_1g() -> Result<ResourceRecord, String> {
    scan_stream("scan_synthetic_1g", GIBIBYTE)
}

/// Parses roughly `bytes` of generated mmCIF, retaining nothing per row.
fn scan_stream(name: &'static str, bytes: u64) -> Result<ResourceRecord, String> {
    let tile = Tile::from_structure(&structure(Sample::Large), Seed::new(1));
    let copies = SyntheticCifSource::copies_for_bytes(&tile, bytes);

    measure_case(name, move || {
        let source = SyntheticCifSource::new(tile, copies);
        let buffer = InputBuffer::from_reader(source, unbounded())
            .map_err(|finding| format!("{name}: input was refused: {finding:?}"))?;
        let (cells, _findings) = pdbiox::cif::parse_events(&buffer, CountingSink { cells: 0 })
            .map_err(|findings| format!("{name}: parse failed: {findings:?}"))?;
        Ok(black_box(cells))
    })
}

/// Read limits with the decompressed-byte ceiling lifted.
fn unbounded() -> Limits {
    Limits {
        decompressed_bytes: u64::MAX,
        ..Limits::default()
    }
}
