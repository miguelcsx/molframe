//! Parallel in-memory XTC decoding.
//!
//! One scan reads only each frame's 56-byte header and, for compressed
//! frames, the payload length word, so it costs `O(frames)` and never touches
//! the compressed payload. Frames are then decoded in fixed blocks on the
//! worker pool, each block with its own decoder over the shared borrowed
//! bytes. Total work is `O(bytes)` divided over the pool; memory is the
//! output plus one decoder workspace per in-flight block.
//!
//! Anything the scan cannot vouch for (a bad header, a length that overruns
//! the input, a changed atom count) is handed to the sequential decoder as a
//! tail, so errors are exactly those of a frame-by-frame read.

use std::io::Cursor;
use std::ops::Range;
use std::panic::{AssertUnwindSafe, catch_unwind};

use molframe_core::ExecutionContext;
use molframe_core::parallel::{BlockPlan, map_blocks_in};

use super::{XtcError, cell_from_box, positions_from_source};
use crate::Timestep;

/// Frames decoded per parallel block.
const FRAMES_PER_BLOCK: usize = 4;

type Decoded = (Timestep, u32, f32);

/// Starts of every frame the scan could size, plus the first frame's atoms.
struct Scan {
    starts: Vec<usize>,
    atoms: usize,
    complete: bool,
}

pub(super) fn decode_all(bytes: &[u8]) -> Result<Vec<Decoded>, XtcError> {
    let scan = scan(bytes);
    let mut edges = scan.starts.clone();
    if !scan.complete || edges.is_empty() {
        // The tail starts where scanning stopped and is read sequentially.
        let resume = scan.starts.last().map_or(0, |_| bytes.len());
        if scan.starts.is_empty() {
            edges.push(0);
        } else if resume != bytes.len() {
            edges.push(resume);
        }
    }
    let units = edges.len();
    edges.push(bytes.len());
    let context = ExecutionContext::default();
    let blocks = map_blocks_in(
        BlockPlan::new(units, FRAMES_PER_BLOCK),
        &context,
        |_, range| decode_block(bytes, &edges, range, scan.atoms),
    )
    .map_err(|_| XtcError::InconsistentRecord)?;
    let mut frames = Vec::with_capacity(units);
    for result in blocks.into_iter().flatten() {
        frames.push(result?);
    }
    Ok(frames)
}

fn scan(bytes: &[u8]) -> Scan {
    let mut starts = Vec::new();
    let mut atoms = 0;
    let mut position = 0;
    while position < bytes.len() {
        let Some(end) = frame_end(bytes, position, atoms) else {
            return Scan {
                starts,
                atoms,
                complete: false,
            };
        };
        if starts.is_empty() {
            atoms = end.atoms;
        }
        starts.push(position);
        position = end.end;
    }
    Scan {
        starts,
        atoms,
        complete: true,
    }
}

struct FrameEnd {
    end: usize,
    atoms: usize,
}

/// End of the frame starting at `start`, or `None` if it cannot be sized.
fn frame_end(bytes: &[u8], start: usize, expected_atoms: usize) -> Option<FrameEnd> {
    let mut input = bytes.get(start..)?;
    let header = catch_unwind(AssertUnwindSafe(|| molly::Header::read(&mut input)))
        .ok()?
        .ok()?;
    let atoms = header.natoms;
    if atoms == 0 || (expected_atoms != 0 && atoms != expected_atoms) {
        return None;
    }
    let body = start.checked_add(molly::Header::SIZE)?;
    let end = if atoms <= 9 {
        body.checked_add(atoms.checked_mul(size_of::<[f32; 3]>())?)?
    } else {
        let count_at = body.checked_add(32)?;
        let mut tail = bytes.get(count_at..)?;
        let count = catch_unwind(AssertUnwindSafe(|| {
            molly::reader::read_nbytes(&mut tail, header.magic)
        }))
        .ok()?
        .ok()?;
        let count_bytes = match header.magic {
            molly::Magic::Xtc1995 => 4,
            molly::Magic::Xtc2023 => 8,
        };
        count_at
            .checked_add(count_bytes)?
            .checked_add(count)?
            .checked_add(molly::padding(count))?
    };
    (end <= bytes.len() && end > start).then_some(FrameEnd { end, atoms })
}

fn decode_block(
    bytes: &[u8],
    edges: &[usize],
    units: Range<usize>,
    atoms: usize,
) -> Vec<Result<Decoded, XtcError>> {
    let mut reader = molly::XTCReader::new(Cursor::new(bytes));
    reader.file.set_position(edges[units.start] as u64);
    let stop = edges[units.end];
    let mut out = Vec::with_capacity(units.len());
    let mut expected = (atoms != 0).then_some(atoms);
    let mut source = molly::Frame::default();
    while usize::try_from(reader.file.position()).is_ok_and(|position| position < stop) {
        let decoded = decode_one(&mut reader, &mut source, &mut expected);
        let failed = decoded.is_err();
        out.push(decoded);
        if failed {
            break;
        }
    }
    out
}

fn decode_one(
    reader: &mut molly::XTCReader<Cursor<&[u8]>>,
    source: &mut molly::Frame,
    expected: &mut Option<usize>,
) -> Result<Decoded, XtcError> {
    catch_unwind(AssertUnwindSafe(|| reader.read_frame(source)))
        .map_err(|_| XtcError::InconsistentRecord)??;
    let positions = positions_from_source(source)?;
    match *expected {
        Some(atoms) if atoms != positions.len() => return Err(XtcError::InvalidFrame),
        None => *expected = Some(positions.len()),
        Some(_) => {}
    }
    let cell = cell_from_box(source.boxvec)?;
    let time = f64::from(source.time);
    if !time.is_finite() || !source.precision.is_finite() || source.precision <= 0.0 {
        return Err(XtcError::InvalidFrame);
    }
    Ok((
        Timestep {
            time: Some(time),
            positions,
            cell,
            ..Timestep::default()
        },
        source.step,
        source.precision,
    ))
}
