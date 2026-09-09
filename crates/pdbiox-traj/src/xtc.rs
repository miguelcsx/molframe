//! GROMACS XTC compressed-coordinate trajectory support.

use std::io::Cursor;
use std::panic::{AssertUnwindSafe, catch_unwind};

use pdbiox_core::structure::UnitCell;

use crate::Timestep;
use crate::cell::{cell_from_vectors, vectors_from_cell};
use crate::numeric::f32_from_f64;

const NM_TO_ANGSTROM: f32 = 10.0;
const ANGSTROM_TO_NM: f32 = 1.0 / NM_TO_ANGSTROM;

#[path = "xtc/reader.rs"]
mod reader;

pub use reader::XtcReader;

/// Parsed XTC frames and the encoding metadata not carried by [`Timestep`].
#[derive(Clone, Debug, PartialEq)]
pub struct XtcTrajectory {
    /// Frames in file order.
    pub frames: Vec<Timestep>,
    /// Simulation step stored in each frame.
    pub steps: Vec<u32>,
    /// Coordinate quantization in inverse nanometres for each frame.
    pub precision: Vec<f32>,
}

/// Controls deterministic XTC encoding.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct XtcWriteOptions {
    /// Coordinate quantization in inverse nanometres.
    pub precision: f32,
}

impl Default for XtcWriteOptions {
    fn default() -> Self {
        Self { precision: 1_000.0 }
    }
}

/// Malformed or unsupported XTC data.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum XtcError {
    /// The XDR stream could not be decoded or encoded.
    #[error("XTC I/O error: {0}")]
    Io(#[from] std::io::Error),
    /// A frame disagrees with the trajectory shape or contains invalid values.
    #[error("invalid XTC frame")]
    InvalidFrame,
    /// The decoder rejected contradictory internal XTC fields.
    #[error("inconsistent XTC record")]
    InconsistentRecord,
    /// A requested memory ceiling is zero or exceeds the process hard limit.
    #[error("invalid XTC memory limit {requested}; it must be at least one byte")]
    InvalidMemoryLimit {
        /// Caller-provided ceiling.
        requested: usize,
    },
    /// Reader workspaces, its index and one output frame exceed the ceiling.
    #[error("XTC frame storage requires {required} bytes, over the {limit} byte limit")]
    MemoryLimit {
        /// Required bytes, or `usize::MAX` when dimensions overflow.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
}

/// Parses all frames from an XTC byte stream.
///
/// Coordinates and periodic cells are converted from nanometres to ångström.
///
/// # Errors
///
/// Returns an error for empty, truncated, contradictory, non-finite or
/// atom-count-changing input. No partial trajectory is published.
pub fn parse_xtc(bytes: &[u8]) -> Result<XtcTrajectory, XtcError> {
    if bytes.is_empty() {
        return Err(XtcError::InvalidFrame);
    }
    let mut reader = molly::XTCReader::new(Cursor::new(bytes));
    let mut frames = Vec::new();
    let mut steps = Vec::new();
    let mut precisions = Vec::new();
    let mut atoms = None;
    while usize::try_from(reader.file.position()).map_err(|_| XtcError::InvalidFrame)? < bytes.len()
    {
        let mut source = molly::Frame::default();
        let decoded = catch_unwind(AssertUnwindSafe(|| reader.read_frame(&mut source)))
            .map_err(|_| XtcError::InconsistentRecord)?;
        decoded?;
        let positions = positions_from_source(&source)?;
        match atoms {
            Some(expected) if expected != positions.len() => return Err(XtcError::InvalidFrame),
            None => atoms = Some(positions.len()),
            Some(_) => {}
        }
        let frame = frames.len();
        let cell = cell_from_box(source.boxvec)?;
        let time = f64::from(source.time);
        if !time.is_finite() || !source.precision.is_finite() || source.precision <= 0.0 {
            return Err(XtcError::InvalidFrame);
        }
        frames.push(Timestep {
            frame,
            time: Some(time),
            positions,
            cell,
            ..Timestep::default()
        });
        steps.push(source.step);
        precisions.push(source.precision);
    }
    populate_dt(&mut frames);
    Ok(XtcTrajectory {
        frames,
        steps,
        precision: precisions,
    })
}

/// Encodes frames as a GROMACS XTC byte stream.
///
/// Positions and cells are converted from ångström to nanometres. The frame's
/// ordinal is used as its simulation step and missing time is rejected.
///
/// # Errors
///
/// Refuses empty trajectories, changing atom counts, absent/non-finite time,
/// invalid cells, coordinates outside XTC step width, or invalid precision.
pub fn write_xtc(frames: &[Timestep], options: XtcWriteOptions) -> Result<Vec<u8>, XtcError> {
    write_xtc_with_precisions(frames, &vec![options.precision; frames.len()])
}

/// Encodes XTC frames while preserving each frame's quantization.
///
/// # Errors
///
/// Returns the same validation errors as [`write_xtc`] and rejects a precision
/// count that differs from the number of frames.
pub fn write_xtc_with_precisions(
    frames: &[Timestep],
    precisions: &[f32],
) -> Result<Vec<u8>, XtcError> {
    let atoms = frames
        .first()
        .map(|frame| frame.positions.len())
        .filter(|atoms| *atoms > 0)
        .ok_or(XtcError::InvalidFrame)?;
    if precisions.len() != frames.len()
        || precisions
            .iter()
            .any(|precision| !precision.is_finite() || *precision <= 0.0)
    {
        return Err(XtcError::InvalidFrame);
    }
    let mut writer = molly::XTCWriter::new(Cursor::new(Vec::new()));
    for (frame, precision) in frames.iter().zip(precisions) {
        validate_frame(frame, atoms, *precision)?;
        let step = u32::try_from(frame.frame).map_err(|_| XtcError::InvalidFrame)?;
        let time = frame
            .time
            .and_then(f32_from_f64)
            .ok_or(XtcError::InvalidFrame)?;
        let boxvec = box_from_cell(frame)?;
        let positions = frame
            .positions
            .iter()
            .flat_map(|position| position.map(|value| value * ANGSTROM_TO_NM))
            .collect();
        writer.write_frame(&molly::Frame {
            step,
            time,
            boxvec,
            precision: *precision,
            positions,
        })?;
    }
    Ok(writer.file.into_inner())
}

fn positions_from_source(source: &molly::Frame) -> Result<Vec<[f32; 3]>, XtcError> {
    if source.positions.is_empty()
        || !source.positions.len().is_multiple_of(3)
        || source.positions.iter().any(|value| !value.is_finite())
    {
        return Err(XtcError::InvalidFrame);
    }
    Ok(source
        .positions
        .chunks_exact(3)
        .map(|values| {
            [
                values[0] * NM_TO_ANGSTROM,
                values[1] * NM_TO_ANGSTROM,
                values[2] * NM_TO_ANGSTROM,
            ]
        })
        .collect())
}

fn cell_from_box(boxvec: [f32; 9]) -> Result<Option<UnitCell>, XtcError> {
    if boxvec.iter().all(|value| *value == 0.0) {
        return Ok(None);
    }
    if boxvec.iter().any(|value| !value.is_finite()) {
        return Err(XtcError::InvalidFrame);
    }
    let mut vectors = [[0.0; 3]; 3];
    for (target, source) in vectors.iter_mut().zip(boxvec.chunks_exact(3)) {
        for (value, input) in target.iter_mut().zip(source) {
            *value = f64::from(*input * NM_TO_ANGSTROM);
        }
    }
    cell_from_vectors(vectors)
        .map(Some)
        .ok_or(XtcError::InvalidFrame)
}

fn box_from_cell(frame: &Timestep) -> Result<[f32; 9], XtcError> {
    let Some(cell) = frame.cell else {
        return Ok([0.0; 9]);
    };
    let vectors = vectors_from_cell(cell).ok_or(XtcError::InvalidFrame)?;
    let mut boxvec = [0.0; 9];
    for (target, value) in boxvec.iter_mut().zip(vectors.into_iter().flatten()) {
        *target = f32_from_f64(value).ok_or(XtcError::InvalidFrame)? * ANGSTROM_TO_NM;
    }
    Ok(boxvec)
}

fn validate_frame(frame: &Timestep, atoms: usize, precision: f32) -> Result<(), XtcError> {
    let valid = frame.positions.len() == atoms
        && frame.time.is_some_and(f64::is_finite)
        && frame.positions.iter().flatten().all(|value| {
            value.is_finite()
                && (f64::from(*value) * f64::from(ANGSTROM_TO_NM) * f64::from(precision)).abs()
                    <= f64::from(i32::MAX)
        });
    if valid {
        Ok(())
    } else {
        Err(XtcError::InvalidFrame)
    }
}

fn populate_dt(frames: &mut [Timestep]) {
    if frames.len() < 2 {
        return;
    }
    let dt = frames[1]
        .time
        .zip(frames[0].time)
        .map(|(right, left)| right - left);
    if dt.is_some_and(f64::is_finite) {
        for frame in frames {
            frame.dt = dt;
        }
    }
}

#[cfg(test)]
#[path = "xtc_tests.rs"]
mod tests;
