//! HOOMD GSD trajectories with explicit caller-provided length units.

use std::path::Path;

use hoomd_gsd::file_layer::{GsdFile, Mode};
use hoomd_gsd::hoomd::HoomdGsdFile;
use pdbiox_core::structure::UnitCell;

use crate::Timestep;
use crate::cell::{cell_from_vectors, vectors_from_cell};
use crate::numeric::f32_from_f64;

const SCHEMA: &str = "hoomd";
const STEP: &str = "configuration/step";
const BOX: &str = "configuration/box";
const POSITION: &str = "particles/position";

/// Required unit conversion for a GSD file, whose schema is unit-agnostic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GsdOptions {
    /// Number of ångström represented by one file length unit.
    pub length_to_angstrom: f64,
}

impl GsdOptions {
    /// Validates and constructs an explicit length conversion.
    ///
    /// # Errors
    ///
    /// Refuses non-finite and non-positive scales.
    pub fn new(length_to_angstrom: f64) -> Result<Self, GsdError> {
        if length_to_angstrom.is_finite() && length_to_angstrom > 0.0 {
            Ok(Self { length_to_angstrom })
        } else {
            Err(GsdError::InvalidUnits)
        }
    }
}

/// Parsed GSD trajectory and its simulation steps.
#[derive(Clone, Debug, PartialEq)]
pub struct GsdTrajectory {
    /// Frames in file order.
    pub frames: Vec<Timestep>,
    /// Simulation step stored in each frame.
    pub steps: Vec<u64>,
}

/// GSD container or HOOMD schema failure.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum GsdError {
    /// Unit-agnostic GSD requires a finite positive caller scale.
    #[error("GSD requires an explicit positive length conversion")]
    InvalidUnits,
    /// The container could not be opened or decoded.
    #[error("GSD container error: {0}")]
    Container(String),
    /// The file does not use the supported HOOMD schema.
    #[error("GSD file is not a HOOMD-schema trajectory")]
    Schema,
    /// A required chunk is missing or has the wrong shape/type.
    #[error("invalid HOOMD GSD frame")]
    InvalidFrame,
}

/// Reads a HOOMD-schema GSD trajectory from `path`.
///
/// HOOMD permits unchanged chunks to be omitted after frame zero; this reader
/// applies that defined inheritance rule. It never guesses the file's units.
///
/// # Errors
///
/// Returns an error for invalid units, schema, chunks, values, cells or changing
/// particle counts.
pub fn parse_gsd(path: &Path, options: GsdOptions) -> Result<GsdTrajectory, GsdError> {
    validate_options(options)?;
    let file = GsdFile::open(path, Mode::Read).map_err(container_error)?;
    if file.schema() != SCHEMA || file.n_frames() == 0 {
        return Err(GsdError::Schema);
    }
    let mut frames = Vec::new();
    let mut steps = Vec::new();
    let mut atoms = None;
    for frame_index in 0..file.n_frames() {
        let step = read_step(&file, frame_index)?;
        let positions = read_positions(&file, frame_index, options)?;
        match atoms {
            Some(expected) if expected != positions.len() => return Err(GsdError::InvalidFrame),
            None => atoms = Some(positions.len()),
            Some(_) => {}
        }
        frames.push(Timestep {
            frame: usize::try_from(frame_index).map_err(|_| GsdError::InvalidFrame)?,
            positions,
            cell: read_cell(&file, frame_index, options)?,
            ..Timestep::default()
        });
        steps.push(step);
    }
    Ok(GsdTrajectory { frames, steps })
}

/// Writes a HOOMD-schema GSD trajectory to `path`.
///
/// # Errors
///
/// Refuses invalid units, empty trajectories, changing atom counts, invalid
/// values/cells and frame ordinals that do not fit the GSD step field.
pub fn write_gsd(path: &Path, frames: &[Timestep], options: GsdOptions) -> Result<(), GsdError> {
    validate_options(options)?;
    let atoms = frames
        .first()
        .map(|frame| frame.positions.len())
        .filter(|count| *count > 0)
        .ok_or(GsdError::InvalidFrame)?;
    let mut file = HoomdGsdFile::create(path).map_err(container_error)?;
    for frame in frames {
        validate_frame(frame, atoms)?;
        let step = u64::try_from(frame.frame).map_err(|_| GsdError::InvalidFrame)?;
        let positions: Vec<(f64, f64, f64)> = frame
            .positions
            .iter()
            .map(|position| {
                (
                    f64::from(position[0]) / options.length_to_angstrom,
                    f64::from(position[1]) / options.length_to_angstrom,
                    f64::from(position[2]) / options.length_to_angstrom,
                )
            })
            .collect();
        let mut output = file.append_frame(step).map_err(container_error)?;
        if let Some(cell) = frame.cell {
            output = output
                .configuration_box(box_parameters(cell, options)?)
                .map_err(container_error)?;
        }
        output
            .particles_position(positions.into_iter().map(Into::into))
            .map_err(container_error)?
            .end()
            .map_err(container_error)?;
    }
    file.sync_all().map_err(container_error)
}

fn read_step(file: &GsdFile, frame: u64) -> Result<u64, GsdError> {
    let values: Vec<_> = file
        .iter_scalars::<u64>(frame, STEP)
        .map_err(container_error)?
        .collect();
    match values.as_slice() {
        [step] => Ok(*step),
        _ => Err(GsdError::InvalidFrame),
    }
}

fn inherited_frame(file: &GsdFile, frame: u64, chunk: &str) -> Option<u64> {
    file.find_chunk(frame, chunk)
        .map(|_| frame)
        .or_else(|| file.find_chunk(0, chunk).map(|_| 0))
}

fn read_positions(
    file: &GsdFile,
    frame: u64,
    options: GsdOptions,
) -> Result<Vec<[f32; 3]>, GsdError> {
    let source = inherited_frame(file, frame, POSITION).ok_or(GsdError::InvalidFrame)?;
    let scale = f32_from_f64(options.length_to_angstrom).ok_or(GsdError::InvalidFrame)?;
    let values: Vec<_> = file
        .iter_arrays::<f32, 3>(source, POSITION)
        .map_err(container_error)?
        .map(|position| position.map(|value| value * scale))
        .collect();
    if values.is_empty() || values.iter().flatten().any(|value| !value.is_finite()) {
        Err(GsdError::InvalidFrame)
    } else {
        Ok(values)
    }
}

fn read_cell(
    file: &GsdFile,
    frame: u64,
    options: GsdOptions,
) -> Result<Option<UnitCell>, GsdError> {
    let Some(source) = inherited_frame(file, frame, BOX) else {
        return Ok(None);
    };
    let values: Vec<_> = file
        .iter_scalars::<f32>(source, BOX)
        .map_err(container_error)?
        .map(f64::from)
        .collect();
    let [lx, ly, lz, xy, xz, yz] = values.try_into().map_err(|_| GsdError::InvalidFrame)?;
    let scale = options.length_to_angstrom;
    cell_from_vectors([
        [lx * scale, 0.0, 0.0],
        [xy * ly * scale, ly * scale, 0.0],
        [xz * lz * scale, yz * lz * scale, lz * scale],
    ])
    .map(Some)
    .ok_or(GsdError::InvalidFrame)
}

fn box_parameters(cell: UnitCell, options: GsdOptions) -> Result<[f64; 6], GsdError> {
    let vectors = vectors_from_cell(cell).ok_or(GsdError::InvalidFrame)?;
    let lx = vectors[0][0] / options.length_to_angstrom;
    let ly = vectors[1][1] / options.length_to_angstrom;
    let lz = vectors[2][2] / options.length_to_angstrom;
    if ly <= 0.0 || lz <= 0.0 {
        return Err(GsdError::InvalidFrame);
    }
    Ok([
        lx,
        ly,
        lz,
        vectors[1][0] / vectors[1][1],
        vectors[2][0] / vectors[2][2],
        vectors[2][1] / vectors[2][2],
    ])
}

fn validate_options(options: GsdOptions) -> Result<(), GsdError> {
    GsdOptions::new(options.length_to_angstrom).map(|_| ())
}

fn validate_frame(frame: &Timestep, atoms: usize) -> Result<(), GsdError> {
    if frame.positions.len() == atoms
        && frame.time.is_none()
        && frame.dt.is_none()
        && frame.velocities.is_none()
        && frame.forces.is_none()
        && frame.data.is_empty()
        && frame
            .positions
            .iter()
            .flatten()
            .all(|value| value.is_finite())
    {
        Ok(())
    } else {
        Err(GsdError::InvalidFrame)
    }
}

fn container_error(error: impl std::fmt::Display) -> GsdError {
    GsdError::Container(error.to_string())
}

#[cfg(test)]
#[path = "gsd_tests.rs"]
mod tests;
