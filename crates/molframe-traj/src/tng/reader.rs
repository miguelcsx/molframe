//! Native TNG reader with canonical-unit conversion.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;

use molframe_core::structure::UnitCell;
use tng_rs::data::Compression;
use tng_rs::gen_block::BlockID;
use tng_rs::trajectory::Trajectory;

use crate::Timestep;
use crate::cell::cell_from_vectors;
use crate::numeric::{f32_triplet, f64_from_i64};

use super::{TngCompression, TngError, TngTrajectory};

const SECONDS_TO_PICOSECONDS: f64 = 1.0e12;
const VECTOR_WIDTH: usize = 3;
const CELL_WIDTH: usize = 9;

/// Reads every coordinate-bearing frame from a TNG file.
///
/// Sparse simulation steps are preserved. Positions, velocities, forces,
/// periodic cells and time are converted to molframe canonical units. The work is
/// linear in decoded values; malformed input never publishes partial results.
///
/// # Errors
///
/// Returns an error for invalid headers, hashes, block shapes, frame ordering,
/// units, non-finite values or native-decoder invariant failures.
pub fn parse_tng(path: &Path) -> Result<TngTrajectory, TngError> {
    catch_unwind(AssertUnwindSafe(|| parse_inner(path))).map_err(|_| TngError::InternalInvariant)?
}

fn parse_inner(path: &Path) -> Result<TngTrajectory, TngError> {
    let mut source = open(path)?;
    let atoms = usize::try_from(source.num_particles_get()).map_err(|_| TngError::InvalidShape)?;
    if atoms == 0 {
        return Err(TngError::InvalidShape);
    }
    let logical_frames = source.num_frames_get()?;
    if logical_frames <= 0 {
        return Err(TngError::InvalidShape);
    }
    let steps = position_steps(&mut source, logical_frames)?;
    let (exponent, length_scale, precision, compression) = parse_metadata(&mut source, steps[0])?;
    let mut frames = read_frames(&mut source, &steps, atoms, length_scale)?;
    populate_dt(&mut frames);
    Ok(TngTrajectory {
        frames,
        steps,
        distance_unit_exponent: exponent,
        compression_precision: precision,
        compression,
    })
}

fn parse_metadata(
    source: &mut Trajectory,
    first_step: i64,
) -> Result<(i64, f64, f64, TngCompression), TngError> {
    let exponent = source.distance_unit_exponential_get();
    let length_scale = length_to_angstrom(exponent)?;
    let precision = source.compression_precision_get();
    if !precision.is_finite() || precision <= 0.0 {
        return Err(TngError::InvalidValue);
    }
    source.frame_set_of_frame_find(first_step)?;
    let (codec, block_precision) =
        source.util_frame_current_compression_get(BlockID::TrajPositions)?;
    let compression = match codec {
        Compression::Uncompressed => TngCompression::Uncompressed,
        Compression::GZip => TngCompression::Lossless,
        Compression::TNG if block_precision.is_finite() && block_precision > 0.0 => {
            TngCompression::Lossy {
                precision: block_precision,
            }
        }
        _ => return Err(TngError::InvalidValue),
    };
    Ok((exponent, length_scale, precision, compression))
}

fn read_frames(
    source: &mut Trajectory,
    steps: &[i64],
    atoms: usize,
    scale: f64,
) -> Result<Vec<Timestep>, TngError> {
    let mut frames = Vec::with_capacity(steps.len());
    let mut start = 0;
    while start < steps.len() {
        let end = frame_set_end(source, steps, start)?;
        read_frame_set(source, steps, start, end, atoms, scale, &mut frames)?;
        start = end;
    }
    Ok(frames)
}

fn frame_set_end(source: &mut Trajectory, steps: &[i64], start: usize) -> Result<usize, TngError> {
    source.frame_set_of_frame_find(steps[start])?;
    let first = source.current_trajectory_frame_set.first_frame;
    let mut end = start + 1;
    while end < steps.len() {
        source.frame_set_of_frame_find(steps[end])?;
        if source.current_trajectory_frame_set.first_frame != first {
            break;
        }
        end += 1;
    }
    Ok(end)
}

fn read_frame_set(
    source: &mut Trajectory,
    steps: &[i64],
    start: usize,
    end: usize,
    atoms: usize,
    scale: f64,
    frames: &mut Vec<Timestep>,
) -> Result<(), TngError> {
    let first_step = steps[start];
    let last_step = steps[end - 1];
    source.frame_set_of_frame_find(first_step)?;
    let first_frame = source.current_trajectory_frame_set.first_frame;
    let first_seconds = source.current_trajectory_frame_set.first_frame_time;
    let time_per_frame = source.time_per_frame;
    let count = end - start;
    let (values, stride) = source.util_pos_read_range(first_step, last_step)?;
    let positions = frame_vectors(&values, stride, atoms, scale, count)?;
    let velocities = optional_frame_vectors(
        source.util_vel_read_range(first_step, last_step),
        atoms,
        scale,
        count,
    )?;
    let forces = optional_frame_vectors(
        source.util_force_read_range(first_step, last_step),
        atoms,
        1.0 / scale,
        count,
    )?;
    let cells = optional_frame_values(
        source.util_box_shape_read_range(first_step, last_step),
        CELL_WIDTH,
        count,
    )?;
    append_frames(
        frames,
        &steps[start..end],
        Timebase {
            offset: start,
            first_frame,
            first_seconds,
            time_per_frame,
            scale,
        },
        AxisData {
            positions,
            velocities,
            forces,
            cells,
        },
    )
}

/// The per-frame axis sequences a TNG frame set contributes, consumed one
/// element each until the steps are exhausted.
struct AxisData {
    positions: Vec<Vec<[f32; 3]>>,
    velocities: Option<Vec<Vec<[f32; 3]>>>,
    forces: Option<Vec<Vec<[f32; 3]>>>,
    cells: Option<Vec<Vec<f64>>>,
}

/// Timestep numbering and clock inputs shared by every appended frame.
struct Timebase {
    offset: usize,
    first_frame: i64,
    first_seconds: f64,
    time_per_frame: f64,
    scale: f64,
}

impl Timebase {
    /// The wall-clock time of one step, absent when the source records none.
    fn time(&self, step: i64) -> Result<Option<f64>, TngError> {
        optional_time(
            self.first_frame,
            self.first_seconds,
            self.time_per_frame,
            step,
        )
    }
}

fn append_frames(
    frames: &mut Vec<Timestep>,
    steps: &[i64],
    timebase: Timebase,
    axes: AxisData,
) -> Result<(), TngError> {
    let mut positions = axes.positions.into_iter();
    let mut velocities = axes.velocities.map(Vec::into_iter);
    let mut forces = axes.forces.map(Vec::into_iter);
    let mut cells = axes.cells.map(Vec::into_iter);
    for (frame, step) in steps.iter().copied().enumerate() {
        let time = timebase.time(step)?;
        let positions = positions.next().ok_or(TngError::InvalidShape)?;
        let velocities = velocities
            .as_mut()
            .map(|values| values.next().ok_or(TngError::InvalidShape))
            .transpose()?;
        let forces = forces
            .as_mut()
            .map(|values| values.next().ok_or(TngError::InvalidShape))
            .transpose()?;
        let cell = cells
            .as_mut()
            .map(|values| {
                optional_cell(
                    Ok((values.next().ok_or(TngError::InvalidShape)?, 1)),
                    timebase.scale,
                )
            })
            .transpose()?
            .flatten();
        frames.push(Timestep {
            frame: timebase.offset + frame,
            time,
            positions,
            velocities,
            forces,
            cell,
            ..Timestep::default()
        });
    }
    Ok(())
}

fn open(path: &Path) -> Result<Trajectory, TngError> {
    let mut source = Trajectory::new();
    source.input_file_set(path);
    source.file_headers_read(true)?;
    Ok(source)
}

fn position_steps(source: &mut Trajectory, logical_frames: i64) -> Result<Vec<i64>, TngError> {
    let last = logical_frames
        .checked_sub(1)
        .ok_or(TngError::InvalidShape)?;
    let mut steps = Vec::new();
    let mut current = -1;
    loop {
        let found = source.util_trajectory_next_frame_present_data_blocks_find(
            current,
            1,
            &[BlockID::TrajPositions],
        );
        let (next, blocks) = match found {
            Ok(found) => found,
            Err(tng_rs::TngError::Constraint(_) | tng_rs::TngError::NotFound(_))
                if !steps.is_empty() =>
            {
                break;
            }
            Err(error) => return Err(error.into()),
        };
        if next > last && !steps.is_empty() {
            break;
        }
        if blocks != 1 || next <= current || next > last {
            return Err(TngError::InvalidSteps);
        }
        steps.push(next);
        current = next;
        if current == last {
            break;
        }
    }
    if steps.is_empty() {
        Err(TngError::InvalidShape)
    } else {
        Ok(steps)
    }
}

fn required_vectors(values: &[f64], atoms: usize, scale: f64) -> Result<Vec<[f32; 3]>, TngError> {
    let expected = atoms
        .checked_mul(VECTOR_WIDTH)
        .ok_or(TngError::InvalidShape)?;
    if values.len() != expected {
        return Err(TngError::AtomCountMismatch {
            expected: atoms,
            found: values.len() / VECTOR_WIDTH,
        });
    }
    values_to_vectors(values, scale)
}

fn frame_values(
    values: &[f64],
    stride: i64,
    width: usize,
    frame_count: usize,
) -> Result<Vec<Vec<f64>>, TngError> {
    if stride <= 0 {
        return Err(TngError::InvalidShape);
    }
    let expected = width
        .checked_mul(frame_count)
        .ok_or(TngError::InvalidShape)?;
    if values.len() != expected {
        return Err(TngError::InvalidShape);
    }
    Ok(values.chunks_exact(width).map(ToOwned::to_owned).collect())
}

fn frame_vectors(
    values: &[f64],
    stride: i64,
    atoms: usize,
    scale: f64,
    frame_count: usize,
) -> Result<Vec<Vec<[f32; 3]>>, TngError> {
    let width = atoms
        .checked_mul(VECTOR_WIDTH)
        .ok_or(TngError::InvalidShape)?;
    if stride <= 0 {
        return Err(TngError::InvalidShape);
    }
    let expected = width
        .checked_mul(frame_count)
        .ok_or(TngError::InvalidShape)?;
    if values.len() != expected {
        return Err(TngError::InvalidShape);
    }
    values
        .chunks_exact(width)
        .map(|values| required_vectors(values, atoms, scale))
        .collect()
}

fn optional_frame_vectors(
    result: Result<(Vec<f64>, i64), tng_rs::TngError>,
    atoms: usize,
    scale: f64,
    frame_count: usize,
) -> Result<Option<Vec<Vec<[f32; 3]>>>, TngError> {
    match result {
        Ok((values, stride)) => frame_vectors(&values, stride, atoms, scale, frame_count).map(Some),
        Err(tng_rs::TngError::NotFound(_)) => Ok(None),
        Err(error) => Err(error.into()),
    }
}

fn optional_frame_values(
    result: Result<(Vec<f64>, i64), tng_rs::TngError>,
    width: usize,
    frame_count: usize,
) -> Result<Option<Vec<Vec<f64>>>, TngError> {
    match result {
        Ok((values, stride)) => frame_values(&values, stride, width, frame_count).map(Some),
        Err(tng_rs::TngError::NotFound(_)) => Ok(None),
        Err(error) => Err(error.into()),
    }
}

fn values_to_vectors(values: &[f64], scale: f64) -> Result<Vec<[f32; 3]>, TngError> {
    let (chunks, _) = values.as_chunks::<VECTOR_WIDTH>();
    chunks
        .iter()
        .map(|value| {
            f32_triplet([value[0] * scale, value[1] * scale, value[2] * scale])
                .ok_or(TngError::InvalidValue)
        })
        .collect()
}

fn optional_cell(
    result: Result<(Vec<f64>, i64), tng_rs::TngError>,
    scale: f64,
) -> Result<Option<UnitCell>, TngError> {
    let values = match result {
        Ok((values, _)) => values,
        Err(tng_rs::TngError::NotFound(_)) => return Ok(None),
        Err(error) => return Err(error.into()),
    };
    let values: [f64; CELL_WIDTH] = values.try_into().map_err(|_| TngError::InvalidShape)?;
    let vectors = [
        [values[0] * scale, values[1] * scale, values[2] * scale],
        [values[3] * scale, values[4] * scale, values[5] * scale],
        [values[6] * scale, values[7] * scale, values[8] * scale],
    ];
    cell_from_vectors(vectors)
        .map(Some)
        .ok_or(TngError::InvalidValue)
}

fn optional_time(
    first_frame: i64,
    first_seconds: f64,
    time_per_frame: f64,
    step: i64,
) -> Result<Option<f64>, TngError> {
    let explicit = (first_frame == step && first_seconds.is_finite() && first_seconds >= 0.0)
        .then_some(first_seconds * SECONDS_TO_PICOSECONDS);
    if time_per_frame <= 0.0 {
        return Ok(explicit);
    }
    let delta = f64_from_i64(
        step.checked_sub(first_frame)
            .ok_or(TngError::InvalidValue)?,
    )
    .ok_or(TngError::InvalidValue)?;
    let seconds = first_seconds + time_per_frame * delta;
    if seconds.is_finite() {
        Ok(Some(seconds * SECONDS_TO_PICOSECONDS))
    } else {
        Err(TngError::InvalidValue)
    }
}

fn length_to_angstrom(exponent: i64) -> Result<f64, TngError> {
    let power = i32::try_from(exponent.checked_add(10).ok_or(TngError::InvalidValue)?)
        .map_err(|_| TngError::InvalidValue)?;
    let scale = 10.0_f64.powi(power);
    if scale.is_finite() && scale > 0.0 {
        Ok(scale)
    } else {
        Err(TngError::InvalidValue)
    }
}

fn populate_dt(frames: &mut [Timestep]) {
    for index in 1..frames.len() {
        let dt = frames[index]
            .time
            .zip(frames[index - 1].time)
            .map(|(right, left)| right - left);
        if dt.is_some_and(f64::is_finite) {
            frames[index].dt = dt;
        }
    }
}
