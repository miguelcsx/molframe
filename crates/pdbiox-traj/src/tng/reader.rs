//! Native TNG reader with canonical-unit conversion.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;

use pdbiox_core::structure::UnitCell;
use tng_rs::data::Compression;
use tng_rs::gen_block::BlockID;
use tng_rs::trajectory::Trajectory;

use crate::Timestep;
use crate::cell::cell_from_vectors;
use crate::numeric::f32_triplet;

use super::{TngCompression, TngError, TngTrajectory};

const SECONDS_TO_PICOSECONDS: f64 = 1.0e12;
const VECTOR_WIDTH: usize = 3;
const CELL_WIDTH: usize = 9;

/// Reads every coordinate-bearing frame from a TNG file.
///
/// Sparse simulation steps are preserved. Positions, velocities, forces,
/// periodic cells and time are converted to pdbiox canonical units. The work is
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
    let expected = steps.len();
    let exponent = source.distance_unit_exponential_get();
    let length_scale = length_to_angstrom(exponent)?;
    let precision = source.compression_precision_get();
    if !precision.is_finite() || precision <= 0.0 {
        return Err(TngError::InvalidValue);
    }
    source.frame_set_of_frame_find(steps[0])?;
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

    let frame_width = atoms
        .checked_mul(VECTOR_WIDTH)
        .ok_or(TngError::InvalidShape)?;
    let bulk_value_count = expected
        .checked_mul(frame_width)
        .ok_or(TngError::InvalidShape)?;
    let bulk_frames = match source.util_pos_read() {
        Ok((values, _)) if values.len() == bulk_value_count => Some(
            values
                .chunks_exact(frame_width)
                .map(|frame| f32_values_to_vectors(frame, length_scale))
                .collect::<Result<Vec<_>, _>>()?,
        ),
        Ok(_) => None,
        Err(_) => {
            source = open(path)?;
            None
        }
    };
    let mut bulk_frames = bulk_frames.map(Vec::into_iter);
    let mut frames = Vec::with_capacity(steps.len());
    for (frame, step) in steps.iter().copied().enumerate() {
        source.frame_set_of_frame_find(step)?;
        let positions = if let Some(values) = &mut bulk_frames {
            values.next().ok_or(TngError::InvalidShape)?
        } else {
            let (values, _) = source.util_pos_read_range(step, step)?;
            required_vectors(&values, atoms, length_scale)?
        };
        let velocities =
            optional_vectors(source.util_vel_read_range(step, step), atoms, length_scale)?;
        let forces = optional_vectors(
            source.util_force_read_range(step, step),
            atoms,
            1.0 / length_scale,
        )?;
        let cell = optional_cell(source.util_box_shape_read_range(step, step), length_scale)?;
        let time = optional_time(&mut source, step)?;
        frames.push(Timestep {
            frame,
            time,
            positions,
            velocities,
            forces,
            cell,
            ..Timestep::default()
        });
    }
    populate_dt(&mut frames);
    Ok(TngTrajectory {
        frames,
        steps,
        distance_unit_exponent: exponent,
        compression_precision: precision,
        compression,
    })
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

fn f32_values_to_vectors(values: &[f32], scale: f64) -> Result<Vec<[f32; 3]>, TngError> {
    let mut chunks = values.chunks_exact(VECTOR_WIDTH);
    let vectors = chunks
        .by_ref()
        .map(|value| {
            f32_triplet([
                f64::from(value[0]) * scale,
                f64::from(value[1]) * scale,
                f64::from(value[2]) * scale,
            ])
            .ok_or(TngError::InvalidValue)
        })
        .collect::<Result<Vec<_>, _>>()?;
    if chunks.remainder().is_empty() {
        Ok(vectors)
    } else {
        Err(TngError::InvalidShape)
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

fn optional_vectors(
    result: Result<(Vec<f64>, i64), tng_rs::TngError>,
    atoms: usize,
    scale: f64,
) -> Result<Option<Vec<[f32; 3]>>, TngError> {
    match result {
        Ok((values, _)) => required_vectors(&values, atoms, scale).map(Some),
        Err(tng_rs::TngError::NotFound(_)) => Ok(None),
        Err(error) => Err(error.into()),
    }
}

fn values_to_vectors(values: &[f64], scale: f64) -> Result<Vec<[f32; 3]>, TngError> {
    values
        .chunks_exact(VECTOR_WIDTH)
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

fn optional_time(source: &mut Trajectory, step: i64) -> Result<Option<f64>, TngError> {
    let frame_set = &source.current_trajectory_frame_set;
    let explicit = (frame_set.first_frame == step
        && frame_set.first_frame_time.is_finite()
        && frame_set.first_frame_time >= 0.0)
        .then_some(frame_set.first_frame_time * SECONDS_TO_PICOSECONDS);
    match source.util_time_of_frame_get(step) {
        Ok(seconds) if seconds.is_finite() => Ok(Some(seconds * SECONDS_TO_PICOSECONDS)),
        Ok(_) => Err(TngError::InvalidValue),
        Err(tng_rs::TngError::Constraint(_) | tng_rs::TngError::NotFound(_)) => Ok(explicit),
        Err(error) => Err(error.into()),
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
