//! Deterministic TNG block writer.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;

use tng_rs::data::{Compression, DataType};
use tng_rs::gen_block::BlockID;
use tng_rs::trajectory::{BlockType, Trajectory};

use crate::Timestep;
use crate::cell::vectors_from_cell;
use crate::numeric::{f32_from_f64, f64_from_i64, f64_from_usize};

use super::{TngCompression, TngError, TngWriteOptions};

const PICOSECONDS_TO_SECONDS: f64 = 1.0e-12;
const VALUES_PER_VECTOR: i64 = 3;
const VALUES_PER_CELL: i64 = 9;
const TIME_TOLERANCE: f64 = 1.0e-9;

/// Writes a complete TNG trajectory to `path`.
///
/// The `frame` field of each timestep is the simulation step and must increase
/// strictly. Positions are mandatory; velocities, forces and cells are emitted
/// independently when present. pdbiox canonical units are converted to the
/// explicitly selected file length unit.
///
/// # Errors
///
/// Refuses empty/changing atom sets, mixed auxiliary shapes, absent or negative
/// time, non-increasing steps, invalid cells/units/precision, I/O failures and
/// native-writer invariant failures.
pub fn write_tng(
    path: &Path,
    frames: &[Timestep],
    options: TngWriteOptions,
) -> Result<(), TngError> {
    catch_unwind(AssertUnwindSafe(|| write_inner(path, frames, options)))
        .map_err(|_| TngError::InternalInvariant)?
}

fn write_inner(path: &Path, frames: &[Timestep], options: TngWriteOptions) -> Result<(), TngError> {
    let (atoms, stride, time_per_step) = validate_frames(frames)?;
    let length_scale = length_to_angstrom(options.distance_unit_exponent)?;
    let (compression, precision) = compression(options.compression)?;
    let mut output = Trajectory::new();
    output.output_file_set(path);
    output.distance_unit_exponential = options.distance_unit_exponent;
    output.compression_precision = precision;
    output.set_first_program_name(env!("CARGO_PKG_NAME"));
    output.set_time_per_frame(time_per_step * PICOSECONDS_TO_SECONDS)?;
    add_anonymous_particles(&mut output, atoms);
    output.file_headers_write(options.hashes)?;
    let first = frames.first().ok_or(TngError::InvalidShape)?;
    let last = frames.last().ok_or(TngError::InvalidShape)?;
    let first_step = i64::try_from(first.frame).map_err(|_| TngError::InvalidSteps)?;
    let last_step = i64::try_from(last.frame).map_err(|_| TngError::InvalidSteps)?;
    let span = last_step
        .checked_sub(first_step)
        .and_then(|value| value.checked_add(1))
        .ok_or(TngError::InvalidSteps)?;
    let first_time = first.time.ok_or(TngError::InvalidValue)? * PICOSECONDS_TO_SECONDS;
    output.frame_set_with_time_new(first_step, span, first_time)?;
    let particle_blocks = ParticleBlockContext {
        atoms,
        span,
        stride,
    };

    let positions = collect_vectors(frames, |frame| Some(&frame.positions), 1.0 / length_scale)?;
    particle_blocks.write(
        &mut output,
        BlockID::TrajPositions,
        "POSITIONS",
        &positions,
        compression,
    )?;
    if first.velocities.is_some() {
        let velocities = collect_vectors(
            frames,
            |frame| frame.velocities.as_ref(),
            1.0 / length_scale,
        )?;
        particle_blocks.write(
            &mut output,
            BlockID::TrajVelocities,
            "VELOCITIES",
            &velocities,
            compression,
        )?;
    }
    if first.forces.is_some() {
        let forces = collect_vectors(frames, |frame| frame.forces.as_ref(), length_scale)?;
        particle_blocks.write(
            &mut output,
            BlockID::TrajForces,
            "FORCES",
            &forces,
            Compression::GZip,
        )?;
    }
    if first.cell.is_some() {
        write_cells(&mut output, frames, length_scale, span, stride)?;
    }
    output.frame_set_write(options.hashes)?;
    Ok(())
}

fn validate_frames(frames: &[Timestep]) -> Result<(usize, i64, f64), TngError> {
    let atoms = frames
        .first()
        .map(|frame| frame.positions.len())
        .filter(|count| *count > 0)
        .ok_or(TngError::InvalidShape)?;
    let stride = step_stride(frames)?;
    let reference = frames.first().ok_or(TngError::InvalidShape)?;
    for frame in frames {
        if frame.velocities.is_some() != reference.velocities.is_some()
            || frame.forces.is_some() != reference.forces.is_some()
            || frame.cell.is_some() != reference.cell.is_some()
        {
            return Err(TngError::InvalidShape);
        }
        validate_values(&frame.positions, atoms)?;
        if let Some(values) = &frame.velocities {
            validate_values(values, atoms)?;
        }
        if let Some(values) = &frame.forces {
            validate_values(values, atoms)?;
        }
        if !frame
            .time
            .is_some_and(|time| time.is_finite() && time >= 0.0)
        {
            return Err(TngError::InvalidValue);
        }
    }
    let time_per_step = time_per_step(frames, stride)?;
    Ok((atoms, stride, time_per_step))
}

fn step_stride(frames: &[Timestep]) -> Result<i64, TngError> {
    if frames.len() == 1 {
        return Ok(1);
    }
    let first = frames[1]
        .frame
        .checked_sub(frames[0].frame)
        .and_then(|value| i64::try_from(value).ok())
        .filter(|value| *value > 0)
        .ok_or(TngError::InvalidSteps)?;
    if frames
        .windows(2)
        .all(|pair| pair[1].frame.checked_sub(pair[0].frame) == usize::try_from(first).ok())
    {
        Ok(first)
    } else {
        Err(TngError::InvalidSteps)
    }
}

fn time_per_step(frames: &[Timestep], stride: i64) -> Result<f64, TngError> {
    let rate = if frames.len() == 1 {
        frames[0].dt.ok_or(TngError::InvalidTimeAxis)?
    } else {
        let first = frames[0].time.ok_or(TngError::InvalidTimeAxis)?;
        let second = frames[1].time.ok_or(TngError::InvalidTimeAxis)?;
        let stride = f64_from_i64(stride).ok_or(TngError::InvalidTimeAxis)?;
        (second - first) / stride
    };
    if !rate.is_finite() || rate <= 0.0 {
        return Err(TngError::InvalidTimeAxis);
    }
    let origin_step = f64_from_usize(frames[0].frame).ok_or(TngError::InvalidTimeAxis)?;
    let origin_time = frames[0].time.ok_or(TngError::InvalidTimeAxis)?;
    for frame in frames {
        let time = frame.time.ok_or(TngError::InvalidTimeAxis)?;
        let step = f64_from_usize(frame.frame).ok_or(TngError::InvalidTimeAxis)?;
        let expected = (step - origin_step).mul_add(rate, origin_time);
        let scale = expected.abs().max(time.abs()).max(1.0);
        if (time - expected).abs() > TIME_TOLERANCE * scale {
            return Err(TngError::InvalidTimeAxis);
        }
    }
    Ok(rate)
}

fn validate_values(values: &[[f32; 3]], atoms: usize) -> Result<(), TngError> {
    if values.len() != atoms {
        return Err(TngError::AtomCountMismatch {
            expected: atoms,
            found: values.len(),
        });
    }
    if values.iter().flatten().all(|value| value.is_finite()) {
        Ok(())
    } else {
        Err(TngError::InvalidValue)
    }
}

fn add_anonymous_particles(output: &mut Trajectory, atoms: usize) {
    let molecule = output.add_molecule("");
    let chain = output.add_chain(molecule, "");
    let residue = output.chain_residue_add(molecule, chain, "");
    for _ in 0..atoms {
        output.residue_atom_add(molecule, residue, "", "");
    }
    output.molecule_cnt_set(molecule, 1);
}

#[derive(Clone, Copy)]
struct ParticleBlockContext {
    atoms: usize,
    span: i64,
    stride: i64,
}

impl ParticleBlockContext {
    fn write(
        self,
        output: &mut Trajectory,
        id: BlockID,
        name: &str,
        values: &[f32],
        compression: Compression,
    ) -> Result<(), TngError> {
        if values.iter().any(|value| !value.is_finite()) {
            return Err(TngError::InvalidValue);
        }
        let bytes = native_bytes(values);
        output.particle_data_block_add(
            id,
            name,
            DataType::Float,
            &BlockType::Trajectory,
            self.span,
            VALUES_PER_VECTOR,
            self.stride,
            0,
            i64::try_from(self.atoms).map_err(|_| TngError::InvalidShape)?,
            compression,
            Some(&bytes),
        )?;
        Ok(())
    }
}

fn collect_vectors<'a>(
    frames: &'a [Timestep],
    select: impl Fn(&'a Timestep) -> Option<&'a Vec<[f32; 3]>>,
    scale: f64,
) -> Result<Vec<f32>, TngError> {
    let mut output = Vec::new();
    for frame in frames {
        let values = select(frame).ok_or(TngError::InvalidShape)?;
        let converted = values
            .iter()
            .flatten()
            .map(|value| f32_from_f64(f64::from(*value) * scale).ok_or(TngError::InvalidValue))
            .collect::<Result<Vec<_>, _>>()?;
        output.extend(converted);
    }
    if output.iter().all(|value| value.is_finite()) {
        Ok(output)
    } else {
        Err(TngError::InvalidValue)
    }
}

fn write_cells(
    output: &mut Trajectory,
    frames: &[Timestep],
    length_scale: f64,
    span: i64,
    stride: i64,
) -> Result<(), TngError> {
    let mut values = Vec::with_capacity(
        frames.len() * usize::try_from(VALUES_PER_CELL).map_err(|_| TngError::InvalidShape)?,
    );
    for frame in frames {
        let cell = frame.cell.ok_or(TngError::InvalidShape)?;
        let vectors = vectors_from_cell(cell).ok_or(TngError::InvalidValue)?;
        let converted = vectors
            .into_iter()
            .flatten()
            .map(|value| f32_from_f64(value / length_scale).ok_or(TngError::InvalidValue))
            .collect::<Result<Vec<_>, _>>()?;
        values.extend(converted);
    }
    let bytes = native_bytes(&values);
    output.data_block_add(
        BlockID::TrajBoxShape,
        "BOX SHAPE",
        DataType::Float,
        &BlockType::Trajectory,
        span,
        VALUES_PER_CELL,
        stride,
        Compression::GZip,
        Some(&bytes),
    )?;
    Ok(())
}

fn native_bytes(values: &[f32]) -> Vec<u8> {
    values
        .iter()
        .flat_map(|value| value.to_ne_bytes())
        .collect()
}

fn compression(value: TngCompression) -> Result<(Compression, f64), TngError> {
    match value {
        TngCompression::Uncompressed => Ok((Compression::Uncompressed, 1.0)),
        TngCompression::Lossless => Ok((Compression::GZip, 1_000.0)),
        TngCompression::Lossy { precision } if precision.is_finite() && precision > 0.0 => {
            Ok((Compression::TNG, precision))
        }
        TngCompression::Lossy { .. } => Err(TngError::InvalidValue),
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
