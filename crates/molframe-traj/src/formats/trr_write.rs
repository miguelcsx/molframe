//! Deterministic XDR writer for GROMACS TRR trajectories.

use crate::cell::vectors_from_cell;
use crate::numeric::{f32_from_f64, f64_from_usize};
use crate::{Timestep, TrrError, TrrPrecision};

const MAGIC: i32 = 1993;
const VERSION: &str = "GMX_trn_file";

/// Controls the numeric representation emitted by [`write_trr`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TrrWriteOptions {
    /// Precision used for all real-valued blocks.
    pub precision: TrrPrecision,
}

impl Default for TrrWriteOptions {
    fn default() -> Self {
        Self {
            precision: TrrPrecision::Single,
        }
    }
}

/// Writes concatenated TRR frames with optional box, velocity and force blocks.
///
/// # Errors
///
/// Refuses empty trajectories, changing atom counts, mismatched auxiliary
/// arrays, invalid unit cells and values that exceed TRR integer fields.
pub fn write_trr(frames: &[Timestep], options: TrrWriteOptions) -> Result<Vec<u8>, TrrError> {
    write_trr_with_precisions(frames, &vec![options.precision; frames.len()])
}

/// Writes concatenated TRR frames while preserving each frame's representation.
///
/// # Errors
///
/// Returns the same validation errors as [`write_trr`], and rejects a precision
/// count that differs from the number of frames.
pub fn write_trr_with_precisions(
    frames: &[Timestep],
    precisions: &[TrrPrecision],
) -> Result<Vec<u8>, TrrError> {
    let atoms = frames
        .first()
        .map(|frame| frame.positions.len())
        .ok_or(TrrError::InvalidSizes)?;
    let atoms_i32 = i32::try_from(atoms).map_err(|_| TrrError::InvalidSizes)?;
    if atoms == 0 {
        return Err(TrrError::InvalidSizes);
    }
    if precisions.len() != frames.len() {
        return Err(TrrError::InvalidSizes);
    }
    let mut output = Vec::new();
    for (index, (frame, precision)) in frames.iter().zip(precisions).enumerate() {
        validate_frame(frame, atoms)?;
        let box_size = usize::from(frame.cell.is_some()) * 9 * precision.bytes();
        let positions_size = atoms * 3 * precision.bytes();
        let velocities_size = usize::from(frame.velocities.is_some()) * positions_size;
        let forces_size = usize::from(frame.forces.is_some()) * positions_size;
        put_i32(&mut output, MAGIC);
        put_i32(
            &mut output,
            i32::try_from(VERSION.len() + 1).map_err(|_| TrrError::InvalidSizes)?,
        );
        put_string(&mut output, VERSION)?;
        for size in [
            0,
            0,
            box_size,
            0,
            0,
            0,
            0,
            positions_size,
            velocities_size,
            forces_size,
        ] {
            put_i32(
                &mut output,
                i32::try_from(size).map_err(|_| TrrError::InvalidSizes)?,
            );
        }
        put_i32(&mut output, atoms_i32);
        put_i32(
            &mut output,
            i32::try_from(frame.frame).map_err(|_| TrrError::InvalidSizes)?,
        );
        put_i32(&mut output, 0);
        let time = match frame.time {
            Some(time) => time,
            None => f64_from_usize(index).ok_or(TrrError::InvalidSizes)?,
        };
        put_real(&mut output, time, *precision)?;
        put_real(&mut output, 0.0, *precision)?;
        if let Some(cell) = frame.cell {
            let vectors = vectors_from_cell(cell).ok_or(TrrError::InvalidValue)?;
            for vector in vectors {
                for value in vector {
                    put_real(&mut output, value / 10.0, *precision)?;
                }
            }
        }
        put_vectors(&mut output, &frame.positions, 0.1, *precision)?;
        if let Some(velocities) = &frame.velocities {
            put_vectors(&mut output, velocities, 0.1, *precision)?;
        }
        if let Some(forces) = &frame.forces {
            put_vectors(&mut output, forces, 10.0, *precision)?;
        }
    }
    Ok(output)
}

fn validate_frame(frame: &Timestep, atoms: usize) -> Result<(), TrrError> {
    if frame.positions.len() != atoms
        || frame
            .velocities
            .as_ref()
            .is_some_and(|values| values.len() != atoms)
        || frame
            .forces
            .as_ref()
            .is_some_and(|values| values.len() != atoms)
        || frame
            .positions
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
        || frame
            .velocities
            .as_ref()
            .into_iter()
            .flatten()
            .flatten()
            .any(|value| !value.is_finite())
        || frame
            .forces
            .as_ref()
            .into_iter()
            .flatten()
            .flatten()
            .any(|value| !value.is_finite())
        || frame.time.is_some_and(|time| !time.is_finite())
    {
        return Err(TrrError::InvalidValue);
    }
    Ok(())
}

fn put_vectors(
    output: &mut Vec<u8>,
    vectors: &[[f32; 3]],
    scale: f64,
    precision: TrrPrecision,
) -> Result<(), TrrError> {
    for vector in vectors {
        for value in vector {
            put_real(output, f64::from(*value) * scale, precision)?;
        }
    }
    Ok(())
}

fn put_real(output: &mut Vec<u8>, value: f64, precision: TrrPrecision) -> Result<(), TrrError> {
    match precision {
        TrrPrecision::Single => output.extend(
            f32_from_f64(value)
                .ok_or(TrrError::InvalidValue)?
                .to_be_bytes(),
        ),
        TrrPrecision::Double => output.extend(value.to_be_bytes()),
    }
    Ok(())
}

fn put_i32(output: &mut Vec<u8>, value: i32) {
    output.extend(value.to_be_bytes());
}

fn put_string(output: &mut Vec<u8>, value: &str) -> Result<(), TrrError> {
    put_i32(
        output,
        i32::try_from(value.len()).map_err(|_| TrrError::InvalidSizes)?,
    );
    output.extend(value.as_bytes());
    output.resize(output.len() + (4 - value.len() % 4) % 4, 0);
    Ok(())
}
