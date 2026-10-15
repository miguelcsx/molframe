//! Deterministic writer for IBIsCO/YASP TRZ trajectories.

use crate::cell::vectors_from_cell;
use crate::{FrameValue, Timestep, TrzError, TrzTrajectory};

/// Writes a complete little-endian TRZ trajectory without synthesizing metadata.
///
/// Every frame must provide velocity data and the seven metadata keys emitted by
/// [`crate::parse_trz`]. Force presence must match `trajectory.has_forces`.
///
/// # Errors
///
/// Refuses absent metadata, inconsistent atom counts, invalid cells, non-finite
/// values and titles longer than the fixed 80-byte field.
pub fn write_trz(trajectory: &TrzTrajectory) -> Result<Vec<u8>, TrzError> {
    if trajectory.title.len() > 80 || trajectory.frames.is_empty() {
        return Err(TrzError::InvalidHeader);
    }
    let atoms = trajectory.frames[0].positions.len();
    if atoms == 0 {
        return Err(invalid(0));
    }
    let mut output = Vec::new();
    put_i32(&mut output, 80);
    output.extend(trajectory.title.as_bytes());
    output.resize(output.len() + 80 - trajectory.title.len(), b' ');
    put_i32(&mut output, 80);
    put_i32(&mut output, 4);
    put_i32(&mut output, if trajectory.has_forces { 20 } else { 10 });
    put_i32(&mut output, 4);
    for (ordinal, frame) in trajectory.frames.iter().enumerate() {
        write_frame(&mut output, frame, ordinal, atoms, trajectory.has_forces)?;
    }
    Ok(output)
}

fn write_frame(
    output: &mut Vec<u8>,
    frame: &Timestep,
    ordinal: usize,
    atoms: usize,
    has_forces: bool,
) -> Result<(), TrzError> {
    validate(frame, atoms, has_forces)?;
    let frame_number = i32::try_from(ordinal + 1).map_err(|_| invalid(output.len()))?;
    let atom_count = i32::try_from(atoms).map_err(|_| invalid(output.len()))?;
    put_i32(output, 20);
    put_i32(output, frame_number);
    put_i32(output, integer(frame, "trajectory_step")?);
    put_i32(output, atom_count);
    put_f64(
        output,
        frame
            .time
            .ok_or(TrzError::MissingMetadata { name: "time" })?,
    );
    put_i32(output, 20);
    put_i32(output, 72);
    let vectors = match frame.cell {
        Some(cell) => vectors_from_cell(cell).ok_or_else(|| invalid(output.len()))?,
        None => [[0.0; 3]; 3],
    };
    for vector in vectors {
        for value in vector {
            put_f64(output, value / 10.0);
        }
    }
    put_i32(output, 72);
    put_i32(output, 56);
    put_f64(output, scalar(frame, "pressure")?);
    for value in vector(frame, "pressure_tensor", 6)? {
        put_f64(output, value);
    }
    put_i32(output, 56);
    put_i32(output, 52);
    put_i32(output, 6);
    for name in [
        "total_energy",
        "potential_energy",
        "kinetic_energy",
        "temperature",
    ] {
        put_f64(output, scalar(frame, name)?);
    }
    put_f64(output, 0.0);
    put_f64(output, 0.0);
    put_i32(output, 52);
    write_vectors(output, &frame.positions, 0.1)?;
    write_vectors(
        output,
        frame
            .velocities
            .as_deref()
            .ok_or(TrzError::MissingMetadata { name: "velocities" })?,
        0.1,
    )?;
    if has_forces {
        write_vectors(
            output,
            frame
                .forces
                .as_deref()
                .ok_or(TrzError::MissingMetadata { name: "forces" })?,
            10.0,
        )?;
    }
    Ok(())
}

fn validate(frame: &Timestep, atoms: usize, has_forces: bool) -> Result<(), TrzError> {
    if frame.positions.len() != atoms
        || frame
            .velocities
            .as_ref()
            .is_some_and(|values| values.len() != atoms)
        || frame
            .forces
            .as_ref()
            .is_some_and(|values| values.len() != atoms)
        || frame.forces.is_some() != has_forces
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
    {
        return Err(invalid(0));
    }
    Ok(())
}

fn write_vectors(output: &mut Vec<u8>, vectors: &[[f32; 3]], scale: f32) -> Result<(), TrzError> {
    let marker = i32::try_from(vectors.len() * 4).map_err(|_| invalid(output.len()))?;
    for axis in 0..3 {
        put_i32(output, marker);
        for vector in vectors {
            output.extend((vector[axis] * scale).to_le_bytes());
        }
        put_i32(output, marker);
    }
    Ok(())
}

fn scalar(frame: &Timestep, name: &'static str) -> Result<f64, TrzError> {
    match frame.data.get(name) {
        Some(FrameValue::Float(value)) if value.is_finite() => Ok(*value),
        _ => Err(TrzError::MissingMetadata { name }),
    }
}

fn integer(frame: &Timestep, name: &'static str) -> Result<i32, TrzError> {
    match frame.data.get(name) {
        Some(FrameValue::Integer(value)) => i32::try_from(*value).map_err(|_| invalid(0)),
        _ => Err(TrzError::MissingMetadata { name }),
    }
}

fn vector(frame: &Timestep, name: &'static str, length: usize) -> Result<Vec<f64>, TrzError> {
    match frame.data.get(name) {
        Some(FrameValue::Floats(values))
            if values.len() == length && values.iter().all(|value| value.is_finite()) =>
        {
            Ok(values.clone())
        }
        _ => Err(TrzError::MissingMetadata { name }),
    }
}

fn invalid(offset: usize) -> TrzError {
    TrzError::InvalidFrame { offset }
}
fn put_i32(output: &mut Vec<u8>, value: i32) {
    output.extend(value.to_le_bytes());
}
fn put_f64(output: &mut Vec<u8>, value: f64) {
    output.extend(value.to_le_bytes());
}
