//! Deterministic CHARMM-compatible DCD writing.

use crate::Timestep;
use crate::dcd::{DcdEndian, DcdError};

/// Controls stored in a newly written DCD header.
#[derive(Clone, Debug, PartialEq)]
pub struct DcdWriteOptions {
    /// File byte order.
    pub endian: DcdEndian,
    /// One title record, limited to the fixed 80-byte field.
    pub title: Box<str>,
    /// Integrator step of the first frame.
    pub start_step: i32,
    /// Integrator steps between frames.
    pub save_interval: i32,
    /// Timestep in native AKMA units.
    pub delta_akma: f32,
}

impl Default for DcdWriteOptions {
    fn default() -> Self {
        Self {
            endian: DcdEndian::Little,
            title: "molframe trajectory".into(),
            start_step: 0,
            save_interval: 1,
            delta_akma: 1.0,
        }
    }
}

/// Writes a CHARMM-style DCD with optional per-frame periodic cells.
///
/// # Errors
///
/// Returns a controls, count or cell error before publishing partial bytes.
pub fn write_dcd(frames: &[Timestep], options: &DcdWriteOptions) -> Result<Vec<u8>, DcdError> {
    if frames.is_empty()
        || options.save_interval <= 0
        || !options.delta_akma.is_finite()
        || options.delta_akma <= 0.0
        || options.title.len() > 80
        || options.title.contains(['\n', '\r'])
    {
        return Err(DcdError::InvalidControls);
    }
    let atom_count = frames[0].positions.len();
    let has_cell = frames[0].cell.is_some();
    if atom_count == 0
        || frames
            .iter()
            .any(|frame| frame.positions.len() != atom_count || frame.cell.is_some() != has_cell)
    {
        return Err(DcdError::CoordinateCount);
    }
    let frame_count = i32::try_from(frames.len()).map_err(|_| DcdError::InvalidControls)?;
    let atom_count_i32 = i32::try_from(atom_count).map_err(|_| DcdError::InvalidControls)?;
    let mut output = Vec::new();
    let mut header = vec![0u8; 84];
    header[0..4].copy_from_slice(b"CORD");
    put_i32(&mut header, 4, frame_count, options.endian);
    put_i32(&mut header, 8, options.start_step, options.endian);
    put_i32(&mut header, 12, options.save_interval, options.endian);
    let last_step = options
        .start_step
        .checked_add(
            (frame_count - 1)
                .checked_mul(options.save_interval)
                .ok_or(DcdError::InvalidControls)?,
        )
        .ok_or(DcdError::InvalidControls)?;
    put_i32(&mut header, 16, last_step, options.endian);
    put_f32(&mut header, 40, options.delta_akma, options.endian);
    put_i32(&mut header, 44, i32::from(has_cell), options.endian);
    put_i32(&mut header, 80, 24, options.endian);
    record(&mut output, &header, options.endian)?;
    let mut title = vec![0u8; 84];
    put_i32(&mut title, 0, 1, options.endian);
    title[4..].fill(b' ');
    let bytes = options.title.as_bytes();
    title[4..4 + bytes.len()].copy_from_slice(bytes);
    record(&mut output, &title, options.endian)?;
    record(
        &mut output,
        &integer_bytes(atom_count_i32, options.endian),
        options.endian,
    )?;
    for frame in frames {
        if let Some(cell) = frame.cell {
            let values = [
                cell.lengths[0],
                cell.angles[2].to_radians().cos(),
                cell.lengths[1],
                cell.angles[1].to_radians().cos(),
                cell.angles[0].to_radians().cos(),
                cell.lengths[2],
            ];
            if values.iter().any(|value| !value.is_finite())
                || cell.lengths.iter().any(|length| *length <= 0.0)
            {
                return Err(DcdError::InvalidCell);
            }
            let mut bytes = Vec::with_capacity(48);
            for value in values {
                bytes.extend(float64_bytes(value, options.endian));
            }
            record(&mut output, &bytes, options.endian)?;
        }
        for axis in 0..3 {
            let mut bytes = Vec::with_capacity(atom_count * 4);
            for position in &frame.positions {
                bytes.extend(float32_bytes(position[axis], options.endian));
            }
            record(&mut output, &bytes, options.endian)?;
        }
    }
    Ok(output)
}

fn record(output: &mut Vec<u8>, payload: &[u8], endian: DcdEndian) -> Result<(), DcdError> {
    let length = i32::try_from(payload.len()).map_err(|_| DcdError::InvalidControls)?;
    output.extend(integer_bytes(length, endian));
    output.extend(payload);
    output.extend(integer_bytes(length, endian));
    Ok(())
}

fn put_i32(bytes: &mut [u8], offset: usize, value: i32, endian: DcdEndian) {
    bytes[offset..offset + 4].copy_from_slice(&integer_bytes(value, endian));
}

fn put_f32(bytes: &mut [u8], offset: usize, value: f32, endian: DcdEndian) {
    bytes[offset..offset + 4].copy_from_slice(&float32_bytes(value, endian));
}

fn integer_bytes(value: i32, endian: DcdEndian) -> [u8; 4] {
    match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
    }
}

fn float32_bytes(value: f32, endian: DcdEndian) -> [u8; 4] {
    match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
    }
}

fn float64_bytes(value: f64, endian: DcdEndian) -> [u8; 8] {
    match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
    }
}
