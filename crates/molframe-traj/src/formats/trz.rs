//! Little-endian IBIsCO/YASP TRZ trajectory reader.

use crate::cell::cell_from_vectors;
use crate::{FrameValue, Timestep};
use std::collections::BTreeMap;

const HEADER_BYTES: usize = 100;
const NM_TO_ANGSTROM: f32 = 10.0;
const FORCE_TO_CANONICAL: f32 = 0.1;

/// TRZ header plus all coordinate frames.
#[derive(Clone, Debug, PartialEq)]
pub struct TrzTrajectory {
    /// Fixed-width file title without trailing spaces.
    pub title: Box<str>,
    /// Whether each frame carries force arrays.
    pub has_forces: bool,
    /// Frames in stream order.
    pub frames: Vec<Timestep>,
}

/// Malformed or incomplete TRZ data.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum TrzError {
    /// Header markers or force controls are invalid.
    #[error("invalid TRZ header")]
    InvalidHeader,
    /// Frame markers, atom counts or numeric values are invalid.
    #[error("invalid TRZ frame at byte {offset}")]
    InvalidFrame {
        /// Frame start or first unavailable byte.
        offset: usize,
    },
    /// A writer input omits a required TRZ metadata field.
    #[error("TRZ frame is missing required metadata: {name}")]
    MissingMetadata {
        /// Required `Timestep::data` key.
        name: &'static str,
    },
}

/// Parses a complete TRZ stream, including velocities, optional forces,
/// pressure tensors and energy metadata.
///
/// # Errors
///
/// Returns the first invalid header, marker, count or numeric error.
pub fn parse_trz(bytes: &[u8]) -> Result<TrzTrajectory, TrzError> {
    let mut input = Input::new(bytes);
    input.marker(80)?;
    let title = String::from_utf8_lossy(input.take(80)?).trim_end().into();
    input.marker(80)?;
    input.marker(4)?;
    let has_forces = match input.i32()? {
        10 => false,
        20 => true,
        _ => return Err(TrzError::InvalidHeader),
    };
    input.marker(4)?;
    if input.cursor != HEADER_BYTES {
        return Err(TrzError::InvalidHeader);
    }
    let mut frames = Vec::new();
    while !input.is_empty() {
        frames.push(read_frame(&mut input, frames.len(), has_forces)?);
    }
    populate_dt(&mut frames);
    Ok(TrzTrajectory {
        title,
        has_forces,
        frames,
    })
}

fn read_frame(
    input: &mut Input<'_>,
    ordinal: usize,
    has_forces: bool,
) -> Result<Timestep, TrzError> {
    let start = input.cursor;
    input.marker(20)?;
    let number = input.i32()?;
    let trajectory_step = input.i32()?;
    let atoms = usize::try_from(input.i32()?).map_err(|_| invalid(start))?;
    let time = input.f64()?;
    input.marker(20)?;
    if atoms == 0 || number <= 0 || !time.is_finite() {
        return Err(invalid(start));
    }
    let cell = read_cell(input, start)?;
    input.marker(56)?;
    let pressure = input.f64()?;
    let pressure_tensor = read_f64_array::<6>(input)?;
    input.marker(56)?;
    let energy_marker = input.marker_one_of([52, 60])?;
    if input.i32()? != 6 {
        return Err(invalid(start));
    }
    let energies = read_f64_array::<4>(input)?;
    let blanks = read_f64_array::<2>(input)?;
    input.marker(energy_marker)?;
    if [pressure, time]
        .into_iter()
        .chain(pressure_tensor)
        .chain(energies)
        .chain(blanks)
        .any(|value| !value.is_finite())
    {
        return Err(invalid(start));
    }
    let x = read_axis(input, atoms)?;
    let y = read_axis(input, atoms)?;
    let z = read_axis(input, atoms)?;
    let vx = read_axis(input, atoms)?;
    let vy = read_axis(input, atoms)?;
    let vz = read_axis(input, atoms)?;
    let positions = join_axes(&x, &y, &z, NM_TO_ANGSTROM);
    let velocities = Some(join_axes(&vx, &vy, &vz, NM_TO_ANGSTROM));
    let forces = if has_forces {
        let fx = read_axis(input, atoms)?;
        let fy = read_axis(input, atoms)?;
        let fz = read_axis(input, atoms)?;
        Some(join_axes(&fx, &fy, &fz, FORCE_TO_CANONICAL))
    } else {
        None
    };
    let mut data = BTreeMap::new();
    data.insert(
        "trajectory_step".into(),
        FrameValue::Integer(i64::from(trajectory_step)),
    );
    data.insert("pressure".into(), FrameValue::Float(pressure));
    data.insert(
        "pressure_tensor".into(),
        FrameValue::Floats(pressure_tensor.into()),
    );
    for (name, value) in [
        ("total_energy", energies[0]),
        ("potential_energy", energies[1]),
        ("kinetic_energy", energies[2]),
        ("temperature", energies[3]),
    ] {
        data.insert(name.into(), FrameValue::Float(value));
    }
    Ok(Timestep {
        frame: ordinal,
        time: Some(time),
        positions,
        velocities,
        forces,
        cell,
        data,
        ..Timestep::default()
    })
}

fn read_cell(
    input: &mut Input<'_>,
    start: usize,
) -> Result<Option<molframe_core::structure::UnitCell>, TrzError> {
    input.marker(72)?;
    let raw = read_f64_array::<9>(input)?;
    input.marker(72)?;
    if raw.iter().all(|value| value.abs() <= f64::EPSILON) {
        return Ok(None);
    }
    let vectors = std::array::from_fn(|row| {
        std::array::from_fn(|column| raw[row * 3 + column] * f64::from(NM_TO_ANGSTROM))
    });
    cell_from_vectors(vectors)
        .ok_or_else(|| invalid(start))
        .map(Some)
}

fn read_axis(input: &mut Input<'_>, atoms: usize) -> Result<Vec<f32>, TrzError> {
    let bytes = i32::try_from(atoms.checked_mul(4).ok_or_else(|| invalid(input.cursor))?)
        .map_err(|_| invalid(input.cursor))?;
    input.marker(bytes)?;
    let values = (0..atoms)
        .map(|_| input.f32())
        .collect::<Result<Vec<_>, _>>()?;
    input.marker(bytes)?;
    if values.iter().any(|value| !value.is_finite()) {
        return Err(invalid(input.cursor));
    }
    Ok(values)
}

fn join_axes(x: &[f32], y: &[f32], z: &[f32], scale: f32) -> Vec<[f32; 3]> {
    (0..x.len())
        .map(|index| [x[index] * scale, y[index] * scale, z[index] * scale])
        .collect()
}

fn read_f64_array<const N: usize>(input: &mut Input<'_>) -> Result<[f64; N], TrzError> {
    let mut values = [0.0; N];
    for value in &mut values {
        *value = input.f64()?;
    }
    Ok(values)
}

fn populate_dt(frames: &mut [Timestep]) {
    for index in 1..frames.len() {
        if let (Some(previous), Some(current)) = (frames[index - 1].time, frames[index].time) {
            frames[index].dt = Some(current - previous);
        }
    }
    if frames.len() > 1 {
        frames[0].dt = frames[1].dt;
    }
}

fn invalid(offset: usize) -> TrzError {
    TrzError::InvalidFrame { offset }
}

struct Input<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Input<'a> {
    const fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }
    fn is_empty(&self) -> bool {
        self.cursor == self.bytes.len()
    }
    fn take(&mut self, length: usize) -> Result<&'a [u8], TrzError> {
        let start = self.cursor;
        let end = start.checked_add(length).ok_or_else(|| invalid(start))?;
        let value = self.bytes.get(start..end).ok_or_else(|| invalid(start))?;
        self.cursor = end;
        Ok(value)
    }
    fn i32(&mut self) -> Result<i32, TrzError> {
        let offset = self.cursor;
        let raw = self.take(4)?.try_into().map_err(|_| invalid(offset))?;
        Ok(i32::from_le_bytes(raw))
    }
    fn f32(&mut self) -> Result<f32, TrzError> {
        let offset = self.cursor;
        let raw = self.take(4)?.try_into().map_err(|_| invalid(offset))?;
        Ok(f32::from_le_bytes(raw))
    }
    fn f64(&mut self) -> Result<f64, TrzError> {
        let offset = self.cursor;
        let raw = self.take(8)?.try_into().map_err(|_| invalid(offset))?;
        Ok(f64::from_le_bytes(raw))
    }
    fn marker(&mut self, expected: i32) -> Result<(), TrzError> {
        let offset = self.cursor;
        if self.i32()? == expected {
            Ok(())
        } else {
            Err(invalid(offset))
        }
    }
    fn marker_one_of(&mut self, expected: [i32; 2]) -> Result<i32, TrzError> {
        let offset = self.cursor;
        let marker = self.i32()?;
        if expected.contains(&marker) {
            Ok(marker)
        } else {
            Err(invalid(offset))
        }
    }
}

#[cfg(test)]
#[path = "trz_tests.rs"]
mod tests;
