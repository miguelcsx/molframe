//! Portable DCD reader for CHARMM, NAMD, X-PLOR and LAMMPS output.

use crate::Timestep;
use pdbiox_core::structure::UnitCell;
use std::collections::BTreeSet;

#[path = "dcd/reader.rs"]
mod reader;

pub use reader::DcdReader;

pub(crate) const AKMA_TO_PS: f64 = 0.048_888_21;

/// Integer and floating-point byte order used by a DCD stream.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DcdEndian {
    /// Least-significant byte first.
    Little,
    /// Most-significant byte first.
    Big,
}

/// Parsed DCD header controls.
#[derive(Clone, Debug, PartialEq)]
pub struct DcdHeader {
    /// Declared number of coordinate sets.
    pub frame_count: usize,
    /// Number of atoms in a complete frame.
    pub atom_count: usize,
    /// Integrator step of the first frame.
    pub start_step: i32,
    /// Integrator steps between saved frames.
    pub save_interval: i32,
    /// Native DCD timestep in AKMA units.
    pub delta_akma: f64,
    /// Atoms held fixed after the first frame.
    pub fixed_atom_count: usize,
    /// Title records without fixed-width padding.
    pub titles: Vec<Box<str>>,
    /// File byte order.
    pub endian: DcdEndian,
    /// Whether the CHARMM version marker is present.
    pub charmm: bool,
    /// Whether each frame declares an extra unit-cell block.
    pub has_unit_cell: bool,
    /// Whether a fourth coordinate block follows Z and is ignored.
    pub has_fourth_dimension: bool,
}

/// A complete DCD trajectory in canonical pdbiox units.
#[derive(Clone, Debug, PartialEq)]
pub struct DcdTrajectory {
    /// Format controls and titles.
    pub header: DcdHeader,
    /// Coordinate frames.
    pub frames: Vec<Timestep>,
}

/// Malformed or unsupported DCD bytes.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum DcdError {
    /// The file could not be opened or read.
    #[error("DCD I/O failed: {kind:?}")]
    Io {
        /// Portable operating-system error category.
        kind: std::io::ErrorKind,
    },
    /// Initial record marker is neither endian representation of 84.
    #[error("DCD header marker is invalid")]
    InvalidHeader,
    /// A Fortran record is truncated or has unequal markers.
    #[error("invalid DCD Fortran record at byte {offset}")]
    InvalidRecord {
        /// Byte offset of the leading marker.
        offset: usize,
    },
    /// Header counts or flags are inconsistent.
    #[error("invalid DCD header controls")]
    InvalidControls,
    /// Coordinate or free-index record has the wrong length.
    #[error("DCD coordinate record length is inconsistent")]
    CoordinateCount,
    /// A free-atom index is zero, duplicated or beyond the atom table.
    #[error("DCD free-atom index is invalid")]
    InvalidFreeIndex,
    /// Unit-cell values cannot define a crystallographic cell.
    #[error("DCD unit-cell record is invalid")]
    InvalidCell,
    /// A requested memory ceiling is zero or exceeds the process hard limit.
    #[error("invalid DCD memory limit {requested}; it must be at least one byte")]
    InvalidMemoryLimit {
        /// Caller-provided ceiling.
        requested: usize,
    },
    /// Reader workspaces and two reusable frame buffers exceed the ceiling.
    #[error("DCD frame storage requires {required} bytes, over the {limit} byte limit")]
    MemoryLimit {
        /// Required bytes, or `usize::MAX` when dimensions overflow.
        required: usize,
        /// Caller-provided ceiling.
        limit: usize,
    },
}

/// Parses a complete DCD byte stream including fixed atoms and periodic cells.
///
/// # Errors
///
/// Returns a header, record, count, index or unit-cell error without publishing
/// a partial trajectory.
pub fn parse_dcd(bytes: &[u8]) -> Result<DcdTrajectory, DcdError> {
    let endian = detect_endian(bytes)?;
    let mut records = Records::new(bytes, endian);
    let (header, free_indices) = parse_header(&mut records, endian)?;
    let frames = parse_frames(&mut records, &header, &free_indices)?;
    if records.next()?.is_some() {
        return Err(DcdError::InvalidRecord {
            offset: records.cursor,
        });
    }
    Ok(DcdTrajectory { header, frames })
}

fn parse_header(
    records: &mut Records<'_>,
    endian: DcdEndian,
) -> Result<(DcdHeader, Vec<usize>), DcdError> {
    let header_record = records.next()?.ok_or(DcdError::InvalidHeader)?;
    if header_record.len() != 84 || header_record.get(0..4) != Some(b"CORD") {
        return Err(DcdError::InvalidHeader);
    }
    let frame_count = nonnegative(i32_at(header_record, 4, endian)?)?;
    let start_step = i32_at(header_record, 8, endian)?;
    let save_interval = i32_at(header_record, 12, endian)?;
    let fixed_atom_count = nonnegative(i32_at(header_record, 36, endian)?)?;
    let charmm = i32_at(header_record, 80, endian)? != 0;
    let has_unit_cell = charmm && i32_at(header_record, 44, endian)? != 0;
    let has_fourth_dimension = charmm && i32_at(header_record, 48, endian)? == 1;
    let delta_akma = if charmm {
        f64::from(f32_at(header_record, 40, endian)?)
    } else {
        f64_at(header_record, 40, endian)?
    };
    if save_interval <= 0 || !delta_akma.is_finite() || delta_akma <= 0.0 {
        return Err(DcdError::InvalidControls);
    }
    let title_record = records.next()?.ok_or(DcdError::InvalidHeader)?;
    let titles = parse_titles(title_record, endian)?;
    let atom_record = records.next()?.ok_or(DcdError::InvalidHeader)?;
    if atom_record.len() != 4 {
        return Err(DcdError::InvalidControls);
    }
    let atom_count = nonnegative(i32_at(atom_record, 0, endian)?)?;
    if atom_count == 0 || fixed_atom_count > atom_count {
        return Err(DcdError::InvalidControls);
    }
    let free_indices = if fixed_atom_count == 0 {
        Vec::new()
    } else {
        parse_free_indices(
            records.next()?.ok_or(DcdError::InvalidControls)?,
            endian,
            atom_count,
            atom_count - fixed_atom_count,
        )?
    };
    Ok((
        DcdHeader {
            frame_count,
            atom_count,
            start_step,
            save_interval,
            delta_akma,
            fixed_atom_count,
            titles,
            endian,
            charmm,
            has_unit_cell,
            has_fourth_dimension,
        },
        free_indices,
    ))
}

fn parse_frames(
    records: &mut Records<'_>,
    header: &DcdHeader,
    free_indices: &[usize],
) -> Result<Vec<Timestep>, DcdError> {
    let mut frames = Vec::with_capacity(header.frame_count);
    let mut previous: Option<Vec<[f32; 3]>> = None;
    for frame in 0..header.frame_count {
        let cell = if header.has_unit_cell {
            Some(parse_cell(
                records.next()?.ok_or(DcdError::CoordinateCount)?,
                header.endian,
            )?)
        } else {
            None
        };
        let coordinate_count = if frame == 0 || header.fixed_atom_count == 0 {
            header.atom_count
        } else {
            header.atom_count - header.fixed_atom_count
        };
        let x = parse_axis(
            records.next()?.ok_or(DcdError::CoordinateCount)?,
            header.endian,
            coordinate_count,
        )?;
        let y = parse_axis(
            records.next()?.ok_or(DcdError::CoordinateCount)?,
            header.endian,
            coordinate_count,
        )?;
        let z = parse_axis(
            records.next()?.ok_or(DcdError::CoordinateCount)?,
            header.endian,
            coordinate_count,
        )?;
        if header.has_fourth_dimension {
            let fourth = records.next()?.ok_or(DcdError::CoordinateCount)?;
            if fourth.len() != coordinate_count * 4 {
                return Err(DcdError::CoordinateCount);
            }
        }
        let positions = assemble_positions(
            header.atom_count,
            free_indices,
            previous.as_deref(),
            &x,
            &y,
            &z,
        )?;
        previous = Some(positions.clone());
        let frame_number = i32::try_from(frame).map_err(|_| DcdError::InvalidControls)?;
        let step = f64::from(header.start_step)
            + f64::from(frame_number) * f64::from(header.save_interval);
        let dt = f64::from(header.save_interval) * header.delta_akma * AKMA_TO_PS;
        frames.push(Timestep {
            frame,
            time: Some(step * header.delta_akma * AKMA_TO_PS),
            dt: Some(dt),
            positions,
            cell,
            ..Timestep::default()
        });
    }
    Ok(frames)
}

fn detect_endian(bytes: &[u8]) -> Result<DcdEndian, DcdError> {
    let marker = bytes.get(0..4).ok_or(DcdError::InvalidHeader)?;
    if marker == 84_i32.to_le_bytes() {
        Ok(DcdEndian::Little)
    } else if marker == 84_i32.to_be_bytes() {
        Ok(DcdEndian::Big)
    } else {
        Err(DcdError::InvalidHeader)
    }
}

struct Records<'a> {
    bytes: &'a [u8],
    cursor: usize,
    endian: DcdEndian,
}

impl<'a> Records<'a> {
    const fn new(bytes: &'a [u8], endian: DcdEndian) -> Self {
        Self {
            bytes,
            cursor: 0,
            endian,
        }
    }

    fn next(&mut self) -> Result<Option<&'a [u8]>, DcdError> {
        if self.cursor == self.bytes.len() {
            return Ok(None);
        }
        let offset = self.cursor;
        let length = usize::try_from(i32_at(self.bytes, self.cursor, self.endian)?)
            .map_err(|_| DcdError::InvalidRecord { offset })?;
        let start = self.cursor + 4;
        let end = start
            .checked_add(length)
            .filter(|end| {
                end.checked_add(4)
                    .is_some_and(|limit| limit <= self.bytes.len())
            })
            .ok_or(DcdError::InvalidRecord { offset })?;
        let trailing = usize::try_from(i32_at(self.bytes, end, self.endian)?)
            .map_err(|_| DcdError::InvalidRecord { offset })?;
        if trailing != length {
            return Err(DcdError::InvalidRecord { offset });
        }
        self.cursor = end + 4;
        Ok(self.bytes.get(start..end))
    }
}

fn parse_titles(record: &[u8], endian: DcdEndian) -> Result<Vec<Box<str>>, DcdError> {
    if record.len() < 4 {
        return Err(DcdError::InvalidControls);
    }
    let count = nonnegative(i32_at(record, 0, endian)?)?;
    if record.len() != 4 + count * 80 {
        return Err(DcdError::InvalidControls);
    }
    Ok(record[4..]
        .chunks_exact(80)
        .map(|title| String::from_utf8_lossy(title).trim_end().into())
        .collect())
}

fn parse_free_indices(
    record: &[u8],
    endian: DcdEndian,
    atoms: usize,
    expected: usize,
) -> Result<Vec<usize>, DcdError> {
    if record.len() != expected * 4 {
        return Err(DcdError::CoordinateCount);
    }
    let mut unique = BTreeSet::new();
    let mut indices = Vec::with_capacity(expected);
    for offset in (0..record.len()).step_by(4) {
        let raw = i32_at(record, offset, endian)?;
        let index = usize::try_from(raw - 1).map_err(|_| DcdError::InvalidFreeIndex)?;
        if index >= atoms || !unique.insert(index) {
            return Err(DcdError::InvalidFreeIndex);
        }
        indices.push(index);
    }
    Ok(indices)
}

fn parse_axis(record: &[u8], endian: DcdEndian, count: usize) -> Result<Vec<f32>, DcdError> {
    if record.len() != count * 4 {
        return Err(DcdError::CoordinateCount);
    }
    (0..count)
        .map(|index| f32_at(record, index * 4, endian))
        .collect()
}

fn assemble_positions(
    atom_count: usize,
    free_indices: &[usize],
    previous: Option<&[[f32; 3]]>,
    x: &[f32],
    y: &[f32],
    z: &[f32],
) -> Result<Vec<[f32; 3]>, DcdError> {
    if let Some(previous) = previous
        && !free_indices.is_empty()
    {
        let mut positions = previous.to_vec();
        for (coordinate, atom) in free_indices.iter().copied().enumerate() {
            positions[atom] = [x[coordinate], y[coordinate], z[coordinate]];
        }
        return Ok(positions);
    }
    if x.len() != atom_count {
        return Err(DcdError::CoordinateCount);
    }
    Ok((0..atom_count)
        .map(|index| [x[index], y[index], z[index]])
        .collect())
}

fn parse_cell(record: &[u8], endian: DcdEndian) -> Result<UnitCell, DcdError> {
    let values: Vec<f64> = match record.len() {
        48 => (0..6)
            .map(|index| f64_at(record, index * 8, endian))
            .collect::<Result<_, _>>()?,
        24 => (0..6)
            .map(|index| f32_at(record, index * 4, endian).map(f64::from))
            .collect::<Result<_, _>>()?,
        _ => return Err(DcdError::InvalidCell),
    };
    let mut angles = [values[4], values[3], values[1]];
    if angles.iter().all(|angle| angle.abs() <= 1.0) {
        angles = angles.map(|cosine| cosine.clamp(-1.0, 1.0).acos().to_degrees());
    }
    let cell = UnitCell {
        lengths: [values[0], values[2], values[5]],
        angles,
    };
    if cell
        .lengths
        .iter()
        .chain(&cell.angles)
        .any(|value| !value.is_finite())
        || cell.lengths.iter().any(|length| *length <= 0.0)
    {
        return Err(DcdError::InvalidCell);
    }
    Ok(cell)
}

fn nonnegative(value: i32) -> Result<usize, DcdError> {
    usize::try_from(value).map_err(|_| DcdError::InvalidControls)
}

fn i32_at(bytes: &[u8], offset: usize, endian: DcdEndian) -> Result<i32, DcdError> {
    let raw: [u8; 4] = bytes
        .get(offset..offset + 4)
        .and_then(|bytes| bytes.try_into().ok())
        .ok_or(DcdError::InvalidRecord { offset })?;
    Ok(match endian {
        DcdEndian::Little => i32::from_le_bytes(raw),
        DcdEndian::Big => i32::from_be_bytes(raw),
    })
}

fn f32_at(bytes: &[u8], offset: usize, endian: DcdEndian) -> Result<f32, DcdError> {
    let raw: [u8; 4] = bytes
        .get(offset..offset + 4)
        .and_then(|bytes| bytes.try_into().ok())
        .ok_or(DcdError::InvalidRecord { offset })?;
    Ok(match endian {
        DcdEndian::Little => f32::from_le_bytes(raw),
        DcdEndian::Big => f32::from_be_bytes(raw),
    })
}

fn f64_at(bytes: &[u8], offset: usize, endian: DcdEndian) -> Result<f64, DcdError> {
    let raw: [u8; 8] = bytes
        .get(offset..offset + 8)
        .and_then(|bytes| bytes.try_into().ok())
        .ok_or(DcdError::InvalidRecord { offset })?;
    Ok(match endian {
        DcdEndian::Little => f64::from_le_bytes(raw),
        DcdEndian::Big => f64::from_be_bytes(raw),
    })
}

#[cfg(test)]
#[path = "dcd_tests.rs"]
mod tests;
