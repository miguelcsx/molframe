//! Bounded pull-based DCD decoding from a file.
//!
//! Decoding is linear in the number of coordinates. Storage is proportional
//! to one input record, two coordinate frames and the fixed/free atom tables;
//! it does not depend on the number of frames in the file.

use std::fs::File;
use std::io::{BufReader, Read};
use std::mem::size_of;
use std::path::Path;

use crate::{
    RandomAccess, Timestep, TrajectoryError, TrajectoryReader, TrajectoryReaderOptions, Units,
};

use super::{
    AKMA_TO_PS, DcdEndian, DcdError, DcdHeader, f32_at, i32_at, nonnegative, parse_cell,
    parse_free_indices, parse_titles,
};

#[path = "reader_support.rs"]
mod support;

use support::{
    build_fixed_indices, estimated_bytes, header_heap_bytes, io_error, marker_length, memory_error,
    read_marker, reserve_exact, trajectory_error,
};

const FORMAT: &str = "dcd";
const INPUT_BUFFER_BYTES: usize = 64 * 1024;

/// File-backed DCD source retaining bounded reusable buffers.
#[derive(Debug)]
pub struct DcdReader {
    records: FileRecords,
    header: DcdHeader,
    free_indices: Vec<usize>,
    fixed_indices: Vec<usize>,
    fixed_positions: Vec<[f32; 3]>,
    decoded: Vec<[f32; 3]>,
    next_frame: usize,
    memory_limit_bytes: usize,
    eof_checked: bool,
}

impl DcdReader {
    /// Opens a DCD file under the default 100 MB reader ceiling.
    ///
    /// # Errors
    ///
    /// Returns an I/O, header, allocation or memory-ceiling error.
    pub fn open(path: &Path) -> Result<Self, DcdError> {
        Self::open_with_memory_limit(path, TrajectoryReaderOptions::DEFAULT_MEMORY_LIMIT_BYTES)
    }

    /// Opens a DCD file under an explicit ceiling of at most 500 MB.
    ///
    /// The ceiling includes the input record buffer, format metadata, fixed
    /// atom state, the private decode frame and one caller-owned output frame.
    ///
    /// # Errors
    ///
    /// Returns an I/O, header, allocation or invalid-memory-limit error.
    pub fn open_with_memory_limit(
        path: &Path,
        memory_limit_bytes: usize,
    ) -> Result<Self, DcdError> {
        if memory_limit_bytes == 0 {
            return Err(DcdError::InvalidMemoryLimit {
                requested: memory_limit_bytes,
            });
        }
        let file = File::open(path).map_err(|error| io_error(&error))?;
        let mut records = FileRecords::new(file, memory_limit_bytes)?;
        let (header, free_indices) = read_header(&mut records)?;
        let fixed_count = header.fixed_atom_count;
        let atom_count = header.atom_count;
        let axis_bytes = atom_count
            .checked_mul(size_of::<f32>())
            .ok_or(memory_error(usize::MAX, memory_limit_bytes))?;
        let requested = estimated_bytes(
            &records,
            &header,
            free_indices.capacity(),
            fixed_count,
            atom_count,
            axis_bytes,
        );
        if requested > memory_limit_bytes {
            return Err(memory_error(requested, memory_limit_bytes));
        }
        records.reserve_payload(axis_bytes)?;
        let fixed_indices = build_fixed_indices(atom_count, &free_indices, memory_limit_bytes)?;
        let mut fixed_positions = Vec::new();
        reserve_exact(&mut fixed_positions, fixed_count, memory_limit_bytes)?;
        let mut decoded = Vec::new();
        reserve_exact(&mut decoded, atom_count, memory_limit_bytes)?;
        let reader = Self {
            records,
            header,
            free_indices,
            fixed_indices,
            fixed_positions,
            decoded,
            next_frame: 0,
            memory_limit_bytes,
            eof_checked: false,
        };
        let actual = reader.total_bytes(atom_count);
        if actual > memory_limit_bytes {
            return Err(memory_error(actual, memory_limit_bytes));
        }
        Ok(reader)
    }

    /// Parsed controls and title metadata.
    #[must_use]
    pub const fn header(&self) -> &DcdHeader {
        &self.header
    }

    /// Reader-owned allocation capacity, excluding the caller's timestep.
    #[must_use]
    pub fn workspace_bytes(&self) -> usize {
        self.reader_bytes()
    }

    /// Reads one frame into caller-owned storage while retaining allocations.
    ///
    /// # Errors
    ///
    /// Returns an I/O, malformed-record, allocation or memory-ceiling error
    /// without publishing a partially decoded frame.
    pub fn read_next_frame(&mut self, timestep: &mut Timestep) -> Result<bool, DcdError> {
        if self.next_frame == self.header.frame_count {
            if !self.eof_checked {
                if self.records.next()?.is_some() {
                    return Err(DcdError::InvalidRecord {
                        offset: self.records.offset(),
                    });
                }
                self.eof_checked = true;
            }
            return Ok(false);
        }
        self.prepare_decode(timestep.positions.capacity())?;
        let cell = if self.header.has_unit_cell {
            let record = self.records.next()?.ok_or(DcdError::CoordinateCount)?;
            Some(parse_cell(record, self.header.endian)?)
        } else {
            None
        };
        let coordinate_count = if self.next_frame == 0 || self.header.fixed_atom_count == 0 {
            self.header.atom_count
        } else {
            self.header.atom_count - self.header.fixed_atom_count
        };
        self.seed_fixed_positions();
        let axis_indices = if coordinate_count == self.header.atom_count {
            &[][..]
        } else {
            self.free_indices.as_slice()
        };
        for axis in 0..3 {
            let record = self.records.next()?.ok_or(DcdError::CoordinateCount)?;
            decode_axis(
                record,
                self.header.endian,
                coordinate_count,
                axis,
                axis_indices,
                &mut self.decoded,
            )?;
        }
        if self.header.has_fourth_dimension {
            let fourth = self.records.next()?.ok_or(DcdError::CoordinateCount)?;
            let expected = coordinate_count
                .checked_mul(size_of::<f32>())
                .ok_or(DcdError::CoordinateCount)?;
            if fourth.len() != expected {
                return Err(DcdError::CoordinateCount);
            }
        }
        if self.next_frame == 0 {
            for &index in &self.fixed_indices {
                self.fixed_positions.push(self.decoded[index]);
            }
        }
        let frame_number = i32::try_from(self.next_frame).map_err(|_| DcdError::InvalidControls)?;
        let step = f64::from(self.header.start_step)
            + f64::from(frame_number) * f64::from(self.header.save_interval);
        let dt = f64::from(self.header.save_interval) * self.header.delta_akma * AKMA_TO_PS;
        std::mem::swap(&mut timestep.positions, &mut self.decoded);
        timestep.frame = self.next_frame;
        timestep.time = Some(step * self.header.delta_akma * AKMA_TO_PS);
        timestep.dt = Some(dt);
        timestep.velocities = None;
        timestep.forces = None;
        timestep.cell = cell;
        timestep.data.clear();
        self.next_frame = self
            .next_frame
            .checked_add(1)
            .ok_or(DcdError::InvalidControls)?;
        Ok(true)
    }

    fn prepare_decode(&mut self, output_capacity: usize) -> Result<(), DcdError> {
        let output_capacity = output_capacity.max(self.header.atom_count);
        let required = self.total_bytes(output_capacity);
        if required > self.memory_limit_bytes {
            return Err(memory_error(required, self.memory_limit_bytes));
        }
        reserve_exact(
            &mut self.decoded,
            self.header.atom_count,
            self.memory_limit_bytes,
        )?;
        self.decoded.clear();
        self.decoded.resize(self.header.atom_count, [0.0; 3]);
        let actual = self.total_bytes(output_capacity);
        if actual > self.memory_limit_bytes {
            return Err(memory_error(actual, self.memory_limit_bytes));
        }
        Ok(())
    }

    fn seed_fixed_positions(&mut self) {
        if self.next_frame == 0 {
            return;
        }
        for (&index, &position) in self.fixed_indices.iter().zip(&self.fixed_positions) {
            self.decoded[index] = position;
        }
    }

    fn reader_bytes(&self) -> usize {
        size_of::<Self>()
            .saturating_add(self.records.allocated_bytes())
            .saturating_add(header_heap_bytes(&self.header))
            .saturating_add(
                self.free_indices
                    .capacity()
                    .saturating_mul(size_of::<usize>()),
            )
            .saturating_add(
                self.fixed_indices
                    .capacity()
                    .saturating_mul(size_of::<usize>()),
            )
            .saturating_add(
                self.fixed_positions
                    .capacity()
                    .saturating_mul(size_of::<[f32; 3]>()),
            )
            .saturating_add(
                self.decoded
                    .capacity()
                    .saturating_mul(size_of::<[f32; 3]>()),
            )
    }

    fn total_bytes(&self, output_capacity: usize) -> usize {
        self.reader_bytes()
            .saturating_add(size_of::<Timestep>())
            .saturating_add(output_capacity.saturating_mul(size_of::<[f32; 3]>()))
    }
}

impl TrajectoryReader for DcdReader {
    fn format(&self) -> &'static str {
        FORMAT
    }

    fn n_atoms(&self) -> usize {
        self.header.atom_count
    }

    fn n_frames(&self) -> Option<usize> {
        Some(self.header.frame_count)
    }

    fn units(&self) -> Units {
        Units::CANONICAL
    }

    fn random_access(&self) -> RandomAccess {
        RandomAccess::None
    }

    fn read_next(&mut self, timestep: &mut Timestep) -> Result<bool, TrajectoryError> {
        self.read_next_frame(timestep)
            .map_err(|error| trajectory_error(&error))
    }

    fn read_next_bounded(
        &mut self,
        timestep: &mut Timestep,
        bytes: usize,
    ) -> Result<bool, TrajectoryError> {
        let previous = self.memory_limit_bytes;
        self.memory_limit_bytes = previous.min(bytes);
        let result = self.read_next(timestep);
        self.memory_limit_bytes = previous;
        result
    }

    fn seek(&mut self, _frame: usize) -> Result<(), TrajectoryError> {
        Err(TrajectoryError::RandomAccessUnavailable)
    }
}

#[derive(Debug)]
struct FileRecords {
    input: BufReader<File>,
    endian: DcdEndian,
    buffer: Vec<u8>,
    pending_length: Option<usize>,
    offset: usize,
    payload_limit: usize,
}

impl FileRecords {
    fn new(file: File, memory_limit_bytes: usize) -> Result<Self, DcdError> {
        let mut input = BufReader::with_capacity(INPUT_BUFFER_BYTES, file);
        let mut marker = [0_u8; 4];
        input
            .read_exact(&mut marker)
            .map_err(|error| io_error(&error))?;
        let endian = if marker == 84_i32.to_le_bytes() {
            DcdEndian::Little
        } else if marker == 84_i32.to_be_bytes() {
            DcdEndian::Big
        } else {
            return Err(DcdError::InvalidHeader);
        };
        let fixed = size_of::<Self>().saturating_add(input.capacity());
        let payload_limit = memory_limit_bytes.saturating_sub(fixed);
        Ok(Self {
            input,
            endian,
            buffer: Vec::new(),
            pending_length: Some(84),
            offset: 0,
            payload_limit,
        })
    }

    fn next(&mut self) -> Result<Option<&[u8]>, DcdError> {
        let record_offset = self.offset;
        let length = if let Some(length) = self.pending_length.take() {
            length
        } else {
            let Some(marker) = read_marker(&mut self.input, record_offset)? else {
                return Ok(None);
            };
            marker_length(marker, self.endian, record_offset)?
        };
        if length > self.payload_limit {
            return Err(memory_error(length, self.payload_limit));
        }
        reserve_exact(&mut self.buffer, length, self.payload_limit)?;
        self.buffer.clear();
        self.buffer.resize(length, 0);
        self.input
            .read_exact(&mut self.buffer)
            .map_err(|_| DcdError::InvalidRecord {
                offset: record_offset,
            })?;
        let mut trailing = [0_u8; 4];
        self.input
            .read_exact(&mut trailing)
            .map_err(|_| DcdError::InvalidRecord {
                offset: record_offset,
            })?;
        if marker_length(trailing, self.endian, record_offset)? != length {
            return Err(DcdError::InvalidRecord {
                offset: record_offset,
            });
        }
        self.offset =
            self.offset
                .checked_add(length.saturating_add(8))
                .ok_or(DcdError::InvalidRecord {
                    offset: record_offset,
                })?;
        Ok(Some(&self.buffer))
    }

    fn reserve_payload(&mut self, bytes: usize) -> Result<(), DcdError> {
        if bytes > self.payload_limit {
            return Err(memory_error(bytes, self.payload_limit));
        }
        reserve_exact(&mut self.buffer, bytes, self.payload_limit)
    }

    const fn offset(&self) -> usize {
        self.offset
    }

    fn allocated_bytes(&self) -> usize {
        self.input.capacity().saturating_add(self.buffer.capacity())
    }
}

fn read_header(records: &mut FileRecords) -> Result<(DcdHeader, Vec<usize>), DcdError> {
    let endian = records.endian;
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
        super::f64_at(header_record, 40, endian)?
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
        let record = records.next()?.ok_or(DcdError::InvalidControls)?;
        parse_free_indices(record, endian, atom_count, atom_count - fixed_atom_count)?
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

fn decode_axis(
    record: &[u8],
    endian: DcdEndian,
    count: usize,
    axis: usize,
    free_indices: &[usize],
    positions: &mut [[f32; 3]],
) -> Result<(), DcdError> {
    let expected = count
        .checked_mul(size_of::<f32>())
        .ok_or(DcdError::CoordinateCount)?;
    if record.len() != expected {
        return Err(DcdError::CoordinateCount);
    }
    for coordinate in 0..count {
        let atom = if free_indices.is_empty() {
            coordinate
        } else {
            *free_indices
                .get(coordinate)
                .ok_or(DcdError::CoordinateCount)?
        };
        positions[atom][axis] = f32_at(record, coordinate * size_of::<f32>(), endian)?;
    }
    Ok(())
}

#[cfg(test)]
#[path = "reader_tests.rs"]
mod tests;
