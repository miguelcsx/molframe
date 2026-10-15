//! Allocation accounting and record helpers for the DCD pull reader.

use std::fs::File;
use std::io::{BufReader, Read};
use std::mem::size_of;

use crate::{Timestep, TrajectoryError};

use super::{DcdReader, FORMAT, FileRecords};
use crate::{DcdEndian, DcdError, DcdHeader};

pub(super) fn build_fixed_indices(
    atom_count: usize,
    free_indices: &[usize],
    memory_limit_bytes: usize,
) -> Result<Vec<usize>, DcdError> {
    if free_indices.is_empty() {
        return Ok(Vec::new());
    }
    let mut free = Vec::new();
    reserve_exact(&mut free, atom_count, memory_limit_bytes)?;
    free.resize(atom_count, false);
    for &index in free_indices {
        free[index] = true;
    }
    let fixed_count = atom_count.saturating_sub(free_indices.len());
    let mut fixed = Vec::new();
    reserve_exact(&mut fixed, fixed_count, memory_limit_bytes)?;
    fixed.extend(
        free.into_iter()
            .enumerate()
            .filter_map(|(index, is_free)| (!is_free).then_some(index)),
    );
    Ok(fixed)
}

pub(super) fn estimated_bytes(
    records: &FileRecords,
    header: &DcdHeader,
    free_capacity: usize,
    fixed_count: usize,
    atom_count: usize,
    record_capacity: usize,
) -> usize {
    size_of::<DcdReader>()
        .saturating_add(records.input.capacity())
        .saturating_add(records.buffer.capacity().max(record_capacity))
        .saturating_add(header_heap_bytes(header))
        .saturating_add(free_capacity.saturating_mul(size_of::<usize>()))
        .saturating_add(fixed_count.saturating_mul(size_of::<usize>()))
        .saturating_add(fixed_count.saturating_mul(size_of::<[f32; 3]>()))
        .saturating_add(
            atom_count
                .saturating_mul(size_of::<[f32; 3]>())
                .saturating_mul(2),
        )
        .saturating_add(size_of::<Timestep>())
}

pub(super) fn header_heap_bytes(header: &DcdHeader) -> usize {
    header
        .titles
        .capacity()
        .saturating_mul(size_of::<Box<str>>())
        .saturating_add(header.titles.iter().map(|title| title.len()).sum::<usize>())
}

pub(super) fn read_marker(
    input: &mut BufReader<File>,
    offset: usize,
) -> Result<Option<[u8; 4]>, DcdError> {
    let mut marker = [0_u8; 4];
    let read = input.read(&mut marker).map_err(|error| io_error(&error))?;
    if read == 0 {
        return Ok(None);
    }
    input
        .read_exact(&mut marker[read..])
        .map_err(|_| DcdError::InvalidRecord { offset })?;
    Ok(Some(marker))
}

pub(super) fn marker_length(
    marker: [u8; 4],
    endian: DcdEndian,
    offset: usize,
) -> Result<usize, DcdError> {
    let signed = match endian {
        DcdEndian::Little => i32::from_le_bytes(marker),
        DcdEndian::Big => i32::from_be_bytes(marker),
    };
    usize::try_from(signed).map_err(|_| DcdError::InvalidRecord { offset })
}

pub(super) fn reserve_exact<T>(
    values: &mut Vec<T>,
    capacity: usize,
    limit: usize,
) -> Result<(), DcdError> {
    if capacity <= values.capacity() {
        return Ok(());
    }
    values
        .try_reserve_exact(capacity.saturating_sub(values.len()))
        .map_err(|_| memory_error(usize::MAX, limit))
}

pub(super) const fn memory_error(required: usize, limit: usize) -> DcdError {
    DcdError::MemoryLimit { required, limit }
}

pub(super) fn io_error(error: &std::io::Error) -> DcdError {
    DcdError::Io { kind: error.kind() }
}

pub(super) fn trajectory_error(error: &DcdError) -> TrajectoryError {
    match error {
        DcdError::Io { kind } => TrajectoryError::SourceIo {
            format: FORMAT,
            kind: *kind,
        },
        DcdError::MemoryLimit { required, limit } => TrajectoryError::MemoryLimit {
            required: *required,
            limit: *limit,
        },
        _ => TrajectoryError::InvalidSource { format: FORMAT },
    }
}
