//! Bounded-memory, range-addressable MRC/CCP4 subvolume reads.
//!
//! Construction reads only the fixed header. A block read performs one
//! contiguous read per stored row and uses memory proportional to the returned
//! subvolume plus one encoded row, independent of the complete map size.

use std::fmt;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use pdbiox_core::structure::UnitCell;

use super::{Endian, HEADER_BYTES, MrcError, product};

#[path = "mrc_block_support.rs"]
mod support;
use support::{
    axis_index, clamp_padding, decode_row_value, encoded_span, maximum_encoded_row_bytes,
    parse_header, read_error, resize_bytes_fallibly, resize_fallibly, stored_linear, to_stored,
    validate_limit, validate_region,
};

/// Default maximum for a returned block and its reusable decode workspace.
pub const DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES: usize = 100_000_000;

/// Working-memory policy for range-addressable MRC/CCP4 reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MrcBlockOptions {
    /// Maximum bytes used by the output values and encoded-row workspace.
    pub memory_limit_bytes: usize,
}

impl MrcBlockOptions {
    /// Returns the policy with a caller-selected working-memory ceiling.
    #[must_use]
    pub const fn with_memory_limit(mut self, bytes: usize) -> Self {
        self.memory_limit_bytes = bytes;
        self
    }
}

impl Default for MrcBlockOptions {
    fn default() -> Self {
        Self {
            memory_limit_bytes: DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES,
        }
    }
}

/// Header metadata needed to plan block requests without materialising values.
#[derive(Clone, Debug, PartialEq)]
pub struct MrcMapDescriptor {
    /// Grid dimensions in canonical crystallographic X, Y, Z order.
    pub dimensions: [usize; 3],
    /// First stored grid point in canonical crystallographic X, Y, Z order.
    pub starts: [i32; 3],
    /// Unit-cell sampling intervals along crystallographic X, Y, Z.
    pub sampling: [usize; 3],
    /// Unit cell in ångströms and degrees.
    pub cell: UnitCell,
    /// Real-space origin in ångströms.
    pub origin: [f64; 3],
    /// International Tables space-group number, or zero when absent.
    pub space_group: i32,
    /// Header labels with trailing padding removed.
    pub labels: Vec<Box<str>>,
    /// Byte length of the uninterpreted extended header.
    pub extended_header_bytes: usize,
    /// Scalar storage mode declared by the file.
    pub mode: i32,
    /// Declared whole-map minimum and maximum scalar values.
    pub value_range: [f32; 2],
}

/// Pull reader for canonical subvolumes backed by a seekable source.
pub struct MrcBlockReader<R> {
    source: R,
    descriptor: MrcMapDescriptor,
    stored_dimensions: [usize; 3],
    axes: [i32; 3],
    endian: Endian,
    data_offset: u64,
    memory_limit_bytes: usize,
    encoded_row: Vec<u8>,
}

impl<R> fmt::Debug for MrcBlockReader<R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("MrcBlockReader")
            .field("descriptor", &self.descriptor)
            .field("stored_dimensions", &self.stored_dimensions)
            .field("axes", &self.axes)
            .field("memory_limit_bytes", &self.memory_limit_bytes)
            .finish_non_exhaustive()
    }
}

impl MrcBlockReader<File> {
    /// Opens a local map without reading its value payload.
    ///
    /// # Errors
    ///
    /// Returns an error when the file cannot be opened or its header is invalid.
    pub fn open(path: impl AsRef<Path>, options: MrcBlockOptions) -> Result<Self, MrcError> {
        let source = File::open(path).map_err(MrcError::from)?;
        Self::new(source, options)
    }
}

impl<R: Read + Seek> MrcBlockReader<R> {
    /// Reads and validates the fixed header from a seekable source.
    ///
    /// # Errors
    ///
    /// Returns an error for malformed metadata, unsupported scalar modes, I/O
    /// failure, or a memory limit outside the accepted range.
    pub fn new(mut source: R, options: MrcBlockOptions) -> Result<Self, MrcError> {
        validate_limit(options.memory_limit_bytes)?;
        source.seek(SeekFrom::Start(0)).map_err(MrcError::from)?;
        let mut header = [0_u8; HEADER_BYTES];
        source.read_exact(&mut header).map_err(read_error)?;
        let parsed = parse_header(&header)?;
        Ok(Self {
            source,
            descriptor: parsed.descriptor,
            stored_dimensions: parsed.stored_dimensions,
            axes: parsed.axes,
            endian: parsed.endian,
            data_offset: parsed.data_offset,
            memory_limit_bytes: options.memory_limit_bytes,
            encoded_row: Vec::new(),
        })
    }

    /// Returns immutable metadata for request planning.
    #[must_use]
    pub const fn descriptor(&self) -> &MrcMapDescriptor {
        &self.descriptor
    }

    /// Returns the underlying source after all buffered state is discarded.
    #[must_use]
    pub fn into_inner(self) -> R {
        self.source
    }

    /// Decodes one canonical X-fastest block into a caller-owned reusable buffer.
    ///
    /// `origin` and `dimensions` use canonical X, Y, Z grid coordinates. The
    /// output is cleared and then filled in X-fastest order. Existing capacity
    /// is retained, so equally sized requests do not allocate after warm-up.
    ///
    /// # Errors
    ///
    /// Returns an error when the region is empty or outside the map, arithmetic
    /// overflows, the configured memory budget is exceeded, or source I/O fails.
    pub fn read_block_into(
        &mut self,
        origin: [usize; 3],
        dimensions: [usize; 3],
        output: &mut Vec<f32>,
    ) -> Result<(), MrcError> {
        validate_region(origin, dimensions, self.descriptor.dimensions)?;
        let value_count = product(dimensions)?;
        let output_bytes = value_count
            .checked_mul(size_of::<f32>())
            .ok_or(MrcError::SizeOverflow)?;
        let stored_shape = to_stored(dimensions, self.axes)?;
        let row_bytes = maximum_encoded_row_bytes(self.descriptor.mode, stored_shape[0])?;
        let required = output_bytes
            .checked_add(row_bytes)
            .ok_or(MrcError::SizeOverflow)?;
        if required > self.memory_limit_bytes {
            return Err(MrcError::MemoryLimit {
                required,
                limit: self.memory_limit_bytes,
            });
        }
        resize_fallibly(output, value_count)?;
        resize_bytes_fallibly(&mut self.encoded_row, row_bytes)?;
        self.read_region_into(origin, dimensions, [0; 3], dimensions, output)
    }

    /// Decodes an interior block plus a symmetric clamped halo.
    ///
    /// Values outside the stored map repeat the nearest boundary voxel. The
    /// returned allocation is the only scalar payload; the encoded-row buffer
    /// is retained by the reader and bounded by the same memory policy.
    ///
    /// # Errors
    ///
    /// Returns an error for an invalid interior, shape arithmetic overflow,
    /// budget exhaustion, allocation failure, or source I/O failure.
    pub fn read_block_with_halo_into(
        &mut self,
        origin: [usize; 3],
        interior: [usize; 3],
        halo: usize,
        output: &mut Vec<f32>,
    ) -> Result<[usize; 3], MrcError> {
        validate_region(origin, interior, self.descriptor.dimensions)?;
        let doubled = halo.checked_mul(2).ok_or(MrcError::SizeOverflow)?;
        let [stored_x, stored_y, stored_z] =
            interior.map(|value| value.checked_add(doubled).ok_or(MrcError::SizeOverflow));
        let stored = [stored_x?, stored_y?, stored_z?];
        let value_count = product(stored)?;
        let output_bytes = value_count
            .checked_mul(size_of::<f32>())
            .ok_or(MrcError::SizeOverflow)?;
        let stored_axes = to_stored(stored, self.axes)?;
        let row_bytes = maximum_encoded_row_bytes(self.descriptor.mode, stored_axes[0])?;
        let required = output_bytes
            .checked_add(row_bytes)
            .ok_or(MrcError::SizeOverflow)?;
        if required > self.memory_limit_bytes {
            return Err(MrcError::MemoryLimit {
                required,
                limit: self.memory_limit_bytes,
            });
        }
        resize_fallibly(output, value_count)?;
        resize_bytes_fallibly(&mut self.encoded_row, row_bytes)?;

        let mut source_origin = [0; 3];
        let mut source_shape = [0; 3];
        let mut destination = [0; 3];
        for axis in 0..3 {
            let left = halo.min(origin[axis]);
            source_origin[axis] = origin[axis] - left;
            destination[axis] = halo - left;
            let requested_end = origin[axis]
                .checked_add(interior[axis])
                .and_then(|value| value.checked_add(halo))
                .ok_or(MrcError::SizeOverflow)?;
            let source_end = requested_end.min(self.descriptor.dimensions[axis]);
            source_shape[axis] = source_end - source_origin[axis];
        }
        self.read_region_into(source_origin, source_shape, destination, stored, output)?;
        clamp_padding(output, stored, destination, source_shape)?;
        Ok(stored)
    }

    fn read_region_into(
        &mut self,
        origin: [usize; 3],
        dimensions: [usize; 3],
        destination: [usize; 3],
        output_shape: [usize; 3],
        output: &mut [f32],
    ) -> Result<(), MrcError> {
        let stored_origin = to_stored(origin, self.axes)?;
        let stored_shape = to_stored(dimensions, self.axes)?;
        let strides = [
            1,
            output_shape[0],
            output_shape[0]
                .checked_mul(output_shape[1])
                .ok_or(MrcError::SizeOverflow)?,
        ];
        let column_stride = strides[axis_index(self.axes[0])?];
        let row_stride = strides[axis_index(self.axes[1])?];
        let section_stride = strides[axis_index(self.axes[2])?];
        let destination_base = destination[0]
            .checked_add(
                destination[1]
                    .checked_mul(strides[1])
                    .ok_or(MrcError::SizeOverflow)?,
            )
            .and_then(|value| {
                destination[2]
                    .checked_mul(strides[2])
                    .and_then(|offset| value.checked_add(offset))
            })
            .ok_or(MrcError::SizeOverflow)?;
        for section in 0..stored_shape[2] {
            for row in 0..stored_shape[1] {
                let stored_index = stored_linear(
                    [
                        stored_origin[0],
                        stored_origin[1] + row,
                        stored_origin[2] + section,
                    ],
                    self.stored_dimensions,
                )?;
                self.read_encoded_row(stored_index, stored_shape[0])?;
                let output_base = destination_base
                    .checked_add(
                        row.checked_mul(row_stride)
                            .and_then(|value| {
                                section
                                    .checked_mul(section_stride)
                                    .and_then(|section| value.checked_add(section))
                            })
                            .ok_or(MrcError::SizeOverflow)?,
                    )
                    .ok_or(MrcError::SizeOverflow)?;
                for column in 0..stored_shape[0] {
                    let output_index = column
                        .checked_mul(column_stride)
                        .and_then(|value| output_base.checked_add(value))
                        .ok_or(MrcError::SizeOverflow)?;
                    output[output_index] = decode_row_value(
                        &self.encoded_row,
                        self.descriptor.mode,
                        column,
                        stored_index & 1,
                        self.endian,
                    )?;
                }
            }
        }
        Ok(())
    }

    fn read_encoded_row(&mut self, first_value: usize, count: usize) -> Result<(), MrcError> {
        let (relative_offset, byte_count) = encoded_span(self.descriptor.mode, first_value, count)?;
        resize_bytes_fallibly(&mut self.encoded_row, byte_count)?;
        let offset = self
            .data_offset
            .checked_add(u64::try_from(relative_offset).map_err(|_| MrcError::SizeOverflow)?)
            .ok_or(MrcError::SizeOverflow)?;
        self.source
            .seek(SeekFrom::Start(offset))
            .map_err(MrcError::from)?;
        self.source
            .read_exact(&mut self.encoded_row)
            .map_err(read_error)
    }
}

#[cfg(test)]
#[path = "mrc_block_tests.rs"]
mod tests;
