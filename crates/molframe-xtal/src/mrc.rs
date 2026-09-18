//! MRC2014/CCP4 density maps and coordinate-space sampling.

use crate::numeric::i32_to_usize;
use molframe_core::structure::UnitCell;
const HEADER_BYTES: usize = 1024;
#[path = "mrc_block.rs"]
mod block;
#[path = "mrc_brick.rs"]
mod brick;
#[path = "mrc_sampling.rs"]
mod sampling;
pub use sampling::DensitySampler;
#[path = "mrc_values.rs"]
mod values;
use values::{
    f64_triplet_to_f32, half_to_f32, statistics, to_i32, usize_triplet_to_i32, write_f32,
    write_f32_triplet, write_i32, write_i32_triplet,
};

/// Boundary handling for density-map interpolation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MapBoundary {
    /// Coordinates outside the stored grid have no value.
    Missing,
    /// Wrap indices over the stored grid dimensions.
    Periodic,
}

/// A canonical X-fastest scalar density grid.
#[derive(Clone, Debug, PartialEq)]
pub struct DensityMap {
    /// Stored grid dimensions along crystallographic X, Y, Z.
    pub dimensions: [usize; 3],
    /// First stored grid point along crystallographic X, Y, Z.
    pub starts: [i32; 3],
    /// Unit-cell sampling intervals along crystallographic X, Y, Z.
    pub sampling: [usize; 3],
    /// Unit cell in ångströms and degrees.
    pub cell: UnitCell,
    /// MRC real-space origin in ångströms.
    pub origin: [f64; 3],
    /// International Tables space-group number, or zero when absent.
    pub space_group: i32,
    /// Header labels with trailing padding removed.
    pub labels: Vec<Box<str>>,
    /// Uninterpreted symmetry or vendor-specific extended header.
    pub extended_header: Vec<u8>,
    /// Scalar values in X-fastest, then Y, then Z order.
    pub values: Vec<f32>,
}

/// Invalid or unsupported MRC/CCP4 data.
#[derive(Clone, Debug, thiserror::Error, PartialEq, Eq)]
#[non_exhaustive]
pub enum MrcError {
    /// The file is shorter than its declared content.
    #[error("truncated MRC/CCP4 file")]
    Truncated,
    /// Header dimensions, sampling, cell, axes, or identifiers are invalid.
    #[error("invalid MRC/CCP4 header")]
    InvalidHeader,
    /// Complex transforms cannot be represented by a scalar density map.
    #[error("unsupported scalar-map mode {0}")]
    UnsupportedMode(i32),
    /// Grid size arithmetic exceeded the addressable range.
    #[error("MRC/CCP4 grid size overflow")]
    SizeOverflow,
    /// The host could not reserve the validated output buffer.
    #[error("MRC/CCP4 output allocation exceeds host resources")]
    ResourceLimit,
    /// A value cannot be written as a deterministic scalar map.
    #[error("MRC/CCP4 map contains a non-finite density")]
    NonFiniteDensity,
    /// A requested subvolume lies outside the canonical map dimensions.
    #[error("MRC/CCP4 block is empty or outside the stored grid")]
    InvalidRegion,
    /// A caller requested a zero or excessively large working-memory budget.
    #[error("invalid MRC/CCP4 memory limit {requested}; it must be at least one byte")]
    InvalidMemoryLimit {
        /// Requested limit in bytes.
        requested: usize,
    },
    /// A requested subvolume cannot fit in the configured working-memory budget.
    #[error("MRC/CCP4 block requires {required} bytes; limit is {limit}")]
    MemoryLimit {
        /// Output plus reusable encoded-row workspace in bytes.
        required: usize,
        /// Configured limit in bytes.
        limit: usize,
    },
    /// Seeking or reading the backing source failed.
    #[error("MRC/CCP4 I/O failed: {0:?}")]
    Io(std::io::ErrorKind),
}

#[derive(Clone, Copy, Debug)]
enum Endian {
    Little,
    Big,
}

pub use block::{
    DEFAULT_MRC_BLOCK_MEMORY_LIMIT_BYTES, MrcBlockOptions, MrcBlockReader, MrcMapDescriptor,
};
pub use brick::{
    DEFAULT_MRC_BRICK_PAYLOAD_BYTES, DEFAULT_MRC_BRICK_WORKING_SET_BYTES, MapBrickAddress,
    MapBrickId, MapBrickShape, MrcBrickBudget, MrcBrickDescriptor, MrcBrickError, MrcBrickOptions,
    MrcBrickProvider, ScalarBrickMetadata, ScalarBrickPayload,
};

impl DensityMap {
    /// Reads an MRC2014 or compatible CCP4 scalar map.
    ///
    /// Modes 0, 1, 2, 6, 12 and 101 are converted to `f32`; complex Fourier
    /// modes are refused because silently discarding their imaginary channel
    /// would change the data model.
    ///
    /// # Errors
    ///
    /// Returns an error for truncated data, invalid headers, size overflow, or
    /// a complex/unknown mode.
    pub fn from_mrc_bytes(bytes: &[u8]) -> Result<Self, MrcError> {
        if bytes.len() < HEADER_BYTES {
            return Err(MrcError::Truncated);
        }
        let endian = detect_endian(bytes)?;
        if &bytes[208..212] != b"MAP " {
            return Err(MrcError::InvalidHeader);
        }
        let stored_dimensions = read_positive_triplet(bytes, 0, endian)?;
        let mode = read_i32(bytes, 12, endian)?;
        let stored_starts = read_i32_triplet(bytes, 16, endian)?;
        let sampling = read_positive_triplet(bytes, 28, endian)?;
        let lengths = read_f32_triplet(bytes, 40, endian)?.map(f64::from);
        let angles = read_f32_triplet(bytes, 52, endian)?.map(f64::from);
        let axes = read_i32_triplet(bytes, 64, endian)?;
        if !is_axis_permutation(axes)
            || lengths
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            || angles
                .iter()
                .any(|value| !value.is_finite() || !(0.0..180.0).contains(value))
        {
            return Err(MrcError::InvalidHeader);
        }
        let extended_length =
            usize::try_from(read_i32(bytes, 92, endian)?).map_err(|_| MrcError::InvalidHeader)?;
        let data_offset = HEADER_BYTES
            .checked_add(extended_length)
            .ok_or(MrcError::SizeOverflow)?;
        if data_offset > bytes.len() {
            return Err(MrcError::Truncated);
        }
        let count = product(stored_dimensions)?;
        let dimensions = permute_usize(stored_dimensions, axes);
        let starts = permute_i32(stored_starts, axes);
        let values = decode_values(
            &bytes[data_offset..],
            mode,
            count,
            endian,
            stored_dimensions,
            dimensions,
            axes,
        )?;
        let label_count = usize::try_from(read_i32(bytes, 220, endian)?)
            .map_err(|_| MrcError::InvalidHeader)?
            .min(10);
        let labels = (0..label_count)
            .filter_map(|index| {
                let start = 224 + index * 80;
                let label = String::from_utf8_lossy(&bytes[start..start + 80]);
                let trimmed = label.trim_end_matches(['\0', ' ']);
                (!trimmed.is_empty()).then(|| Box::<str>::from(trimmed))
            })
            .collect();
        Ok(Self {
            dimensions,
            starts,
            sampling,
            cell: UnitCell { lengths, angles },
            origin: read_f32_triplet(bytes, 196, endian)?.map(f64::from),
            space_group: read_i32(bytes, 88, endian)?,
            labels,
            extended_header: bytes[HEADER_BYTES..data_offset].to_vec(),
            values,
        })
    }

    /// Writes a portable little-endian MRC2014 mode-2 scalar map.
    ///
    /// # Errors
    ///
    /// Returns an error when metadata and grid lengths disagree, dimensions
    /// cannot be encoded, or any density is non-finite.
    pub fn to_mrc_bytes(&self) -> Result<Vec<u8>, MrcError> {
        self.validate()?;
        let data_bytes = self
            .values
            .len()
            .checked_mul(4)
            .ok_or(MrcError::SizeOverflow)?;
        let capacity = HEADER_BYTES
            .checked_add(self.extended_header.len())
            .and_then(|value| value.checked_add(data_bytes))
            .ok_or(MrcError::SizeOverflow)?;
        let mut output = vec![0; HEADER_BYTES];
        write_i32_triplet(&mut output, 0, usize_triplet_to_i32(self.dimensions)?);
        write_i32(&mut output, 12, 2);
        write_i32_triplet(&mut output, 16, self.starts);
        write_i32_triplet(&mut output, 28, usize_triplet_to_i32(self.sampling)?);
        write_f32_triplet(&mut output, 40, f64_triplet_to_f32(self.cell.lengths)?);
        write_f32_triplet(&mut output, 52, f64_triplet_to_f32(self.cell.angles)?);
        write_i32_triplet(&mut output, 64, [1, 2, 3]);
        let (minimum, maximum, mean, rms) = statistics(&self.values);
        write_f32_triplet(&mut output, 76, [minimum, maximum, mean]);
        write_i32(&mut output, 88, self.space_group);
        write_i32(&mut output, 92, to_i32(self.extended_header.len())?);
        write_i32(&mut output, 108, 20_141);
        write_f32_triplet(&mut output, 196, f64_triplet_to_f32(self.origin)?);
        output[208..212].copy_from_slice(b"MAP ");
        output[212..216].copy_from_slice(&[0x44, 0x44, 0, 0]);
        write_f32(&mut output, 216, rms);
        let label_count = self.labels.len().min(10);
        write_i32(&mut output, 220, to_i32(label_count)?);
        for (index, label) in self.labels.iter().take(10).enumerate() {
            let source = label.as_bytes();
            let length = source.len().min(80);
            output[224 + index * 80..224 + index * 80 + length].copy_from_slice(&source[..length]);
        }
        let additional = capacity
            .checked_sub(output.len())
            .ok_or(MrcError::SizeOverflow)?;
        output
            .try_reserve_exact(additional)
            .map_err(|_| MrcError::ResourceLimit)?;
        output.extend_from_slice(&self.extended_header);
        for value in &self.values {
            output.extend_from_slice(&value.to_le_bytes());
        }
        Ok(output)
    }

    fn validate(&self) -> Result<(), MrcError> {
        if self.dimensions.contains(&0)
            || self.sampling.contains(&0)
            || product(self.dimensions)? != self.values.len()
            || self.extended_header.len() > i32::MAX as usize
            || self
                .cell
                .lengths
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            || self
                .cell
                .angles
                .iter()
                .any(|value| !value.is_finite() || !(0.0..180.0).contains(value))
        {
            return Err(MrcError::InvalidHeader);
        }
        if self.values.iter().any(|value| !value.is_finite()) {
            return Err(MrcError::NonFiniteDensity);
        }
        Ok(())
    }
}

fn detect_endian(bytes: &[u8]) -> Result<Endian, MrcError> {
    match &bytes[212..214] {
        [0x11, 0x11] => Ok(Endian::Big),
        [0x44, 0x44 | 0x41] => Ok(Endian::Little),
        _ => {
            let little_mode =
                i32::from_le_bytes(bytes[12..16].try_into().map_err(|_| MrcError::Truncated)?);
            let big_mode =
                i32::from_be_bytes(bytes[12..16].try_into().map_err(|_| MrcError::Truncated)?);
            match (known_mode(little_mode), known_mode(big_mode)) {
                (true, false) => Ok(Endian::Little),
                (false, true) => Ok(Endian::Big),
                _ => Err(MrcError::InvalidHeader),
            }
        }
    }
}

fn known_mode(mode: i32) -> bool {
    matches!(mode, 0 | 1 | 2 | 3 | 4 | 6 | 12 | 101)
}

fn read_i32(bytes: &[u8], offset: usize, endian: Endian) -> Result<i32, MrcError> {
    let raw: [u8; 4] = bytes
        .get(offset..offset + 4)
        .ok_or(MrcError::Truncated)?
        .try_into()
        .map_err(|_| MrcError::Truncated)?;
    Ok(match endian {
        Endian::Little => i32::from_le_bytes(raw),
        Endian::Big => i32::from_be_bytes(raw),
    })
}

fn read_f32(bytes: &[u8], offset: usize, endian: Endian) -> Result<f32, MrcError> {
    let raw: [u8; 4] = bytes
        .get(offset..offset + 4)
        .ok_or(MrcError::Truncated)?
        .try_into()
        .map_err(|_| MrcError::Truncated)?;
    Ok(match endian {
        Endian::Little => f32::from_le_bytes(raw),
        Endian::Big => f32::from_be_bytes(raw),
    })
}

fn read_i32_triplet(bytes: &[u8], offset: usize, endian: Endian) -> Result<[i32; 3], MrcError> {
    Ok([
        read_i32(bytes, offset, endian)?,
        read_i32(bytes, offset + 4, endian)?,
        read_i32(bytes, offset + 8, endian)?,
    ])
}
fn read_f32_triplet(bytes: &[u8], offset: usize, endian: Endian) -> Result<[f32; 3], MrcError> {
    Ok([
        read_f32(bytes, offset, endian)?,
        read_f32(bytes, offset + 4, endian)?,
        read_f32(bytes, offset + 8, endian)?,
    ])
}
fn read_positive_triplet(
    bytes: &[u8],
    offset: usize,
    endian: Endian,
) -> Result<[usize; 3], MrcError> {
    read_i32_triplet(bytes, offset, endian)?
        .map(|value| usize::try_from(value).map_err(|_| MrcError::InvalidHeader))
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?
        .try_into()
        .map_err(|_| MrcError::InvalidHeader)
}

fn decode_values(
    bytes: &[u8],
    mode: i32,
    count: usize,
    endian: Endian,
    stored_dimensions: [usize; 3],
    dimensions: [usize; 3],
    axes: [i32; 3],
) -> Result<Vec<f32>, MrcError> {
    let width = match mode {
        0 => 1,
        1 | 6 | 12 => 2,
        2 => 4,
        101 => 0,
        _ => return Err(MrcError::UnsupportedMode(mode)),
    };
    let required = if mode == 101 {
        count.div_ceil(2)
    } else {
        count.checked_mul(width).ok_or(MrcError::SizeOverflow)?
    };
    if bytes.len() < required {
        return Err(MrcError::Truncated);
    }
    let mut output = Vec::new();
    output
        .try_reserve_exact(count)
        .map_err(|_| MrcError::ResourceLimit)?;
    output.resize(count, 0.0);
    let identity = axes == [1, 2, 3];
    let mut stored_index = 0_usize;
    for section in 0..stored_dimensions[2] {
        for row in 0..stored_dimensions[1] {
            for column in 0..stored_dimensions[0] {
                let offset = stored_index
                    .checked_mul(width)
                    .ok_or(MrcError::SizeOverflow)?;
                let value = match mode {
                    0 => f32::from(bytes[stored_index].cast_signed()),
                    1 => f32::from(read_i16(bytes, offset, endian)?),
                    2 => read_f32(bytes, offset, endian)?,
                    6 => f32::from(read_u16(bytes, offset, endian)?),
                    12 => half_to_f32(read_u16(bytes, offset, endian)?),
                    101 => {
                        let byte = bytes[stored_index / 2];
                        f32::from(if stored_index.is_multiple_of(2) {
                            byte & 0x0f
                        } else {
                            byte >> 4
                        })
                    }
                    _ => return Err(MrcError::UnsupportedMode(mode)),
                };
                let target = if identity {
                    stored_index
                } else {
                    let mut xyz = [0; 3];
                    xyz[i32_to_usize(axes[0] - 1)] = column;
                    xyz[i32_to_usize(axes[1] - 1)] = row;
                    xyz[i32_to_usize(axes[2] - 1)] = section;
                    linear(xyz, dimensions)
                };
                output[target] = value;
                stored_index += 1;
            }
        }
    }
    Ok(output)
}

fn read_i16(bytes: &[u8], offset: usize, endian: Endian) -> Result<i16, MrcError> {
    let raw: [u8; 2] = bytes
        .get(offset..offset + 2)
        .ok_or(MrcError::Truncated)?
        .try_into()
        .map_err(|_| MrcError::Truncated)?;
    Ok(match endian {
        Endian::Little => i16::from_le_bytes(raw),
        Endian::Big => i16::from_be_bytes(raw),
    })
}
fn read_u16(bytes: &[u8], offset: usize, endian: Endian) -> Result<u16, MrcError> {
    let raw: [u8; 2] = bytes
        .get(offset..offset + 2)
        .ok_or(MrcError::Truncated)?
        .try_into()
        .map_err(|_| MrcError::Truncated)?;
    Ok(match endian {
        Endian::Little => u16::from_le_bytes(raw),
        Endian::Big => u16::from_be_bytes(raw),
    })
}
fn is_axis_permutation(axes: [i32; 3]) -> bool {
    axes.into_iter().all(|axis| (1..=3).contains(&axis))
        && axes[0] != axes[1]
        && axes[0] != axes[2]
        && axes[1] != axes[2]
}
fn permute_usize(values: [usize; 3], axes: [i32; 3]) -> [usize; 3] {
    let mut output = [0; 3];
    for stored in 0..3 {
        output[i32_to_usize(axes[stored] - 1)] = values[stored];
    }
    output
}
fn permute_i32(values: [i32; 3], axes: [i32; 3]) -> [i32; 3] {
    let mut output = [0; 3];
    for stored in 0..3 {
        output[i32_to_usize(axes[stored] - 1)] = values[stored];
    }
    output
}
fn linear(index: [usize; 3], dimensions: [usize; 3]) -> usize {
    index[0] + dimensions[0] * (index[1] + dimensions[1] * index[2])
}
fn product(dimensions: [usize; 3]) -> Result<usize, MrcError> {
    dimensions
        .into_iter()
        .try_fold(1_usize, |total, dimension| {
            total.checked_mul(dimension).ok_or(MrcError::SizeOverflow)
        })
}
fn resolve_index(
    raw: [i64; 3],
    dimensions: [usize; 3],
    boundary: MapBoundary,
) -> Option<[usize; 3]> {
    let mut output = [0; 3];
    for axis in 0..3 {
        let dimension = i64::try_from(dimensions[axis]).ok()?;
        let value = match boundary {
            MapBoundary::Missing if !(0..dimension).contains(&raw[axis]) => return None,
            MapBoundary::Missing => raw[axis],
            MapBoundary::Periodic => raw[axis].rem_euclid(dimension),
        };
        output[axis] = usize::try_from(value).ok()?;
    }
    Some(output)
}
#[cfg(test)]
#[path = "mrc_tests.rs"]
mod tests;
