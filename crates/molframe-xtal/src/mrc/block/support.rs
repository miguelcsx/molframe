//! Header parsing, scalar codecs, and checked layout helpers for block reads.

use std::io;

use molframe_core::structure::UnitCell;

use super::MrcMapDescriptor;
use crate::mrc::{
    Endian, HEADER_BYTES, MrcError, detect_endian, half_to_f32, is_axis_permutation, permute_i32,
    permute_usize, product, read_f32, read_f32_triplet, read_i16, read_i32, read_i32_triplet,
    read_positive_triplet, read_u16,
};

pub(super) struct ParsedHeader {
    pub(super) descriptor: MrcMapDescriptor,
    pub(super) stored_dimensions: [usize; 3],
    pub(super) axes: [i32; 3],
    pub(super) endian: Endian,
    pub(super) data_offset: u64,
}

pub(super) fn parse_header(bytes: &[u8; HEADER_BYTES]) -> Result<ParsedHeader, MrcError> {
    let endian = detect_endian(bytes)?;
    if &bytes[208..212] != b"MAP " {
        return Err(MrcError::InvalidHeader);
    }
    let stored_dimensions = read_positive_triplet(bytes, 0, endian)?;
    product(stored_dimensions)?;
    let mode = read_i32(bytes, 12, endian)?;
    scalar_width(mode)?;
    let declared_min = read_f32(bytes, 76, endian)?;
    let declared_max = read_f32(bytes, 80, endian)?;
    if !declared_min.is_finite() || !declared_max.is_finite() || declared_min > declared_max {
        return Err(MrcError::InvalidHeader);
    }
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
    let extended_header_bytes =
        usize::try_from(read_i32(bytes, 92, endian)?).map_err(|_| MrcError::InvalidHeader)?;
    let data_offset = HEADER_BYTES
        .checked_add(extended_header_bytes)
        .ok_or(MrcError::SizeOverflow)?;
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
    Ok(ParsedHeader {
        descriptor: MrcMapDescriptor {
            dimensions: permute_usize(stored_dimensions, axes),
            starts: permute_i32(stored_starts, axes),
            sampling,
            cell: UnitCell { lengths, angles },
            origin: read_f32_triplet(bytes, 196, endian)?.map(f64::from),
            space_group: read_i32(bytes, 88, endian)?,
            labels,
            extended_header_bytes,
            mode,
            value_range: [declared_min, declared_max],
        },
        stored_dimensions,
        axes,
        endian,
        data_offset: u64::try_from(data_offset).map_err(|_| MrcError::SizeOverflow)?,
    })
}

pub(super) fn validate_limit(limit: usize) -> Result<(), MrcError> {
    if limit == 0 {
        return Err(MrcError::InvalidMemoryLimit { requested: limit });
    }
    Ok(())
}

pub(super) fn validate_region(
    origin: [usize; 3],
    shape: [usize; 3],
    map: [usize; 3],
) -> Result<(), MrcError> {
    if shape.contains(&0) {
        return Err(MrcError::InvalidRegion);
    }
    for axis in 0..3 {
        let Some(end) = origin[axis].checked_add(shape[axis]) else {
            return Err(MrcError::SizeOverflow);
        };
        if end > map[axis] {
            return Err(MrcError::InvalidRegion);
        }
    }
    Ok(())
}

pub(super) fn to_stored(canonical: [usize; 3], axes: [i32; 3]) -> Result<[usize; 3], MrcError> {
    Ok([
        canonical[axis_index(axes[0])?],
        canonical[axis_index(axes[1])?],
        canonical[axis_index(axes[2])?],
    ])
}

pub(super) fn axis_index(axis: i32) -> Result<usize, MrcError> {
    usize::try_from(axis - 1).map_err(|_| MrcError::InvalidHeader)
}

pub(super) fn stored_linear(index: [usize; 3], dimensions: [usize; 3]) -> Result<usize, MrcError> {
    index[1]
        .checked_add(
            dimensions[1]
                .checked_mul(index[2])
                .ok_or(MrcError::SizeOverflow)?,
        )
        .and_then(|value| dimensions[0].checked_mul(value))
        .and_then(|value| index[0].checked_add(value))
        .ok_or(MrcError::SizeOverflow)
}

pub(super) fn clamp_padding(
    values: &mut [f32],
    shape: [usize; 3],
    valid_origin: [usize; 3],
    valid_shape: [usize; 3],
) -> Result<(), MrcError> {
    let valid_end = [
        valid_origin[0]
            .checked_add(valid_shape[0])
            .ok_or(MrcError::SizeOverflow)?,
        valid_origin[1]
            .checked_add(valid_shape[1])
            .ok_or(MrcError::SizeOverflow)?,
        valid_origin[2]
            .checked_add(valid_shape[2])
            .ok_or(MrcError::SizeOverflow)?,
    ];
    let stride_z = shape[0]
        .checked_mul(shape[1])
        .ok_or(MrcError::SizeOverflow)?;
    for z in 0..shape[2] {
        let source_z = z.clamp(valid_origin[2], valid_end[2] - 1);
        for y in 0..shape[1] {
            let source_y = y.clamp(valid_origin[1], valid_end[1] - 1);
            for x in 0..shape[0] {
                if (valid_origin[0]..valid_end[0]).contains(&x)
                    && (valid_origin[1]..valid_end[1]).contains(&y)
                    && (valid_origin[2]..valid_end[2]).contains(&z)
                {
                    continue;
                }
                let source_x = x.clamp(valid_origin[0], valid_end[0] - 1);
                let source = source_x
                    .checked_add(source_y * shape[0])
                    .and_then(|value| value.checked_add(source_z * stride_z))
                    .ok_or(MrcError::SizeOverflow)?;
                let target = x
                    .checked_add(y * shape[0])
                    .and_then(|value| value.checked_add(z * stride_z))
                    .ok_or(MrcError::SizeOverflow)?;
                values[target] = values[source];
            }
        }
    }
    Ok(())
}

fn scalar_width(mode: i32) -> Result<usize, MrcError> {
    match mode {
        0 => Ok(1),
        1 | 6 | 12 => Ok(2),
        2 => Ok(4),
        101 => Ok(0),
        _ => Err(MrcError::UnsupportedMode(mode)),
    }
}

pub(super) fn maximum_encoded_row_bytes(mode: i32, count: usize) -> Result<usize, MrcError> {
    let width = scalar_width(mode)?;
    if mode == 101 {
        return count
            .checked_add(1)
            .map(|value| value.div_ceil(2))
            .ok_or(MrcError::SizeOverflow);
    }
    count.checked_mul(width).ok_or(MrcError::SizeOverflow)
}

pub(super) fn encoded_span(
    mode: i32,
    first: usize,
    count: usize,
) -> Result<(usize, usize), MrcError> {
    let width = scalar_width(mode)?;
    if mode == 101 {
        let byte_offset = first / 2;
        let nibble_count = (first & 1)
            .checked_add(count)
            .ok_or(MrcError::SizeOverflow)?;
        return Ok((byte_offset, nibble_count.div_ceil(2)));
    }
    Ok((
        first.checked_mul(width).ok_or(MrcError::SizeOverflow)?,
        count.checked_mul(width).ok_or(MrcError::SizeOverflow)?,
    ))
}

pub(super) fn decode_row_value(
    bytes: &[u8],
    mode: i32,
    index: usize,
    first_parity: usize,
    endian: Endian,
) -> Result<f32, MrcError> {
    match mode {
        0 => bytes
            .get(index)
            .copied()
            .map(|value| f32::from(value.cast_signed()))
            .ok_or(MrcError::Truncated),
        1 => read_i16(
            bytes,
            index.checked_mul(2).ok_or(MrcError::SizeOverflow)?,
            endian,
        )
        .map(f32::from),
        2 => read_f32(
            bytes,
            index.checked_mul(4).ok_or(MrcError::SizeOverflow)?,
            endian,
        ),
        6 => read_u16(
            bytes,
            index.checked_mul(2).ok_or(MrcError::SizeOverflow)?,
            endian,
        )
        .map(f32::from),
        12 => read_u16(
            bytes,
            index.checked_mul(2).ok_or(MrcError::SizeOverflow)?,
            endian,
        )
        .map(half_to_f32),
        101 => {
            let nibble = first_parity
                .checked_add(index)
                .ok_or(MrcError::SizeOverflow)?;
            let byte = bytes.get(nibble / 2).copied().ok_or(MrcError::Truncated)?;
            Ok(f32::from(if nibble & 1 == 0 {
                byte & 0x0f
            } else {
                byte >> 4
            }))
        }
        _ => Err(MrcError::UnsupportedMode(mode)),
    }
}

pub(super) fn resize_fallibly(values: &mut Vec<f32>, length: usize) -> Result<(), MrcError> {
    values.clear();
    if length > values.capacity() {
        values
            .try_reserve_exact(length - values.capacity())
            .map_err(|_| MrcError::ResourceLimit)?;
    }
    values.resize(length, 0.0);
    Ok(())
}

pub(super) fn resize_bytes_fallibly(values: &mut Vec<u8>, length: usize) -> Result<(), MrcError> {
    values.clear();
    if length > values.capacity() {
        values
            .try_reserve_exact(length - values.capacity())
            .map_err(|_| MrcError::ResourceLimit)?;
    }
    values.resize(length, 0);
    Ok(())
}

impl From<io::Error> for MrcError {
    fn from(error: io::Error) -> Self {
        Self::Io(error.kind())
    }
}

pub(super) fn read_error(error: io::Error) -> MrcError {
    if error.kind() == io::ErrorKind::UnexpectedEof {
        MrcError::Truncated
    } else {
        MrcError::from(error)
    }
}
