//! NAMD binary coordinate snapshots (`.coor` / `.namdbin`).

use crate::Timestep;
use crate::numeric::f32_triplet;

/// Byte order of a NAMD binary snapshot.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NamdEndian {
    /// Least-significant byte first.
    Little,
    /// Most-significant byte first.
    Big,
}

/// Parsed NAMD binary coordinates and their detected byte order.
#[derive(Clone, Debug, PartialEq)]
pub struct NamdBinary {
    /// Coordinate snapshot in canonical ångström units.
    pub frame: Timestep,
    /// Byte order inferred from the count and exact file length.
    pub endian: NamdEndian,
}

/// Invalid NAMD binary coordinate data.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum NamdError {
    /// Count, length or numeric values are invalid.
    #[error("invalid NAMD binary coordinate snapshot")]
    InvalidData,
}

/// Parses a complete NAMD binary coordinate snapshot.
///
/// Both native-endian variants are accepted by validating the encoded atom
/// count against the exact stream length.
///
/// # Errors
///
/// Refuses truncated, overlong, zero-atom and non-finite snapshots.
pub fn parse_namd_binary(bytes: &[u8]) -> Result<NamdBinary, NamdError> {
    let count_bytes: [u8; 4] = bytes
        .get(..4)
        .and_then(|raw| raw.try_into().ok())
        .ok_or(NamdError::InvalidData)?;
    let endian = [NamdEndian::Little, NamdEndian::Big]
        .into_iter()
        .find(|endian| valid_length(bytes.len(), count(count_bytes, *endian)))
        .ok_or(NamdError::InvalidData)?;
    let atoms = usize::try_from(count(count_bytes, endian)).map_err(|_| NamdError::InvalidData)?;
    let mut positions = Vec::with_capacity(atoms);
    for atom in 0..atoms {
        let offset = 4 + atom * 24;
        let position = [
            read_f64(bytes, offset, endian)?,
            read_f64(bytes, offset + 8, endian)?,
            read_f64(bytes, offset + 16, endian)?,
        ];
        if position.iter().any(|value| !value.is_finite()) {
            return Err(NamdError::InvalidData);
        }
        positions.push(f32_triplet(position).ok_or(NamdError::InvalidData)?);
    }
    Ok(NamdBinary {
        frame: Timestep {
            positions,
            ..Timestep::default()
        },
        endian,
    })
}

/// Writes one coordinate snapshot in NAMD's count-plus-f64 layout.
///
/// # Errors
///
/// Refuses zero atoms, non-finite coordinates and counts beyond signed 32-bit.
pub fn write_namd_binary(frame: &Timestep, endian: NamdEndian) -> Result<Vec<u8>, NamdError> {
    let count = i32::try_from(frame.positions.len()).map_err(|_| NamdError::InvalidData)?;
    if count == 0
        || frame
            .positions
            .iter()
            .flatten()
            .any(|value| !value.is_finite())
    {
        return Err(NamdError::InvalidData);
    }
    let mut output = Vec::with_capacity(4 + frame.positions.len() * 24);
    output.extend(match endian {
        NamdEndian::Little => count.to_le_bytes(),
        NamdEndian::Big => count.to_be_bytes(),
    });
    for position in &frame.positions {
        for value in position {
            let value = f64::from(*value);
            output.extend(match endian {
                NamdEndian::Little => value.to_le_bytes(),
                NamdEndian::Big => value.to_be_bytes(),
            });
        }
    }
    Ok(output)
}

fn count(bytes: [u8; 4], endian: NamdEndian) -> i32 {
    match endian {
        NamdEndian::Little => i32::from_le_bytes(bytes),
        NamdEndian::Big => i32::from_be_bytes(bytes),
    }
}

fn valid_length(length: usize, atoms: i32) -> bool {
    usize::try_from(atoms)
        .ok()
        .filter(|atoms| *atoms != 0)
        .and_then(|atoms| atoms.checked_mul(24))
        .and_then(|size| size.checked_add(4))
        == Some(length)
}

fn read_f64(bytes: &[u8], offset: usize, endian: NamdEndian) -> Result<f64, NamdError> {
    let raw: [u8; 8] = bytes
        .get(offset..offset + 8)
        .and_then(|raw| raw.try_into().ok())
        .ok_or(NamdError::InvalidData)?;
    Ok(match endian {
        NamdEndian::Little => f64::from_le_bytes(raw),
        NamdEndian::Big => f64::from_be_bytes(raw),
    })
}

#[cfg(test)]
#[path = "namd_tests.rs"]
mod tests;
