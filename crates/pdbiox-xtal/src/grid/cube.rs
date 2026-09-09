//! Gaussian cube reader.
//!
//! A cube is two comment lines, an atom count with the grid origin, three axis
//! lines giving a voxel count and that axis's step vector, one line per atom,
//! and then the values with the third axis varying fastest. A negative atom
//! count marks an orbital cube, which carries an extra line of orbital indices
//! and may store several values per voxel.
//!
//! Geometry is in bohr unless a voxel count is written negative, which is the
//! format's own way of saying angstroms. Both are normalised to angstroms here.

use super::{
    Axes, BOHR_ANGSTROMS, GridError, canonical_value_index, density_map, value_count, zeroed_values,
};
use crate::mrc::DensityMap;

const FORMAT: &str = "cube";

/// One nucleus recorded in a cube's header, in angstroms.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CubeAtom {
    /// Atomic number.
    pub number: i32,
    /// Nuclear charge as written, which for an ECP is not the atomic number.
    pub charge: f64,
    /// Cartesian position in angstroms.
    pub position: [f64; 3],
}

/// A cube file's grid, the nuclei it was computed around, and its orbitals.
#[derive(Clone, Debug, PartialEq)]
pub struct CubeGrid {
    /// The first (or only) scalar field, canonicalised to X-fastest.
    pub map: DensityMap,
    /// Nuclei from the header.
    pub atoms: Vec<CubeAtom>,
    /// Molecular-orbital indices, empty for an ordinary density cube.
    pub orbitals: Vec<i32>,
    /// Fields stored per voxel; greater than one for a multi-orbital cube.
    pub fields: usize,
}

impl CubeGrid {
    /// Extracts one field of a multi-orbital cube as its own map.
    ///
    /// Values interleave per voxel in the file, so a single field is a strided
    /// gather; the result is contiguous and indistinguishable from a
    /// single-field cube. Returns `None` past the last field.
    #[must_use]
    pub fn field(&self, index: usize) -> Option<DensityMap> {
        if index >= self.fields {
            return None;
        }
        if self.fields == 1 {
            return Some(self.map.clone());
        }
        let values = self
            .map
            .values
            .iter()
            .skip(index)
            .step_by(self.fields)
            .copied()
            .collect();
        Some(DensityMap {
            values,
            ..self.map.clone()
        })
    }
}

/// Reads a Gaussian cube.
///
/// # Errors
///
/// Returns [`GridError`] for non-UTF-8 text, a malformed or truncated header,
/// dimensions out of range, or a value that is not a number.
pub fn read_cube(bytes: &[u8]) -> Result<CubeGrid, GridError> {
    let Ok(text) = std::str::from_utf8(bytes) else {
        return Err(GridError::NotUtf8 { format: FORMAT });
    };
    let mut lines = text.lines();
    let header = read_header(&mut lines)?;

    let voxels = value_count(header.axes.counts, FORMAT)?;
    let Some(total) = voxels.checked_mul(header.fields) else {
        return Err(GridError::SizeOverflow { format: FORMAT });
    };
    if total > super::MAX_GRID_VALUES {
        return Err(GridError::SizeOverflow { format: FORMAT });
    }

    let mut values = zeroed_values(total, FORMAT)?;
    let mut parsed = 0_usize;
    for token in lines.flat_map(str::split_ascii_whitespace) {
        if parsed == total {
            break;
        }
        match token.parse::<f32>() {
            Ok(value) => {
                let destination = canonical_value_index(parsed, header.axes.counts, header.fields);
                values[destination] = value;
                parsed += 1;
            }
            Err(_) => {
                return Err(GridError::InvalidNumber {
                    format: FORMAT,
                    index: parsed,
                });
            }
        }
    }
    if parsed != total {
        return Err(GridError::Truncated { format: FORMAT });
    }
    Ok(CubeGrid {
        map: density_map(&header.axes, values),
        atoms: header.atoms,
        orbitals: header.orbitals,
        fields: header.fields,
    })
}

/// Everything above a cube's value block.
struct CubeHeader {
    axes: Axes,
    atoms: Vec<CubeAtom>,
    orbitals: Vec<i32>,
    fields: usize,
}

/// Reads the comment lines, the grid geometry, the nuclei and, for an orbital
/// cube, the orbital indices — leaving the iterator on the first value.
fn read_header<'a>(lines: &mut impl Iterator<Item = &'a str>) -> Result<CubeHeader, GridError> {
    // Two free-form comment lines open every cube.
    for _ in 0..2 {
        if lines.next().is_none() {
            return Err(GridError::Truncated { format: FORMAT });
        }
    }

    let Some(counts_line) = lines.next() else {
        return Err(GridError::Truncated { format: FORMAT });
    };
    let mut counts_fields = counts_line.split_ascii_whitespace();
    let atom_count = read_i32(&mut counts_fields, "atom count")?;
    let mut origin = [0.0_f64; 3];
    for value in &mut origin {
        *value = read_f64(&mut counts_fields, "grid origin")?;
    }
    // A cube written by some programs adds a field count here; when it is
    // absent the cube carries exactly one field.
    let declared_fields = counts_fields
        .next()
        .and_then(|token| token.parse::<i32>().ok())
        .filter(|count| *count > 0);

    let axes = read_axes(lines)?;
    for value in &mut origin {
        *value *= axes.scale;
    }

    let atoms = read_atoms(lines, atom_count, axes.scale)?;
    let (orbitals, fields) = read_orbitals(lines, atom_count, declared_fields)?;
    Ok(CubeHeader {
        axes: Axes {
            origin,
            steps: axes.steps,
            counts: axes.counts,
        },
        atoms,
        orbitals,
        fields,
    })
}

/// Voxel counts, step vectors already scaled to angstroms, and the scale that
/// was applied — the nuclei need the same scale, so it is reported rather than
/// recomputed.
struct AxisBlock {
    counts: [usize; 3],
    steps: [[f64; 3]; 3],
    scale: f64,
}

/// The three axis lines.
fn read_axes<'a>(lines: &mut impl Iterator<Item = &'a str>) -> Result<AxisBlock, GridError> {
    let mut counts = [0_usize; 3];
    let mut steps = [[0.0_f64; 3]; 3];
    let mut angstroms = false;
    for axis in 0..3 {
        let Some(line) = lines.next() else {
            return Err(GridError::Truncated { format: FORMAT });
        };
        let mut fields = line.split_ascii_whitespace();
        let declared = read_i32(&mut fields, "axis voxel count")?;
        // The format signals its length unit through the sign of the count.
        if declared < 0 {
            angstroms = true;
        }
        let Ok(count) = usize::try_from(declared.unsigned_abs()) else {
            return Err(GridError::SizeOverflow { format: FORMAT });
        };
        counts[axis] = count;
        for component in &mut steps[axis] {
            *component = read_f64(&mut fields, "axis step vector")?;
        }
    }
    let scale = if angstroms { 1.0 } else { BOHR_ANGSTROMS };
    for step in &mut steps {
        for component in step {
            *component *= scale;
        }
    }
    Ok(AxisBlock {
        counts,
        steps,
        scale,
    })
}

fn read_atoms<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    atom_count: i32,
    scale: f64,
) -> Result<Vec<CubeAtom>, GridError> {
    let Ok(nuclei) = usize::try_from(atom_count.unsigned_abs()) else {
        return Err(GridError::SizeOverflow { format: FORMAT });
    };
    let mut atoms = Vec::new();
    if atoms.try_reserve_exact(nuclei).is_err() {
        return Err(GridError::ResourceLimit { format: FORMAT });
    }
    for _ in 0..nuclei {
        let Some(line) = lines.next() else {
            return Err(GridError::Truncated { format: FORMAT });
        };
        let mut fields = line.split_ascii_whitespace();
        let number = read_i32(&mut fields, "atomic number")?;
        let charge = read_f64(&mut fields, "nuclear charge")?;
        let mut position = [0.0_f64; 3];
        for value in &mut position {
            *value = read_f64(&mut fields, "nuclear position")? * scale;
        }
        atoms.push(CubeAtom {
            number,
            charge,
            position,
        });
    }
    Ok(atoms)
}

/// A negative atom count marks an orbital cube, whose extra line names the
/// orbitals and therefore how many values each voxel carries.
fn read_orbitals<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    atom_count: i32,
    declared_fields: Option<i32>,
) -> Result<(Vec<i32>, usize), GridError> {
    let fallback = match declared_fields {
        Some(fields) => fields,
        None => 1,
    };
    if atom_count >= 0 {
        let Ok(fields) = usize::try_from(fallback.max(1)) else {
            return Err(GridError::SizeOverflow { format: FORMAT });
        };
        return Ok((Vec::new(), fields));
    }
    let Some(line) = lines.next() else {
        return Err(GridError::Truncated { format: FORMAT });
    };
    let mut tokens = line.split_ascii_whitespace();
    let declared = read_i32(&mut tokens, "orbital count")?;
    let Ok(count) = usize::try_from(declared.max(0)) else {
        return Err(GridError::SizeOverflow { format: FORMAT });
    };
    let mut orbitals = Vec::new();
    for token in tokens.take(count) {
        match token.parse::<i32>() {
            Ok(index) => orbitals.push(index),
            Err(_) => {
                return Err(GridError::InvalidHeader {
                    format: FORMAT,
                    reason: "orbital index",
                });
            }
        }
    }
    let Ok(fields) = usize::try_from(declared.max(1)) else {
        return Err(GridError::SizeOverflow { format: FORMAT });
    };
    Ok((orbitals, fields))
}

fn read_i32<'a>(
    fields: &mut impl Iterator<Item = &'a str>,
    reason: &'static str,
) -> Result<i32, GridError> {
    let Some(token) = fields.next() else {
        return Err(GridError::InvalidHeader {
            format: FORMAT,
            reason,
        });
    };
    match token.parse::<i32>() {
        Ok(value) => Ok(value),
        // Some writers emit a whole-number count as a float.
        Err(_) => match token.parse::<f64>() {
            Ok(value) => Ok(crate::numeric::f64_to_i32(value)),
            Err(_) => Err(GridError::InvalidHeader {
                format: FORMAT,
                reason,
            }),
        },
    }
}

fn read_f64<'a>(
    fields: &mut impl Iterator<Item = &'a str>,
    reason: &'static str,
) -> Result<f64, GridError> {
    let Some(token) = fields.next() else {
        return Err(GridError::InvalidHeader {
            format: FORMAT,
            reason,
        });
    };
    match token.parse::<f64>() {
        Ok(value) => Ok(value),
        Err(_) => Err(GridError::InvalidHeader {
            format: FORMAT,
            reason,
        }),
    }
}
