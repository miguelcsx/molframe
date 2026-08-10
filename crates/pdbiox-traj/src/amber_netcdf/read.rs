//! AMBER `NetCDF` decoding.

use crate::Timestep;
use crate::numeric::f32_from_f64;
use pdbiox_core::structure::UnitCell;

use super::AmberNetcdfError;
use super::AmberNetcdfPrecision;
use super::schema::{
    ANGLE_UNIT, CELL_ANGLES, CELL_LENGTHS, COORDINATES, FORCE_UNIT, FORCES,
    KILOCALORIE_TO_KILOJOULE, LENGTH_UNIT, PROGRAM_ATTRIBUTE, PROGRAM_VERSION_ATTRIBUTE, TIME,
    TIME_UNIT, VELOCITIES, VELOCITY_UNIT, optional_attribute_string, require_variable,
    validate_convention, validate_units,
};

/// Complete decoded AMBER trajectory and source encoding metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct AmberNetcdfTrajectory {
    /// Canonical trajectory frames.
    pub frames: Vec<Timestep>,
    /// Source precision and producer attributes, when declared.
    pub metadata: AmberNetcdfMetadata,
}

/// Source metadata required to reproduce the AMBER container faithfully.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AmberNetcdfMetadata {
    /// Numeric representation of coordinate and vector arrays.
    pub precision: AmberNetcdfPrecision,
    /// Producer name, if the source declares it.
    pub program: Option<Box<str>>,
    /// Producer version, if the source declares it.
    pub program_version: Option<Box<str>>,
}

/// Parses an AMBER `NetCDF` trajectory from memory.
///
/// Coordinates are returned in ångström, time in picoseconds, velocities in
/// ångström/ps, and forces in kJ mol⁻¹ Å⁻¹. Unit metadata is mandatory for
/// every present physical quantity; no implicit-unit compatibility path is
/// applied.
///
/// # Errors
///
/// Returns an error for invalid convention metadata, dimensions, units,
/// non-finite values, or inconsistent optional streams.
pub fn parse_amber_netcdf(bytes: &[u8]) -> Result<Vec<Timestep>, AmberNetcdfError> {
    parse_amber_netcdf_record(bytes).map(|trajectory| trajectory.frames)
}

/// Parses frames together with source precision and producer attributes.
///
/// # Errors
///
/// Returns the same strict convention, shape, unit and numeric errors as
/// [`parse_amber_netcdf`].
pub fn parse_amber_netcdf_record(bytes: &[u8]) -> Result<AmberNetcdfTrajectory, AmberNetcdfError> {
    let file = netcdf_reader::NcFile::from_bytes(bytes)?;
    validate_convention(&file)?;
    let coordinates = require_variable(&file, COORDINATES)?;
    validate_units(coordinates, COORDINATES, LENGTH_UNIT)?;
    validate_dimensions(coordinates, COORDINATES, &["frame", "atom", "spatial"])?;
    let shape = coordinates.shape();
    let (frame_count, atom_count) = coordinate_shape(&shape)?;
    let coordinate_values = read_finite(&file, COORDINATES)?;

    let times = read_optional(&file, TIME, &[frame_count], TIME_UNIT, 1.0)?;
    let velocities = read_optional(
        &file,
        VELOCITIES,
        &[frame_count, atom_count, 3],
        VELOCITY_UNIT,
        1.0,
    )?;
    let forces = read_optional(
        &file,
        FORCES,
        &[frame_count, atom_count, 3],
        FORCE_UNIT,
        KILOCALORIE_TO_KILOJOULE,
    )?;
    let lengths = read_optional(&file, CELL_LENGTHS, &[frame_count, 3], LENGTH_UNIT, 1.0)?;
    let angles = read_optional(&file, CELL_ANGLES, &[frame_count, 3], ANGLE_UNIT, 1.0)?;
    if lengths.is_some() != angles.is_some() {
        return Err(AmberNetcdfError::InconsistentFrames);
    }
    let precision = match coordinates.dtype {
        netcdf_reader::NcType::Float => AmberNetcdfPrecision::Single,
        netcdf_reader::NcType::Double => AmberNetcdfPrecision::Double,
        _ => return Err(AmberNetcdfError::InvalidValue),
    };
    let frames = assemble_frames(&DecodedArrays {
        frame_count,
        atom_count,
        coordinates: &coordinate_values,
        times: times.as_deref(),
        velocities: velocities.as_deref(),
        forces: forces.as_deref(),
        lengths: lengths.as_deref(),
        angles: angles.as_deref(),
    })?;
    Ok(AmberNetcdfTrajectory {
        frames,
        metadata: AmberNetcdfMetadata {
            precision,
            program: optional_attribute_string(&file, PROGRAM_ATTRIBUTE)
                .map(String::into_boxed_str),
            program_version: optional_attribute_string(&file, PROGRAM_VERSION_ATTRIBUTE)
                .map(String::into_boxed_str),
        },
    })
}

fn coordinate_shape(shape: &[u64]) -> Result<(usize, usize), AmberNetcdfError> {
    let [frames, atoms, spatial] = shape else {
        return Err(AmberNetcdfError::InvalidShape {
            variable: COORDINATES,
            expected: "[frame, atom, 3]",
            found: shape.to_vec(),
        });
    };
    if *frames == 0 || *atoms == 0 || *spatial != 3 {
        return Err(AmberNetcdfError::InvalidShape {
            variable: COORDINATES,
            expected: "non-empty [frame, atom, 3]",
            found: shape.to_vec(),
        });
    }
    Ok((to_usize(*frames)?, to_usize(*atoms)?))
}

fn read_optional(
    file: &netcdf_reader::NcFile,
    name: &'static str,
    expected: &[usize],
    units: &'static str,
    scale: f64,
) -> Result<Option<Vec<f64>>, AmberNetcdfError> {
    let Ok(variable) = file.variable(name) else {
        return Ok(None);
    };
    validate_units(variable, name, units)?;
    validate_dimensions(variable, name, expected_dimensions(name))?;
    let found = variable.shape();
    let expected_u64 = expected
        .iter()
        .map(|value| u64::try_from(*value).map_err(|_| AmberNetcdfError::InvalidValue))
        .collect::<Result<Vec<_>, _>>()?;
    if found != expected_u64 {
        return Err(AmberNetcdfError::InvalidShape {
            variable: name,
            expected: expected_description(name),
            found,
        });
    }
    let mut values = read_finite(file, name)?;
    for value in &mut values {
        *value *= scale;
        if !value.is_finite() {
            return Err(AmberNetcdfError::InvalidValue);
        }
    }
    Ok(Some(values))
}

fn expected_description(name: &str) -> &'static str {
    match name {
        TIME => "[frame]",
        CELL_LENGTHS | CELL_ANGLES => "[frame, 3]",
        _ => "[frame, atom, 3]",
    }
}

fn read_finite(
    file: &netcdf_reader::NcFile,
    name: &'static str,
) -> Result<Vec<f64>, AmberNetcdfError> {
    let array = file.read_variable_unpacked(name)?;
    let values: Vec<_> = array.iter().copied().collect();
    if values.iter().all(|value| value.is_finite()) {
        Ok(values)
    } else {
        Err(AmberNetcdfError::InvalidValue)
    }
}

fn expected_dimensions(name: &str) -> &'static [&'static str] {
    match name {
        TIME => &["frame"],
        CELL_LENGTHS => &["frame", "cell_spatial"],
        CELL_ANGLES => &["frame", "cell_angular"],
        _ => &["frame", "atom", "spatial"],
    }
}

fn validate_dimensions(
    variable: &netcdf_reader::NcVariable,
    name: &'static str,
    expected: &[&str],
) -> Result<(), AmberNetcdfError> {
    let matches = variable
        .dimensions()
        .iter()
        .map(|dimension| dimension.name.as_str())
        .eq(expected.iter().copied());
    if matches {
        Ok(())
    } else {
        Err(AmberNetcdfError::InvalidShape {
            variable: name,
            expected: expected_description(name),
            found: variable.shape(),
        })
    }
}

struct DecodedArrays<'a> {
    frame_count: usize,
    atom_count: usize,
    coordinates: &'a [f64],
    times: Option<&'a [f64]>,
    velocities: Option<&'a [f64]>,
    forces: Option<&'a [f64]>,
    lengths: Option<&'a [f64]>,
    angles: Option<&'a [f64]>,
}

fn assemble_frames(arrays: &DecodedArrays<'_>) -> Result<Vec<Timestep>, AmberNetcdfError> {
    let mut frames = Vec::with_capacity(arrays.frame_count);
    for frame in 0..arrays.frame_count {
        let vector_start = frame
            .checked_mul(arrays.atom_count)
            .and_then(|value| value.checked_mul(3))
            .ok_or(AmberNetcdfError::InvalidValue)?;
        let positions = vectors(arrays.coordinates, vector_start, arrays.atom_count)?;
        let velocities = arrays
            .velocities
            .map(|values| vectors(values, vector_start, arrays.atom_count))
            .transpose()?;
        let forces = arrays
            .forces
            .map(|values| vectors(values, vector_start, arrays.atom_count))
            .transpose()?;
        let cell = match (arrays.lengths, arrays.angles) {
            (Some(lengths), Some(angles)) => Some(UnitCell {
                lengths: triple(lengths, frame * 3)?,
                angles: triple(angles, frame * 3)?,
            }),
            _ => None,
        };
        frames.push(Timestep {
            frame,
            time: arrays.times.and_then(|values| values.get(frame)).copied(),
            positions,
            velocities,
            forces,
            cell,
            ..Timestep::default()
        });
    }
    populate_dt(&mut frames);
    Ok(frames)
}

fn vectors(values: &[f64], start: usize, count: usize) -> Result<Vec<[f32; 3]>, AmberNetcdfError> {
    let end = start
        .checked_add(count.checked_mul(3).ok_or(AmberNetcdfError::InvalidValue)?)
        .ok_or(AmberNetcdfError::InvalidValue)?;
    values
        .get(start..end)
        .ok_or(AmberNetcdfError::InvalidValue)?
        .chunks_exact(3)
        .map(|chunk| Ok([to_f32(chunk[0])?, to_f32(chunk[1])?, to_f32(chunk[2])?]))
        .collect()
}

fn triple(values: &[f64], start: usize) -> Result<[f64; 3], AmberNetcdfError> {
    let slice = values
        .get(start..start.checked_add(3).ok_or(AmberNetcdfError::InvalidValue)?)
        .ok_or(AmberNetcdfError::InvalidValue)?;
    Ok([slice[0], slice[1], slice[2]])
}

fn to_f32(value: f64) -> Result<f32, AmberNetcdfError> {
    f32_from_f64(value).ok_or(AmberNetcdfError::InvalidValue)
}

fn to_usize(value: u64) -> Result<usize, AmberNetcdfError> {
    usize::try_from(value).map_err(|_| AmberNetcdfError::InvalidValue)
}

fn populate_dt(frames: &mut [Timestep]) {
    for index in 1..frames.len() {
        if let (Some(previous), Some(current)) = (frames[index - 1].time, frames[index].time) {
            frames[index].dt = Some(current - previous);
        }
    }
}
