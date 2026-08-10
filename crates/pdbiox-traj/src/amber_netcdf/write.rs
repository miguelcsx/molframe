//! AMBER `NetCDF` encoding.

use netcdf_writer::{NcAttrValue, NcFileBuilder, NcWriteOptions};

use crate::Timestep;
use crate::numeric::f32_from_f64;

use super::AmberNetcdfError;
use super::schema::{
    ANGLE_UNIT, ATOM, CELL_ANGLES, CELL_ANGULAR, CELL_LENGTHS, CELL_SPATIAL, CONVENTION,
    CONVENTION_ATTRIBUTE, COORDINATES, FORCE_UNIT, FORCES, FRAME, KILOCALORIE_TO_KILOJOULE, LABEL,
    LENGTH_UNIT, PROGRAM_ATTRIBUTE, PROGRAM_VERSION_ATTRIBUTE, SPATIAL, TIME, TIME_UNIT,
    UNITS_ATTRIBUTE, VELOCITIES, VELOCITY_UNIT, VERSION, VERSION_ATTRIBUTE,
};

/// Numeric precision used for physical arrays in AMBER `NetCDF` output.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AmberNetcdfPrecision {
    /// IEEE-754 binary32, the conventional compact representation.
    #[default]
    Single,
    /// IEEE-754 binary64.
    Double,
}

/// Explicit AMBER `NetCDF` writer configuration.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AmberNetcdfWriteOptions {
    /// Numeric precision for coordinates and optional vector streams.
    pub precision: AmberNetcdfPrecision,
    /// Program version recorded in the file metadata.
    pub program: String,
    /// Program version recorded in the file metadata.
    pub program_version: String,
}

impl Default for AmberNetcdfWriteOptions {
    fn default() -> Self {
        Self {
            precision: AmberNetcdfPrecision::Single,
            program: env!("CARGO_PKG_NAME").into(),
            program_version: env!("CARGO_PKG_VERSION").into(),
        }
    }
}

/// Encodes frames as a canonical AMBER `NetCDF` classic file.
///
/// # Errors
///
/// Returns an error when frames disagree in atom count or optional stream
/// presence, contain invalid numbers/cells, or cannot be represented by the
/// `NetCDF` encoder.
pub fn write_amber_netcdf(
    frames: &[Timestep],
    options: AmberNetcdfWriteOptions,
) -> Result<Vec<u8>, AmberNetcdfError> {
    let schema = FrameSchema::validate(frames)?;
    match options.precision {
        AmberNetcdfPrecision::Single => encode::<f32>(frames, schema, options),
        AmberNetcdfPrecision::Double => encode::<f64>(frames, schema, options),
    }
}

#[derive(Clone, Copy)]
struct FrameSchema {
    atoms: usize,
    time: Presence,
    velocities: Presence,
    forces: Presence,
    cell: Presence,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Presence {
    Absent,
    Present,
}

impl Presence {
    const fn of(value: bool) -> Self {
        if value { Self::Present } else { Self::Absent }
    }

    const fn is_present(self) -> bool {
        matches!(self, Self::Present)
    }
}

impl FrameSchema {
    fn validate(frames: &[Timestep]) -> Result<Self, AmberNetcdfError> {
        let first = frames.first().ok_or(AmberNetcdfError::InconsistentFrames)?;
        if first.positions.is_empty() {
            return Err(AmberNetcdfError::InconsistentFrames);
        }
        let schema = Self {
            atoms: first.positions.len(),
            time: Presence::of(first.time.is_some()),
            velocities: Presence::of(first.velocities.is_some()),
            forces: Presence::of(first.forces.is_some()),
            cell: Presence::of(first.cell.is_some()),
        };
        for frame in frames {
            if frame.positions.len() != schema.atoms
                || Presence::of(frame.time.is_some()) != schema.time
                || !matches_vectors(frame.velocities.as_deref(), schema.velocities, schema.atoms)
                || !matches_vectors(frame.forces.as_deref(), schema.forces, schema.atoms)
                || Presence::of(frame.cell.is_some()) != schema.cell
            {
                return Err(AmberNetcdfError::InconsistentFrames);
            }
            validate_frame(frame)?;
        }
        Ok(schema)
    }
}

fn encode<T: Physical>(
    frames: &[Timestep],
    schema: FrameSchema,
    options: AmberNetcdfWriteOptions,
) -> Result<Vec<u8>, AmberNetcdfError> {
    let mut builder = NcFileBuilder::new();
    builder.add_attribute(CONVENTION_ATTRIBUTE, NcAttrValue::Chars(CONVENTION.into()))?;
    builder.add_attribute(VERSION_ATTRIBUTE, NcAttrValue::Chars(VERSION.into()))?;
    if options.program.is_empty() || options.program_version.is_empty() {
        return Err(AmberNetcdfError::InvalidValue);
    }
    builder.add_attribute(PROGRAM_ATTRIBUTE, NcAttrValue::Chars(options.program))?;
    builder.add_attribute(
        PROGRAM_VERSION_ATTRIBUTE,
        NcAttrValue::Chars(options.program_version),
    )?;
    let frame = builder.add_unlimited_dimension(FRAME)?;
    let atom = builder.add_dimension(ATOM, to_u64(schema.atoms)?)?;
    let spatial = builder.add_dimension(SPATIAL, 3)?;
    let cell_spatial = builder.add_dimension(CELL_SPATIAL, 3)?;
    let cell_angular = builder.add_dimension(CELL_ANGULAR, 3)?;
    let label = builder.add_dimension(LABEL, 5)?;
    let spatial_labels = builder.add_char_variable(SPATIAL, &[spatial])?;
    let cell_spatial_labels = builder.add_char_variable(CELL_SPATIAL, &[cell_spatial])?;
    let cell_angular_labels = builder.add_char_variable(CELL_ANGULAR, &[cell_angular, label])?;
    builder.write_char_variable(spatial_labels, b"xyz")?;
    builder.write_char_variable(cell_spatial_labels, b"abc")?;
    builder.write_char_variable(cell_angular_labels, b"alphabeta gamma")?;
    let variables = add_physical_variables::<T>(
        &mut builder,
        schema,
        frame,
        atom,
        spatial,
        cell_spatial,
        cell_angular,
    )?;
    write_physical::<T>(&mut builder, frames, variables)?;
    let (_, bytes) = builder.to_vec(NcWriteOptions::offset64())?;
    Ok(bytes)
}

#[derive(Clone, Copy)]
struct PhysicalVariables {
    coordinates: netcdf_writer::VariableId,
    time: Option<netcdf_writer::VariableId>,
    velocities: Option<netcdf_writer::VariableId>,
    forces: Option<netcdf_writer::VariableId>,
    cell_lengths: Option<netcdf_writer::VariableId>,
    cell_angles: Option<netcdf_writer::VariableId>,
}

fn add_physical_variables<T: netcdf_writer::NcWriteType>(
    builder: &mut NcFileBuilder,
    schema: FrameSchema,
    frame: netcdf_writer::DimensionId,
    atom: netcdf_writer::DimensionId,
    spatial: netcdf_writer::DimensionId,
    cell_spatial: netcdf_writer::DimensionId,
    cell_angular: netcdf_writer::DimensionId,
) -> Result<PhysicalVariables, AmberNetcdfError> {
    let coordinates =
        add_variable::<T>(builder, COORDINATES, &[frame, atom, spatial], LENGTH_UNIT)?;
    let time = schema
        .time
        .is_present()
        .then(|| add_variable::<T>(builder, TIME, &[frame], TIME_UNIT))
        .transpose()?;
    let velocities = schema
        .velocities
        .is_present()
        .then(|| add_variable::<T>(builder, VELOCITIES, &[frame, atom, spatial], VELOCITY_UNIT))
        .transpose()?;
    let forces = schema
        .forces
        .is_present()
        .then(|| add_variable::<T>(builder, FORCES, &[frame, atom, spatial], FORCE_UNIT))
        .transpose()?;
    let cell_lengths = schema
        .cell
        .is_present()
        .then(|| add_variable::<T>(builder, CELL_LENGTHS, &[frame, cell_spatial], LENGTH_UNIT))
        .transpose()?;
    let cell_angles = schema
        .cell
        .is_present()
        .then(|| add_variable::<T>(builder, CELL_ANGLES, &[frame, cell_angular], ANGLE_UNIT))
        .transpose()?;
    Ok(PhysicalVariables {
        coordinates,
        time,
        velocities,
        forces,
        cell_lengths,
        cell_angles,
    })
}

fn add_variable<T: netcdf_writer::NcWriteType>(
    builder: &mut NcFileBuilder,
    name: &str,
    dimensions: &[netcdf_writer::DimensionId],
    units: &str,
) -> Result<netcdf_writer::VariableId, AmberNetcdfError> {
    let variable = builder.add_variable::<T>(name, dimensions)?;
    builder.add_variable_attribute(variable, UNITS_ATTRIBUTE, NcAttrValue::Chars(units.into()))?;
    Ok(variable)
}

trait Physical: netcdf_writer::NcWriteType + Copy {
    fn convert(value: f64) -> Result<Self, AmberNetcdfError>;
}

impl Physical for f64 {
    fn convert(value: f64) -> Result<Self, AmberNetcdfError> {
        Ok(value)
    }
}

impl Physical for f32 {
    fn convert(value: f64) -> Result<Self, AmberNetcdfError> {
        f32_from_f64(value).ok_or(AmberNetcdfError::InvalidValue)
    }
}

fn write_physical<T: Physical>(
    builder: &mut NcFileBuilder,
    frames: &[Timestep],
    variables: PhysicalVariables,
) -> Result<(), AmberNetcdfError> {
    write_values(
        builder,
        variables.coordinates,
        flatten_positions(frames),
        T::convert,
    )?;
    if let Some(variable) = variables.time {
        write_values(
            builder,
            variable,
            frames.iter().filter_map(|frame| frame.time),
            T::convert,
        )?;
    }
    if let Some(variable) = variables.velocities {
        write_values(
            builder,
            variable,
            flatten_optional(frames, false, 1.0),
            T::convert,
        )?;
    }
    if let Some(variable) = variables.forces {
        write_values(
            builder,
            variable,
            flatten_optional(frames, true, KILOCALORIE_TO_KILOJOULE),
            T::convert,
        )?;
    }
    if let Some(variable) = variables.cell_lengths {
        write_values(builder, variable, flatten_cell(frames, false), T::convert)?;
    }
    if let Some(variable) = variables.cell_angles {
        write_values(builder, variable, flatten_cell(frames, true), T::convert)?;
    }
    Ok(())
}

fn write_values<T: Physical>(
    builder: &mut NcFileBuilder,
    variable: netcdf_writer::VariableId,
    values: impl Iterator<Item = f64>,
    convert: fn(f64) -> Result<T, AmberNetcdfError>,
) -> Result<(), AmberNetcdfError> {
    let values = values.map(convert).collect::<Result<Vec<_>, _>>()?;
    builder.write_variable(variable, &values)?;
    Ok(())
}

fn flatten_positions(frames: &[Timestep]) -> impl Iterator<Item = f64> + '_ {
    frames
        .iter()
        .flat_map(|frame| frame.positions.iter())
        .flat_map(|value| value.map(f64::from))
}

fn flatten_optional(
    frames: &[Timestep],
    forces: bool,
    divisor: f64,
) -> impl Iterator<Item = f64> + '_ {
    frames
        .iter()
        .filter_map(move |frame| {
            if forces {
                frame.forces.as_ref()
            } else {
                frame.velocities.as_ref()
            }
        })
        .flat_map(|values| values.iter())
        .flat_map(|value| value.map(f64::from))
        .map(move |value| value / divisor)
}

fn flatten_cell(frames: &[Timestep], angles: bool) -> impl Iterator<Item = f64> + '_ {
    frames
        .iter()
        .filter_map(move |frame| frame.cell)
        .flat_map(move |cell| if angles { cell.angles } else { cell.lengths })
}

fn validate_frame(frame: &Timestep) -> Result<(), AmberNetcdfError> {
    let finite_positions = frame
        .positions
        .iter()
        .flatten()
        .all(|value| value.is_finite());
    let finite_velocities = frame
        .velocities
        .iter()
        .flatten()
        .flatten()
        .all(|value| value.is_finite());
    let finite_forces = frame
        .forces
        .iter()
        .flatten()
        .flatten()
        .all(|value| value.is_finite());
    let finite_time = frame.time.is_none_or(f64::is_finite);
    let valid_cell = frame.cell.is_none_or(|cell| {
        cell.lengths
            .iter()
            .all(|value| value.is_finite() && *value > 0.0)
            && cell
                .angles
                .iter()
                .all(|value| value.is_finite() && *value > 0.0 && *value < 180.0)
    });
    if finite_positions && finite_velocities && finite_forces && finite_time && valid_cell {
        Ok(())
    } else {
        Err(AmberNetcdfError::InvalidValue)
    }
}

fn matches_vectors(values: Option<&[[f32; 3]]>, present: Presence, atoms: usize) -> bool {
    Presence::of(values.is_some_and(|values| values.len() == atoms)) == present
}

fn to_u64(value: usize) -> Result<u64, AmberNetcdfError> {
    u64::try_from(value).map_err(|_| AmberNetcdfError::InvalidValue)
}
