//! H5MD encoding.

use hdf5_writer::{AttributeBuilder, DatasetBuilder, Hdf5Builder, WriteOptions};

use crate::Timestep;
use crate::cell::vectors_from_cell;

use super::schema::{
    BOUNDARY_ATTRIBUTE, CREATOR_GROUP, CREATOR_NAME_ATTRIBUTE, CREATOR_VERSION_ATTRIBUTE,
    DIMENSION_ATTRIBUTE, FORCE, H5MD_GROUP, H5MD_VERSION, PERIODIC, POSITION, STEP, TIME,
    UNIT_ATTRIBUTE, VALUE, VELOCITY, VERSION_ATTRIBUTE, box_group_path, box_path, element_path,
};
use super::{H5mdError, H5mdMetadata, H5mdOptions};

/// Encodes frames as an H5MD 1.1 file.
///
/// Every physical dataset receives an explicit unit attribute from `options`.
/// Values are converted from molframe canonical units using the supplied scales.
///
/// # Errors
///
/// Returns an error for invalid options, empty/inconsistent frames, absent
/// times, invalid cells/numbers, or HDF5 encoding failures.
pub fn write_h5md(frames: &[Timestep], options: &H5mdOptions) -> Result<Vec<u8>, H5mdError> {
    write_h5md_with_metadata(
        frames,
        &H5mdMetadata {
            options: options.clone(),
            creator_name: env!("CARGO_PKG_NAME").into(),
            creator_version: env!("CARGO_PKG_VERSION").into(),
        },
    )
}

/// Encodes frames while retaining an explicitly supplied H5MD creator identity.
///
/// # Errors
///
/// Returns the same strict schema and numeric errors as [`write_h5md`], and
/// refuses empty creator fields.
pub fn write_h5md_with_metadata(
    frames: &[Timestep],
    metadata: &H5mdMetadata,
) -> Result<Vec<u8>, H5mdError> {
    let options = &metadata.options;
    options.validate()?;
    if metadata.creator_name.is_empty() || metadata.creator_version.is_empty() {
        return Err(H5mdError::InvalidMetadata);
    }
    let schema = FrameSchema::validate(frames)?;
    let steps = collect_steps(frames)?;
    let times = collect_times(frames, options.units.time_to_picosecond)?;
    let mut builder = metadata_builder(metadata)?;
    builder = add_element(
        builder,
        ElementWriteSpec {
            group: &options.particle_group,
            element: POSITION,
            values: flatten_positions(frames, options.units.length_to_angstrom)?,
            shape: schema.shape()?,
            steps: &steps,
            times: &times,
            value_unit: &options.units.length_unit,
            time_unit: &options.units.time_unit,
        },
    )?;
    if schema.velocities {
        builder = add_element(
            builder,
            ElementWriteSpec {
                group: &options.particle_group,
                element: VELOCITY,
                values: flatten_optional(
                    frames,
                    false,
                    options.units.velocity_to_angstrom_per_picosecond,
                )?,
                shape: schema.shape()?,
                steps: &steps,
                times: &times,
                value_unit: &options.units.velocity_unit,
                time_unit: &options.units.time_unit,
            },
        )?;
    }
    if schema.forces {
        builder = add_element(
            builder,
            ElementWriteSpec {
                group: &options.particle_group,
                element: FORCE,
                values: flatten_optional(
                    frames,
                    true,
                    options.units.force_to_kilojoule_per_mole_angstrom,
                )?,
                shape: schema.shape()?,
                steps: &steps,
                times: &times,
                value_unit: &options.units.force_unit,
                time_unit: &options.units.time_unit,
            },
        )?;
    }
    if schema.cell {
        builder = add_box(builder, frames, options, &steps, &times)?;
    }
    builder
        .into_plan()?
        .encode(WriteOptions::default())
        .map_err(Into::into)
}

#[derive(Clone, Copy)]
struct FrameSchema {
    frames: usize,
    atoms: usize,
    velocities: bool,
    forces: bool,
    cell: bool,
}

impl FrameSchema {
    fn validate(frames: &[Timestep]) -> Result<Self, H5mdError> {
        let first = frames.first().ok_or(H5mdError::InconsistentFrames)?;
        if first.positions.is_empty() || first.time.is_none() {
            return Err(H5mdError::InconsistentFrames);
        }
        let schema = Self {
            frames: frames.len(),
            atoms: first.positions.len(),
            velocities: first.velocities.is_some(),
            forces: first.forces.is_some(),
            cell: first.cell.is_some(),
        };
        for frame in frames {
            if frame.positions.len() != schema.atoms
                || frame.time.is_none()
                || !matches_vectors(frame.velocities.as_deref(), schema.velocities, schema.atoms)
                || !matches_vectors(frame.forces.as_deref(), schema.forces, schema.atoms)
                || frame.cell.is_some() != schema.cell
            {
                return Err(H5mdError::InconsistentFrames);
            }
            validate_frame(frame)?;
        }
        Ok(schema)
    }

    fn shape(self) -> Result<Vec<u64>, H5mdError> {
        Ok(vec![to_u64(self.frames)?, to_u64(self.atoms)?, 3])
    }
}

fn metadata_builder(metadata: &H5mdMetadata) -> Result<Hdf5Builder, H5mdError> {
    Ok(Hdf5Builder::new()
        .group_attribute(
            H5MD_GROUP,
            AttributeBuilder::vector(VERSION_ATTRIBUTE, &H5MD_VERSION)?,
        )
        .group_attribute(
            CREATOR_GROUP,
            AttributeBuilder::fixed_string(CREATOR_NAME_ATTRIBUTE, &metadata.creator_name),
        )
        .group_attribute(
            CREATOR_GROUP,
            AttributeBuilder::fixed_string(CREATOR_VERSION_ATTRIBUTE, &metadata.creator_version),
        ))
}

struct ElementWriteSpec<'a> {
    group: &'a str,
    element: &'a str,
    values: Vec<f64>,
    shape: Vec<u64>,
    steps: &'a [i64],
    times: &'a [f64],
    value_unit: &'a str,
    time_unit: &'a str,
}

fn add_element(builder: Hdf5Builder, spec: ElementWriteSpec<'_>) -> Result<Hdf5Builder, H5mdError> {
    let value = DatasetBuilder::typed_data(
        element_path(spec.group, spec.element, VALUE),
        spec.shape,
        &spec.values,
    )?
    .attribute(AttributeBuilder::fixed_string(
        UNIT_ATTRIBUTE,
        spec.value_unit,
    ));
    let step = DatasetBuilder::typed_data(
        element_path(spec.group, spec.element, STEP),
        vec![to_u64(spec.steps.len())?],
        spec.steps,
    )?;
    let time = DatasetBuilder::typed_data(
        element_path(spec.group, spec.element, TIME),
        vec![to_u64(spec.times.len())?],
        spec.times,
    )?
    .attribute(AttributeBuilder::fixed_string(
        UNIT_ATTRIBUTE,
        spec.time_unit,
    ));
    Ok(builder.dataset(value).dataset(step).dataset(time))
}

fn add_box(
    builder: Hdf5Builder,
    frames: &[Timestep],
    options: &H5mdOptions,
    steps: &[i64],
    times: &[f64],
) -> Result<Hdf5Builder, H5mdError> {
    let scale = options.units.length_to_angstrom;
    let capacity = frames.len().checked_mul(9).ok_or(H5mdError::InvalidValue)?;
    let mut values = Vec::with_capacity(capacity);
    for frame in frames {
        let vectors = frame
            .cell
            .and_then(vectors_from_cell)
            .ok_or(H5mdError::InvalidValue)?;
        values.extend(vectors.into_iter().flatten().map(|value| value / scale));
    }
    let path = box_path(&options.particle_group, VALUE);
    let value =
        DatasetBuilder::typed_data(path, vec![to_u64(frames.len())?, 3, 3], &values)?.attribute(
            AttributeBuilder::fixed_string(UNIT_ATTRIBUTE, &options.units.length_unit),
        );
    let step = DatasetBuilder::typed_data(
        box_path(&options.particle_group, STEP),
        vec![to_u64(steps.len())?],
        steps,
    )?;
    let time = DatasetBuilder::typed_data(
        box_path(&options.particle_group, TIME),
        vec![to_u64(times.len())?],
        times,
    )?
    .attribute(AttributeBuilder::fixed_string(
        UNIT_ATTRIBUTE,
        &options.units.time_unit,
    ));
    let boundary = [PERIODIC, PERIODIC, PERIODIC];
    Ok(builder
        .group_attribute(
            box_group_path(&options.particle_group),
            AttributeBuilder::scalar(DIMENSION_ATTRIBUTE, 3_i32)?,
        )
        .group_attribute(
            box_group_path(&options.particle_group),
            AttributeBuilder::fixed_string_vector(BOUNDARY_ATTRIBUTE, &boundary)?,
        )
        .dataset(value)
        .dataset(step)
        .dataset(time))
}

fn collect_steps(frames: &[Timestep]) -> Result<Vec<i64>, H5mdError> {
    let steps = frames
        .iter()
        .map(|frame| i64::try_from(frame.frame).map_err(|_| H5mdError::InvalidValue))
        .collect::<Result<Vec<_>, _>>()?;
    if steps.windows(2).all(|pair| pair[0] < pair[1]) {
        Ok(steps)
    } else {
        Err(H5mdError::InvalidValue)
    }
}

fn collect_times(frames: &[Timestep], scale: f64) -> Result<Vec<f64>, H5mdError> {
    frames
        .iter()
        .map(|frame| {
            frame
                .time
                .map(|value| value / scale)
                .filter(|value| value.is_finite())
                .ok_or(H5mdError::InvalidValue)
        })
        .collect()
}

fn flatten_positions(frames: &[Timestep], scale: f64) -> Result<Vec<f64>, H5mdError> {
    convert_values(
        frames
            .iter()
            .flat_map(|frame| frame.positions.iter())
            .flat_map(|value| value.map(f64::from)),
        scale,
    )
}

fn flatten_optional(frames: &[Timestep], forces: bool, scale: f64) -> Result<Vec<f64>, H5mdError> {
    convert_values(
        frames
            .iter()
            .filter_map(|frame| {
                if forces {
                    frame.forces.as_ref()
                } else {
                    frame.velocities.as_ref()
                }
            })
            .flat_map(|values| values.iter())
            .flat_map(|value| value.map(f64::from)),
        scale,
    )
}

fn convert_values(values: impl Iterator<Item = f64>, scale: f64) -> Result<Vec<f64>, H5mdError> {
    values
        .map(|value| {
            let converted = value / scale;
            converted
                .is_finite()
                .then_some(converted)
                .ok_or(H5mdError::InvalidValue)
        })
        .collect()
}

fn validate_frame(frame: &Timestep) -> Result<(), H5mdError> {
    let vectors = frame
        .positions
        .iter()
        .chain(frame.velocities.iter().flatten())
        .chain(frame.forces.iter().flatten());
    let finite = vectors.flatten().all(|value| value.is_finite());
    let valid_cell = frame
        .cell
        .is_none_or(|cell| vectors_from_cell(cell).is_some());
    if finite && frame.time.is_some_and(f64::is_finite) && valid_cell {
        Ok(())
    } else {
        Err(H5mdError::InvalidValue)
    }
}

fn matches_vectors(values: Option<&[[f32; 3]]>, present: bool, atoms: usize) -> bool {
    values.is_some_and(|values| values.len() == atoms) == present
}

fn to_u64(value: usize) -> Result<u64, H5mdError> {
    u64::try_from(value).map_err(|_| H5mdError::InvalidValue)
}
