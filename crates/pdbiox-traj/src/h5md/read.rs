//! H5MD decoding.

use hdf5_reader::{Dataset, Datatype, Hdf5File};

use crate::Timestep;
use crate::cell::cell_from_vectors;
use crate::numeric::f32_from_f64;

use super::schema::{
    BOUNDARY_ATTRIBUTE, CREATOR_GROUP, CREATOR_NAME_ATTRIBUTE, CREATOR_VERSION_ATTRIBUTE,
    DIMENSION_ATTRIBUTE, FORCE, H5MD_GROUP, H5MD_VERSION, PERIODIC, POSITION, STEP, TIME,
    UNIT_ATTRIBUTE, VALUE, VELOCITY, VERSION_ATTRIBUTE, box_edges_path, box_group_path, box_path,
    element_path,
};
use super::{H5mdError, H5mdOptions};

/// Complete decoded H5MD trajectory and source metadata.
#[derive(Clone, Debug, PartialEq)]
pub struct H5mdTrajectory {
    /// Canonical trajectory frames.
    pub frames: Vec<Timestep>,
    /// Source group, units and creator identity.
    pub metadata: H5mdMetadata,
}

/// H5MD metadata required for faithful rewriting.
#[derive(Clone, Debug, PartialEq)]
pub struct H5mdMetadata {
    /// Selected particle group and explicit unit conversions.
    pub options: H5mdOptions,
    /// Source creator name.
    pub creator_name: String,
    /// Source creator version.
    pub creator_version: String,
}

/// Parses the sole trajectory-bearing particle group in canonical pdbiox units.
///
/// Discovery is unambiguous only when exactly one group has positions. Unit
/// attributes are mandatory and must be canonical. Use
/// [`parse_h5md_with_options`] to declare a different unit system explicitly.
///
/// # Errors
///
/// Returns an error for invalid HDF5/H5MD metadata, paths, shapes, units,
/// clocks, cells, or numeric values.
pub fn parse_h5md(bytes: &[u8]) -> Result<Vec<Timestep>, H5mdError> {
    let file = Hdf5File::from_bytes(bytes)?;
    let options = H5mdOptions {
        particle_group: discover_particle_group(&file)?,
        ..H5mdOptions::default()
    };
    parse_record_file(&file, options).map(|trajectory| trajectory.frames)
}

/// Parses one explicitly selected H5MD particle group and unit system.
///
/// # Errors
///
/// Returns an error for invalid options or any schema violation.
pub fn parse_h5md_with_options(
    bytes: &[u8],
    options: &H5mdOptions,
) -> Result<Vec<Timestep>, H5mdError> {
    let file = Hdf5File::from_bytes(bytes)?;
    parse_record_file(&file, options.clone()).map(|trajectory| trajectory.frames)
}

/// Parses an explicitly selected particle group while retaining creator metadata.
///
/// # Errors
///
/// Returns an error for invalid options or any HDF5/H5MD schema violation.
pub fn parse_h5md_record_with_options(
    bytes: &[u8],
    options: &H5mdOptions,
) -> Result<H5mdTrajectory, H5mdError> {
    let file = Hdf5File::from_bytes(bytes)?;
    parse_record_file(&file, options.clone())
}

fn parse_record_file(file: &Hdf5File, options: H5mdOptions) -> Result<H5mdTrajectory, H5mdError> {
    let (creator_name, creator_version) = validate_metadata(file)?;
    let frames = parse_file(file, &options)?;
    Ok(H5mdTrajectory {
        frames,
        metadata: H5mdMetadata {
            options,
            creator_name,
            creator_version,
        },
    })
}

fn parse_file(file: &Hdf5File, options: &H5mdOptions) -> Result<Vec<Timestep>, H5mdError> {
    options.validate()?;
    let position = Element::read(
        file,
        ElementReadSpec {
            group: &options.particle_group,
            name: POSITION,
            unit: &options.units.length_unit,
            scale: options.units.length_to_angstrom,
            time_unit: &options.units.time_unit,
            time_scale: options.units.time_to_picosecond,
            required: true,
        },
    )?
    .ok_or_else(|| {
        H5mdError::MissingDataset(element_path(&options.particle_group, POSITION, VALUE))
    })?;
    let [frames, atoms, spatial] = position.shape.as_slice() else {
        return Err(shape_error(
            &position.value_path,
            "[frame, atom, 3]",
            &position.shape,
        ));
    };
    if *frames == 0 || *atoms == 0 || *spatial != 3 {
        return Err(shape_error(
            &position.value_path,
            "non-empty [frame, atom, 3]",
            &position.shape,
        ));
    }
    let frame_count = to_usize(*frames)?;
    let atom_count = to_usize(*atoms)?;
    position.validate_clock(frame_count)?;
    let velocity = Element::read(
        file,
        ElementReadSpec {
            group: &options.particle_group,
            name: VELOCITY,
            unit: &options.units.velocity_unit,
            scale: options.units.velocity_to_angstrom_per_picosecond,
            time_unit: &options.units.time_unit,
            time_scale: options.units.time_to_picosecond,
            required: false,
        },
    )?;
    let force = Element::read(
        file,
        ElementReadSpec {
            group: &options.particle_group,
            name: FORCE,
            unit: &options.units.force_unit,
            scale: options.units.force_to_kilojoule_per_mole_angstrom,
            time_unit: &options.units.time_unit,
            time_scale: options.units.time_to_picosecond,
            required: false,
        },
    )?;
    validate_companion(velocity.as_ref(), &position)?;
    validate_companion(force.as_ref(), &position)?;
    let cells = read_cells(file, options, &position)?;
    assemble(
        &position,
        velocity.as_ref(),
        force.as_ref(),
        cells.as_deref(),
        frame_count,
        atom_count,
    )
}

fn discover_particle_group(file: &Hdf5File) -> Result<String, H5mdError> {
    let particles = file.group("particles")?;
    let candidates = particles
        .groups()?
        .into_iter()
        .filter(|group| {
            file.dataset(&element_path(group.name(), POSITION, VALUE))
                .is_ok()
        })
        .map(|group| group.name().to_string())
        .collect::<Vec<_>>();
    match candidates.as_slice() {
        [group] => Ok(group.clone()),
        _ => Err(H5mdError::InvalidMetadata),
    }
}

struct Element {
    value_path: String,
    shape: Vec<u64>,
    values: Vec<f64>,
    steps: Vec<i64>,
    times: Vec<f64>,
}

#[derive(Clone, Copy)]
struct ElementReadSpec<'a> {
    group: &'a str,
    name: &'a str,
    unit: &'a str,
    scale: f64,
    time_unit: &'a str,
    time_scale: f64,
    required: bool,
}

impl Element {
    fn read(file: &Hdf5File, spec: ElementReadSpec<'_>) -> Result<Option<Self>, H5mdError> {
        let value_path = element_path(spec.group, spec.name, VALUE);
        let value = match file.dataset(&value_path) {
            Ok(dataset) => dataset,
            Err(_) if !spec.required => return Ok(None),
            Err(_) => return Err(H5mdError::MissingDataset(value_path)),
        };
        validate_unit(&value, &value_path, spec.unit)?;
        let shape = value.shape().to_vec();
        let mut values = read_float(&value)?;
        scale_values(&mut values, spec.scale)?;
        let step_path = element_path(spec.group, spec.name, STEP);
        let time_path = element_path(spec.group, spec.name, TIME);
        let steps = read_steps(&required_dataset(file, &step_path)?)?;
        let time_dataset = required_dataset(file, &time_path)?;
        validate_unit(&time_dataset, &time_path, spec.time_unit)?;
        let mut times = read_float(&time_dataset)?;
        scale_values(&mut times, spec.time_scale)?;
        Ok(Some(Self {
            value_path,
            shape,
            values,
            steps,
            times,
        }))
    }

    fn validate_clock(&self, frame_count: usize) -> Result<(), H5mdError> {
        let lengths_match = self.steps.len() == frame_count && self.times.len() == frame_count;
        let ordered = self.steps.windows(2).all(|pair| pair[0] < pair[1]);
        let representable = self.steps.iter().all(|step| usize::try_from(*step).is_ok());
        if lengths_match && ordered && representable {
            Ok(())
        } else {
            Err(H5mdError::InconsistentFrames)
        }
    }
}

fn validate_metadata(file: &Hdf5File) -> Result<(String, String), H5mdError> {
    let h5md = file.group(H5MD_GROUP)?;
    let version = h5md.attribute(VERSION_ATTRIBUTE)?.read_1d::<i32>()?;
    let creator = file.group(CREATOR_GROUP)?;
    let creator_name = creator.attribute(CREATOR_NAME_ATTRIBUTE)?.read_string()?;
    let creator_version = creator
        .attribute(CREATOR_VERSION_ATTRIBUTE)?
        .read_string()?;
    if version == H5MD_VERSION && !creator_name.is_empty() && !creator_version.is_empty() {
        Ok((creator_name, creator_version))
    } else {
        Err(H5mdError::InvalidMetadata)
    }
}

fn validate_companion(companion: Option<&Element>, position: &Element) -> Result<(), H5mdError> {
    if let Some(companion) = companion
        && (companion.shape != position.shape
            || companion.steps != position.steps
            || companion.times != position.times)
    {
        return Err(H5mdError::InconsistentFrames);
    }
    Ok(())
}

fn read_cells(
    file: &Hdf5File,
    options: &H5mdOptions,
    position: &Element,
) -> Result<Option<Vec<[[f64; 3]; 3]>>, H5mdError> {
    validate_box_metadata(file, &options.particle_group)?;
    let fixed_path = box_edges_path(&options.particle_group);
    if let Ok(dataset) = file.dataset(&fixed_path) {
        validate_unit(&dataset, &fixed_path, &options.units.length_unit)?;
        if dataset.shape() != [3, 3] {
            return Err(shape_error(&fixed_path, "[3, 3]", dataset.shape()));
        }
        let mut values = read_float(&dataset)?;
        scale_values(&mut values, options.units.length_to_angstrom)?;
        let vectors = matrix(&values)?;
        return Ok(Some(vec![vectors; position.steps.len()]));
    }
    let path = box_path(&options.particle_group, VALUE);
    let Ok(dataset) = file.dataset(&path) else {
        return Ok(None);
    };
    validate_unit(&dataset, &path, &options.units.length_unit)?;
    let frames = u64::try_from(position.steps.len()).map_err(|_| H5mdError::InvalidValue)?;
    if dataset.shape() != [frames, 3, 3] {
        return Err(shape_error(&path, "[frame, 3, 3]", dataset.shape()));
    }
    let step_path = box_path(&options.particle_group, STEP);
    let time_path = box_path(&options.particle_group, TIME);
    let time_dataset = required_dataset(file, &time_path)?;
    validate_unit(&time_dataset, &time_path, &options.units.time_unit)?;
    let mut times = read_float(&time_dataset)?;
    scale_values(&mut times, options.units.time_to_picosecond)?;
    if read_steps(&required_dataset(file, &step_path)?)? != position.steps
        || times != position.times
    {
        return Err(H5mdError::InconsistentFrames);
    }
    let mut values = read_float(&dataset)?;
    scale_values(&mut values, options.units.length_to_angstrom)?;
    values
        .chunks_exact(9)
        .map(matrix)
        .collect::<Result<Vec<_>, H5mdError>>()
        .map(Some)
}

fn validate_box_metadata(file: &Hdf5File, group: &str) -> Result<(), H5mdError> {
    let path = box_group_path(group);
    let Ok(box_group) = file.group(&path) else {
        return Ok(());
    };
    let dimension = box_group
        .attribute(DIMENSION_ATTRIBUTE)?
        .read_scalar::<i32>()?;
    let boundary = box_group.attribute(BOUNDARY_ATTRIBUTE)?.read_strings()?;
    if dimension == 3 && boundary.len() == 3 && boundary.iter().all(|value| value == PERIODIC) {
        Ok(())
    } else {
        Err(H5mdError::InvalidMetadata)
    }
}

fn matrix(values: &[f64]) -> Result<[[f64; 3]; 3], H5mdError> {
    if values.len() != 9 {
        return Err(H5mdError::InvalidValue);
    }
    Ok([
        [values[0], values[1], values[2]],
        [values[3], values[4], values[5]],
        [values[6], values[7], values[8]],
    ])
}

fn assemble(
    position: &Element,
    velocity: Option<&Element>,
    force: Option<&Element>,
    cells: Option<&[[[f64; 3]; 3]]>,
    frame_count: usize,
    atom_count: usize,
) -> Result<Vec<Timestep>, H5mdError> {
    let mut frames = Vec::with_capacity(frame_count);
    for index in 0..frame_count {
        let start = index
            .checked_mul(atom_count)
            .and_then(|value| value.checked_mul(3))
            .ok_or(H5mdError::InvalidValue)?;
        frames.push(Timestep {
            frame: usize::try_from(position.steps[index]).map_err(|_| H5mdError::InvalidValue)?,
            time: position.times.get(index).copied(),
            positions: vectors(&position.values, start, atom_count)?,
            velocities: velocity
                .map(|element| vectors(&element.values, start, atom_count))
                .transpose()?,
            forces: force
                .map(|element| vectors(&element.values, start, atom_count))
                .transpose()?,
            cell: cell_at(cells, index)?,
            ..Timestep::default()
        });
    }
    for (frame, times) in frames.iter_mut().skip(1).zip(position.times.windows(2)) {
        frame.dt = Some(times[1] - times[0]);
    }
    Ok(frames)
}

fn cell_at(
    cells: Option<&[[[f64; 3]; 3]]>,
    index: usize,
) -> Result<Option<pdbiox_core::structure::UnitCell>, H5mdError> {
    match cells.and_then(|values| values.get(index)) {
        Some(vectors) => cell_from_vectors(*vectors)
            .map(Some)
            .ok_or(H5mdError::InvalidValue),
        None if cells.is_some() => Err(H5mdError::InvalidValue),
        None => Ok(None),
    }
}

fn required_dataset(file: &Hdf5File, path: &str) -> Result<Dataset, H5mdError> {
    file.dataset(path)
        .map_err(|_| H5mdError::MissingDataset(path.to_string()))
}

fn validate_unit(dataset: &Dataset, path: &str, expected: &str) -> Result<(), H5mdError> {
    let declared = dataset
        .attribute(UNIT_ATTRIBUTE)
        .and_then(|attribute| attribute.read_string());
    let units = match declared {
        Ok(units) => units,
        Err(_) => "<missing>".to_string(),
    };
    if units == expected {
        Ok(())
    } else {
        Err(H5mdError::InvalidUnits {
            dataset: path.to_string(),
            units,
        })
    }
}

fn read_float(dataset: &Dataset) -> Result<Vec<f64>, H5mdError> {
    let values = match dataset.dtype() {
        Datatype::FloatingPoint { size: 4, .. } => dataset
            .read_array::<f32>()?
            .iter()
            .map(|value| f64::from(*value))
            .collect(),
        Datatype::FloatingPoint { size: 8, .. } => {
            dataset.read_array::<f64>()?.iter().copied().collect()
        }
        _ => return Err(H5mdError::InvalidValue),
    };
    Ok(values)
}

fn read_steps(dataset: &Dataset) -> Result<Vec<i64>, H5mdError> {
    let values = match dataset.dtype() {
        Datatype::FixedPoint {
            size: 4,
            signed: true,
            ..
        } => dataset
            .read_array::<i32>()?
            .iter()
            .map(|value| i64::from(*value))
            .collect(),
        Datatype::FixedPoint {
            size: 8,
            signed: true,
            ..
        } => dataset.read_array::<i64>()?.iter().copied().collect(),
        _ => return Err(H5mdError::InvalidValue),
    };
    Ok(values)
}

fn scale_values(values: &mut [f64], scale: f64) -> Result<(), H5mdError> {
    for value in values {
        *value *= scale;
        if !value.is_finite() {
            return Err(H5mdError::InvalidValue);
        }
    }
    Ok(())
}

fn vectors(values: &[f64], start: usize, atoms: usize) -> Result<Vec<[f32; 3]>, H5mdError> {
    let end = start
        .checked_add(atoms.checked_mul(3).ok_or(H5mdError::InvalidValue)?)
        .ok_or(H5mdError::InvalidValue)?;
    values
        .get(start..end)
        .ok_or(H5mdError::InvalidValue)?
        .chunks_exact(3)
        .map(|chunk| Ok([to_f32(chunk[0])?, to_f32(chunk[1])?, to_f32(chunk[2])?]))
        .collect()
}

fn shape_error(path: &str, expected: &'static str, found: &[u64]) -> H5mdError {
    H5mdError::InvalidShape {
        dataset: path.to_string(),
        expected,
        found: found.to_vec(),
    }
}

fn to_f32(value: f64) -> Result<f32, H5mdError> {
    f32_from_f64(value).ok_or(H5mdError::InvalidValue)
}

fn to_usize(value: u64) -> Result<usize, H5mdError> {
    usize::try_from(value).map_err(|_| H5mdError::InvalidValue)
}
