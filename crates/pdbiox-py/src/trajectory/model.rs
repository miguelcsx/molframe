//! Normalized Python trajectory object with zero-copy immutable array views.

use super::arrays;
use numpy::{
    PyArray1, PyArray2, PyArray3, PyReadonlyArray1, PyReadonlyArray2, PyReadonlyArray3,
    PyUntypedArrayMethods,
};
use pdbiox::UnitCell;
use pdbiox::traj::{
    FormatMetadata, Timestep, TrajectoryData, TrajectoryFormat, TrajectoryMetadata,
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::sync::Arc;

const INCONSISTENT_ATOMS: &str = "trajectory frames must have one atom count";
const INCONSISTENT_STREAM: &str =
    "optional trajectory streams must be present in every frame or none";
const VECTOR_WIDTH: usize = 3;
const CELL_WIDTH: usize = 6;
const INVALID_ARRAY_SHAPE: &str = "trajectory arrays must share the frame and atom dimensions";

#[pyclass(name = "Trajectory", frozen, skip_from_py_object)]
pub(crate) struct PyTrajectory {
    frames: usize,
    atoms: usize,
    positions: Arc<[f32]>,
    velocities: Option<Arc<[f32]>>,
    forces: Option<Arc<[f32]>>,
    time: Option<Arc<[f64]>>,
    dt: Option<Arc<[f64]>>,
    cells: Option<Arc<[f64]>>,
    steps: Option<Arc<[i64]>>,
    precision: Option<Arc<[f64]>>,
    frame_numbers: Box<[usize]>,
    data: Vec<std::collections::BTreeMap<Box<str>, pdbiox::traj::FrameValue>>,
    pub(super) format: Option<TrajectoryFormat>,
    format_metadata: FormatMetadata,
}

impl PyTrajectory {
    pub(crate) fn from_frames(
        frames: &[Timestep],
        format: Option<TrajectoryFormat>,
        metadata: TrajectoryMetadata,
    ) -> PyResult<Self> {
        let frame_count = frames.len();
        let atom_count = frames.first().map_or(0, |frame| frame.positions.len());
        if frames
            .iter()
            .any(|frame| frame.positions.len() != atom_count)
        {
            return Err(PyValueError::new_err(INCONSISTENT_ATOMS));
        }
        let positions = match flatten_vectors(frames.iter().map(|frame| Some(&frame.positions)))? {
            Some(positions) => positions,
            None => Arc::from([]),
        };
        let velocities = flatten_vectors(frames.iter().map(|frame| frame.velocities.as_ref()))?;
        let forces = flatten_vectors(frames.iter().map(|frame| frame.forces.as_ref()))?;
        let time = collect_scalars(frames.iter().map(|frame| frame.time))?;
        let dt = collect_scalars(frames.iter().map(|frame| frame.dt))?;
        let cells = collect_cells(frames.iter().map(|frame| frame.cell))?;
        validate_metadata(frame_count, metadata.steps.as_deref())?;
        let precision = coordinate_precision(&metadata.format, frame_count);
        validate_metadata(frame_count, precision.as_deref())?;
        let frame_numbers = frames.iter().map(|frame| frame.frame).collect();
        let data = frames.iter().map(|frame| frame.data.clone()).collect();
        Ok(Self {
            frames: frame_count,
            atoms: atom_count,
            positions,
            velocities,
            forces,
            time,
            dt,
            cells,
            steps: metadata.steps.map(Arc::from),
            precision: precision.map(Arc::from),
            frame_numbers,
            data,
            format,
            format_metadata: metadata.format,
        })
    }

    pub(crate) fn from_data(data: TrajectoryData) -> PyResult<Self> {
        Self::from_frames(&data.frames, Some(data.format), data.metadata)
    }

    pub(crate) fn to_data(&self, format: TrajectoryFormat) -> TrajectoryData {
        let mut frames = Vec::with_capacity(self.frames);
        for frame in 0..self.frames {
            frames.push(Timestep {
                frame: self.frame_numbers[frame],
                time: self.time.as_ref().map(|values| values[frame]),
                dt: self.dt.as_ref().map(|values| values[frame]),
                positions: vector_values(&self.positions, frame, self.atoms),
                velocities: self
                    .velocities
                    .as_ref()
                    .map(|values| vector_values(values, frame, self.atoms)),
                forces: self
                    .forces
                    .as_ref()
                    .map(|values| vector_values(values, frame, self.atoms)),
                cell: self.cells.as_ref().map(|values| cell_value(values, frame)),
                data: self.data[frame].clone(),
            });
        }
        TrajectoryData {
            format,
            frames,
            metadata: TrajectoryMetadata {
                steps: self.steps.as_deref().map(<[i64]>::to_vec),
                format: self.format_metadata.clone(),
            },
        }
    }
}

#[pymethods]
impl PyTrajectory {
    fn __len__(&self) -> usize {
        self.frames
    }

    fn __repr__(&self) -> String {
        format!(
            "Trajectory(format={:?}, frames={}, atoms={})",
            self.format, self.frames, self.atoms
        )
    }

    #[getter]
    const fn atom_count(&self) -> usize {
        self.atoms
    }

    #[getter]
    fn format(&self) -> PyResult<Option<super::types::PyTrajectoryFormat>> {
        self.format.map(TryInto::try_into).transpose()
    }

    #[getter]
    fn xyz<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f32>>> {
        arrays::vectors(py, Arc::clone(&self.positions), self.frames, self.atoms)
    }

    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray3<f32>>>> {
        self.velocities
            .as_ref()
            .map(|values| arrays::vectors(py, Arc::clone(values), self.frames, self.atoms))
            .transpose()
    }

    #[getter]
    fn forces<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray3<f32>>>> {
        self.forces
            .as_ref()
            .map(|values| arrays::vectors(py, Arc::clone(values), self.frames, self.atoms))
            .transpose()
    }

    #[getter]
    fn time<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        optional_scalars(py, self.time.as_ref())
    }

    #[getter]
    fn dt<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        optional_scalars(py, self.dt.as_ref())
    }

    #[getter]
    fn cells<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray2<f64>>>> {
        self.cells
            .as_ref()
            .map(|values| arrays::cells(py, Arc::clone(values), self.frames))
            .transpose()
    }

    #[getter]
    fn steps<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<i64>>>> {
        self.steps
            .as_ref()
            .map(|values| arrays::integers(py, Arc::clone(values)))
            .transpose()
    }

    #[getter]
    fn precision<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        optional_scalars(py, self.precision.as_ref())
    }

    #[new]
    #[pyo3(signature = (xyz, *, time=None, velocities=None, forces=None, cells=None, steps=None))]
    fn new(
        xyz: PyReadonlyArray3<'_, f32>,
        time: Option<PyReadonlyArray1<'_, f64>>,
        velocities: Option<PyReadonlyArray3<'_, f32>>,
        forces: Option<PyReadonlyArray3<'_, f32>>,
        cells: Option<PyReadonlyArray2<'_, f64>>,
        steps: Option<PyReadonlyArray1<'_, i64>>,
    ) -> PyResult<Self> {
        let shape = xyz.shape();
        if shape.len() != 3 || shape[2] != VECTOR_WIDTH {
            return Err(PyValueError::new_err(INVALID_ARRAY_SHAPE));
        }
        let frame_count = shape[0];
        let atom_count = shape[1];
        validate_vector_shape(velocities.as_ref(), frame_count, atom_count)?;
        validate_vector_shape(forces.as_ref(), frame_count, atom_count)?;
        validate_scalar_shape(time.as_ref(), frame_count)?;
        validate_scalar_shape(steps.as_ref(), frame_count)?;
        validate_cell_shape(cells.as_ref(), frame_count)?;
        let positions = xyz.as_array().to_owned();
        drop(xyz);
        let velocity_values = velocities.map(|values| values.as_array().to_owned());
        let force_values = forces.map(|values| values.as_array().to_owned());
        let time_values = time.map(|values| values.as_array().to_owned());
        let cell_values = cells.map(|values| values.as_array().to_owned());
        let mut frames = Vec::with_capacity(frame_count);
        for frame in 0..frame_count {
            frames.push(Timestep {
                frame,
                time: time_values.as_ref().map(|values| values[frame]),
                positions: (0..atom_count)
                    .map(|atom| {
                        positions
                            .slice(numpy::ndarray::s![frame, atom, ..])
                            .to_vec()
                    })
                    .map(|values| [values[0], values[1], values[2]])
                    .collect(),
                velocities: vector_frame(velocity_values.as_ref(), frame, atom_count),
                forces: vector_frame(force_values.as_ref(), frame, atom_count),
                cell: cell_values.as_ref().map(|values| UnitCell {
                    lengths: [values[[frame, 0]], values[[frame, 1]], values[[frame, 2]]],
                    angles: [values[[frame, 3]], values[[frame, 4]], values[[frame, 5]]],
                }),
                ..Timestep::default()
            });
        }
        let steps = steps.map(|values| values.as_array().iter().copied().collect());
        Self::from_frames(
            &frames,
            None,
            TrajectoryMetadata {
                steps,
                format: FormatMetadata::None,
            },
        )
    }
}

fn optional_scalars<'py>(
    py: Python<'py>,
    values: Option<&Arc<[f64]>>,
) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
    values
        .map(|values| arrays::scalars(py, Arc::clone(values)))
        .transpose()
}

fn flatten_vectors<'a>(
    values: impl Iterator<Item = Option<&'a Vec<[f32; VECTOR_WIDTH]>>>,
) -> PyResult<Option<Arc<[f32]>>> {
    let rows: Vec<_> = values.collect();
    if rows.iter().all(Option::is_none) {
        return Ok(None);
    }
    if rows.iter().any(Option::is_none) {
        return Err(PyValueError::new_err(INCONSISTENT_STREAM));
    }
    Ok(Some(Arc::from(
        rows.into_iter()
            .flatten()
            .flatten()
            .flat_map(|vector| vector.iter().copied())
            .collect::<Vec<_>>(),
    )))
}

fn collect_scalars(values: impl Iterator<Item = Option<f64>>) -> PyResult<Option<Arc<[f64]>>> {
    let values: Vec<_> = values.collect();
    if values.iter().all(Option::is_none) {
        return Ok(None);
    }
    values
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .map(|values| Some(Arc::from(values)))
        .ok_or_else(|| PyValueError::new_err(INCONSISTENT_STREAM))
}

fn collect_cells(values: impl Iterator<Item = Option<UnitCell>>) -> PyResult<Option<Arc<[f64]>>> {
    let values: Vec<_> = values.collect();
    if values.iter().all(Option::is_none) {
        return Ok(None);
    }
    let cells = values
        .into_iter()
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| PyValueError::new_err(INCONSISTENT_STREAM))?;
    let mut result = Vec::with_capacity(cells.len() * CELL_WIDTH);
    for cell in cells {
        result.extend(cell.lengths);
        result.extend(cell.angles);
    }
    Ok(Some(Arc::from(result)))
}

fn validate_metadata<T>(frames: usize, values: Option<&[T]>) -> PyResult<()> {
    if values.is_some_and(|values| values.len() != frames) {
        Err(PyValueError::new_err(INCONSISTENT_STREAM))
    } else {
        Ok(())
    }
}

fn validate_vector_shape(
    values: Option<&PyReadonlyArray3<'_, f32>>,
    frames: usize,
    atoms: usize,
) -> PyResult<()> {
    if values.is_some_and(|values| values.shape() != [frames, atoms, VECTOR_WIDTH]) {
        Err(PyValueError::new_err(INVALID_ARRAY_SHAPE))
    } else {
        Ok(())
    }
}

fn validate_scalar_shape<T: numpy::Element>(
    values: Option<&PyReadonlyArray1<'_, T>>,
    frames: usize,
) -> PyResult<()> {
    if values.is_some_and(|values| values.shape() != [frames]) {
        Err(PyValueError::new_err(INVALID_ARRAY_SHAPE))
    } else {
        Ok(())
    }
}

fn validate_cell_shape(values: Option<&PyReadonlyArray2<'_, f64>>, frames: usize) -> PyResult<()> {
    if values.is_some_and(|values| values.shape() != [frames, CELL_WIDTH]) {
        Err(PyValueError::new_err(INVALID_ARRAY_SHAPE))
    } else {
        Ok(())
    }
}

fn vector_frame(
    values: Option<&numpy::ndarray::Array3<f32>>,
    frame: usize,
    atoms: usize,
) -> Option<Vec<[f32; VECTOR_WIDTH]>> {
    values.map(|values| {
        (0..atoms)
            .map(|atom| {
                [
                    values[[frame, atom, 0]],
                    values[[frame, atom, 1]],
                    values[[frame, atom, 2]],
                ]
            })
            .collect()
    })
}

fn vector_values(values: &[f32], frame: usize, atoms: usize) -> Vec<[f32; VECTOR_WIDTH]> {
    let start = frame * atoms * VECTOR_WIDTH;
    values[start..start + atoms * VECTOR_WIDTH]
        .chunks_exact(VECTOR_WIDTH)
        .map(|values| [values[0], values[1], values[2]])
        .collect()
}

fn cell_value(values: &[f64], frame: usize) -> UnitCell {
    let start = frame * CELL_WIDTH;
    UnitCell {
        lengths: [values[start], values[start + 1], values[start + 2]],
        angles: [values[start + 3], values[start + 4], values[start + 5]],
    }
}

fn coordinate_precision(metadata: &FormatMetadata, frames: usize) -> Option<Vec<f64>> {
    match metadata {
        FormatMetadata::Xtc { precision } => {
            Some(precision.iter().copied().map(f64::from).collect())
        }
        FormatMetadata::Tng {
            compression_precision,
            ..
        } => Some(vec![*compression_precision; frames]),
        _ => None,
    }
}
