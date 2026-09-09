//! Normalized Python trajectory object with zero-copy immutable array views.

use self::storage::{CellStorage, IntegerStorage, ScalarStorage, VectorStorage};
use super::generic::PyFrame;
use super::reader_types::PyTimestep;
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

mod storage;

const INCONSISTENT_ATOMS: &str = "trajectory frames must have one atom count";
const INCONSISTENT_STREAM: &str =
    "optional trajectory streams must be present in every frame or none";
const VECTOR_WIDTH: usize = 3;
const CELL_WIDTH: usize = 6;
const INVALID_ARRAY_SHAPE: &str = "trajectory arrays must share the frame and atom dimensions";

#[pyclass(name = "Trajectory", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyTrajectory {
    frames: usize,
    atoms: usize,
    positions: VectorStorage,
    velocities: Option<VectorStorage>,
    forces: Option<VectorStorage>,
    time: Option<ScalarStorage>,
    dt: Option<ScalarStorage>,
    cells: Option<CellStorage>,
    steps: Option<IntegerStorage>,
    precision: Option<ScalarStorage>,
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
            positions: VectorStorage::owned(positions),
            velocities: velocities.map(VectorStorage::owned),
            forces: forces.map(VectorStorage::owned),
            time: time.map(ScalarStorage::owned),
            dt: dt.map(ScalarStorage::owned),
            cells: cells.map(CellStorage::owned),
            steps: metadata.steps.map(Arc::from).map(IntegerStorage::owned),
            precision: precision.map(Arc::from).map(ScalarStorage::owned),
            frame_numbers,
            data,
            format,
            format_metadata: metadata.format,
        })
    }

    pub(crate) fn from_data(data: TrajectoryData) -> PyResult<Self> {
        Self::from_frames(&data.frames, Some(data.format), data.metadata)
    }

    pub(crate) fn to_data(
        &self,
        py: Python<'_>,
        format: TrajectoryFormat,
    ) -> PyResult<TrajectoryData> {
        let frames = (0..self.frames)
            .map(|frame| self.timestep_value(py, frame))
            .collect::<PyResult<_>>()?;
        Ok(TrajectoryData {
            format,
            frames,
            metadata: TrajectoryMetadata {
                steps: self
                    .steps
                    .as_ref()
                    .map(|values| values.values(py, self.frames))
                    .transpose()?,
                format: self.format_metadata.clone(),
            },
        })
    }

    pub(crate) fn timesteps(&self, py: Python<'_>) -> PyResult<Vec<PyTimestep>> {
        (0..self.frames)
            .map(|frame| self.timestep_value(py, frame)?.try_into())
            .collect()
    }

    pub(crate) fn source_format_metadata(&self) -> FormatMetadata {
        self.format_metadata.clone()
    }

    pub(crate) fn source_steps(&self, py: Python<'_>) -> PyResult<Option<Vec<i64>>> {
        self.steps
            .as_ref()
            .map(|values| values.values(py, self.frames))
            .transpose()
    }

    pub(crate) fn native_trajectory(&self, py: Python<'_>) -> PyResult<pdbiox::traj::Trajectory> {
        let frames = (0..self.frames)
            .map(|frame| {
                Ok(pdbiox::traj::Frame {
                    positions: self.positions.frame(py, frame, self.frames, self.atoms)?,
                })
            })
            .collect::<PyResult<_>>()?;
        pdbiox::traj::Trajectory::from_frames(frames)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    fn timestep_value(&self, py: Python<'_>, frame: usize) -> PyResult<Timestep> {
        Ok(Timestep {
            frame: self.frame_numbers[frame],
            time: self
                .time
                .as_ref()
                .map(|values| values.value(py, frame, self.frames))
                .transpose()?,
            dt: self
                .dt
                .as_ref()
                .map(|values| values.value(py, frame, self.frames))
                .transpose()?,
            positions: self.positions.frame(py, frame, self.frames, self.atoms)?,
            velocities: self
                .velocities
                .as_ref()
                .map(|values| values.frame(py, frame, self.frames, self.atoms))
                .transpose()?,
            forces: self
                .forces
                .as_ref()
                .map(|values| values.frame(py, frame, self.frames, self.atoms))
                .transpose()?,
            cell: self
                .cells
                .as_ref()
                .map(|values| values.cell(py, frame, self.frames))
                .transpose()?,
            data: self.data[frame].clone(),
        })
    }

    fn timestep_at(&self, py: Python<'_>, frame: usize) -> PyResult<Option<Timestep>> {
        (frame < self.frames)
            .then(|| self.timestep_value(py, frame))
            .transpose()
    }

    pub(crate) fn python_format(&self) -> PyResult<Option<super::types::PyTrajectoryFormat>> {
        self.format.map(TryInto::try_into).transpose()
    }

    pub(crate) fn select_indices(&self, py: Python<'_>, indices: Vec<usize>) -> PyResult<Self> {
        let format = match self.format {
            Some(format) => format,
            None => TrajectoryFormat::Xtc,
        };
        let data = self
            .to_data(py, format)?
            .select_frames(&indices)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Self::from_data(data)
    }
}

#[pymethods]
impl PyTrajectory {
    fn __len__(&self) -> usize {
        self.frames
    }

    fn __repr__(&self) -> String {
        format!(
            "Trajectory(format={:?}, frames={}, atoms={}, view={})",
            self.format,
            self.frames,
            self.atoms,
            self.positions.is_retained()
        )
    }

    #[getter]
    const fn atom_count(&self) -> usize {
        self.atoms
    }

    /// Whether coordinate arrays are retained from the caller rather than copied.
    #[getter]
    fn is_view(&self) -> bool {
        self.positions.is_retained()
    }

    #[getter]
    fn format(&self) -> PyResult<Option<super::types::PyTrajectoryFormat>> {
        self.python_format()
    }

    #[getter]
    fn xyz<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray3<f32>>> {
        self.positions.array(py, self.frames, self.atoms)
    }

    #[getter]
    fn velocities<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray3<f32>>>> {
        self.velocities
            .as_ref()
            .map(|values| values.array(py, self.frames, self.atoms))
            .transpose()
    }

    #[getter]
    fn forces<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray3<f32>>>> {
        self.forces
            .as_ref()
            .map(|values| values.array(py, self.frames, self.atoms))
            .transpose()
    }

    #[getter]
    fn time<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        self.time
            .as_ref()
            .map(|values| values.array(py, self.frames))
            .transpose()
    }

    #[getter]
    fn dt<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        self.dt
            .as_ref()
            .map(|values| values.array(py, self.frames))
            .transpose()
    }

    #[getter]
    fn cells<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray2<f64>>>> {
        self.cells
            .as_ref()
            .map(|values| values.array(py, self.frames))
            .transpose()
    }

    #[getter]
    fn steps<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<i64>>>> {
        self.steps
            .as_ref()
            .map(|values| values.array(py))
            .transpose()
    }

    #[getter]
    fn precision<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyArray1<f64>>>> {
        self.precision
            .as_ref()
            .map(|values| values.array(py, self.frames))
            .transpose()
    }

    #[new]
    #[pyo3(signature = (xyz, *, time=None, velocities=None, forces=None, cells=None, steps=None, copy=false))]
    fn new(
        xyz: PyReadonlyArray3<'_, f32>,
        time: Option<PyReadonlyArray1<'_, f64>>,
        velocities: Option<PyReadonlyArray3<'_, f32>>,
        forces: Option<PyReadonlyArray3<'_, f32>>,
        cells: Option<PyReadonlyArray2<'_, f64>>,
        steps: Option<PyReadonlyArray1<'_, i64>>,
        copy: bool,
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
        let positions = VectorStorage::from_array(xyz, copy)?;
        let velocities = velocities
            .map(|values| VectorStorage::from_array(values, copy))
            .transpose()?;
        let forces = forces
            .map(|values| VectorStorage::from_array(values, copy))
            .transpose()?;
        let time = time
            .map(|values| ScalarStorage::from_array(values, copy))
            .transpose()?;
        let cells = cells
            .map(|values| CellStorage::from_array(values, copy))
            .transpose()?;
        let steps = steps
            .map(|values| IntegerStorage::from_array(values, copy))
            .transpose()?;
        Ok(Self {
            frames: frame_count,
            atoms: atom_count,
            positions,
            velocities,
            forces,
            time,
            dt: None,
            cells,
            steps,
            precision: None,
            frame_numbers: (0..frame_count).collect(),
            data: (0..frame_count)
                .map(|_| std::collections::BTreeMap::new())
                .collect(),
            format: None,
            format_metadata: FormatMetadata::None,
        })
    }

    #[staticmethod]
    #[pyo3(name = "from_frames")]
    fn from_frame_objects(frames: Vec<PyFrame>) -> PyResult<Self> {
        let frames = frames
            .into_iter()
            .enumerate()
            .map(|(frame, value)| Timestep {
                frame,
                positions: value.positions,
                ..Timestep::default()
            })
            .collect::<Vec<_>>();
        Self::from_frames(
            &frames,
            None,
            TrajectoryMetadata {
                steps: None,
                format: FormatMetadata::None,
            },
        )
    }

    fn frame(&self, py: Python<'_>, index: usize) -> PyResult<Option<PyFrame>> {
        let frame = self.timestep_at(py, index)?;
        Ok(frame.map(|value| PyFrame {
            positions: value.positions,
        }))
    }

    fn frames(&self, py: Python<'_>) -> PyResult<Vec<PyFrame>> {
        Ok(self
            .timesteps(py)?
            .into_iter()
            .map(|value| PyFrame {
                positions: value.positions,
            })
            .collect())
    }

    fn slice(&self, py: Python<'_>, indices: Vec<usize>) -> PyResult<Self> {
        self.select_indices(py, indices)
    }
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
