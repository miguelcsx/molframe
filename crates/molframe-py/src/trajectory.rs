//! Mechanical adapters for reading whole trajectories into `NumPy` arrays.

use molframe::trajectory::{TrajectoryFormat, TrajectoryReadOptions, read_trajectory_materialized};
use numpy::{PyArray1, PyArray3, PyArrayMethods, PyUntypedArrayMethods};
use pyo3::{exceptions::PyValueError, prelude::*};
use std::path::PathBuf;

/// Frames of one trajectory as read-only arrays, in ångström and picoseconds.
#[derive(Debug)]
#[pyclass(
    name = "Trajectory",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
struct PyTrajectory {
    format: &'static str,
    positions: Py<PyArray3<f32>>,
    times: Py<PyArray1<f64>>,
}

#[pymethods]
impl PyTrajectory {
    /// Container the frames were read from.
    #[getter]
    fn format(&self) -> &'static str {
        self.format
    }

    /// Coordinates with shape `(frames, atoms, 3)`.
    #[getter]
    fn positions(&self, py: Python<'_>) -> Py<PyArray3<f32>> {
        self.positions.clone_ref(py)
    }

    /// Simulation time of each frame in picoseconds; `nan` where unrecorded.
    #[getter]
    fn times(&self, py: Python<'_>) -> Py<PyArray1<f64>> {
        self.times.clone_ref(py)
    }

    #[getter]
    fn n_frames(&self, py: Python<'_>) -> usize {
        self.positions.bind(py).shape()[0]
    }

    #[getter]
    fn n_atoms(&self, py: Python<'_>) -> usize {
        self.positions.bind(py).shape()[1]
    }

    fn __len__(&self, py: Python<'_>) -> usize {
        self.n_frames(py)
    }

    fn __repr__(&self, py: Python<'_>) -> String {
        format!(
            "Trajectory(format={}, frames={}, atoms={})",
            self.format,
            self.n_frames(py),
            self.n_atoms(py)
        )
    }
}

fn format_name(format: TrajectoryFormat) -> &'static str {
    match format {
        TrajectoryFormat::Xtc => "xtc",
        TrajectoryFormat::Trr => "trr",
        TrajectoryFormat::Dcd => "dcd",
        TrajectoryFormat::Tng => "tng",
        TrajectoryFormat::Gro => "gro",
        TrajectoryFormat::Xyz => "xyz",
        TrajectoryFormat::LammpsDump => "lammps_dump",
        TrajectoryFormat::AmberNetcdf => "netcdf",
        _ => "other",
    }
}

fn parse_format(name: &str) -> PyResult<TrajectoryFormat> {
    match name {
        "xtc" => Ok(TrajectoryFormat::Xtc),
        "trr" => Ok(TrajectoryFormat::Trr),
        "dcd" => Ok(TrajectoryFormat::Dcd),
        "tng" => Ok(TrajectoryFormat::Tng),
        "gro" => Ok(TrajectoryFormat::Gro),
        "xyz" => Ok(TrajectoryFormat::Xyz),
        "lammps_dump" => Ok(TrajectoryFormat::LammpsDump),
        "netcdf" => Ok(TrajectoryFormat::AmberNetcdf),
        _ => Err(PyValueError::new_err(
            "format must be xtc, trr, dcd, tng, gro, xyz, lammps_dump or netcdf",
        )),
    }
}

/// Reads every frame of a self-describing trajectory file.
#[pyfunction]
#[pyo3(signature = (path, *, format=None))]
fn read(py: Python<'_>, path: PathBuf, format: Option<&str>) -> PyResult<PyTrajectory> {
    let options = TrajectoryReadOptions {
        format: format.map(parse_format).transpose()?,
        ..TrajectoryReadOptions::default()
    };
    let data = py
        .detach(move || read_trajectory_materialized(&path, &options))
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let atoms = data.frames.first().map_or(0, |frame| frame.positions.len());
    if let Some(frame) = data
        .frames
        .iter()
        .find(|frame| frame.positions.len() != atoms)
    {
        return Err(PyValueError::new_err(format!(
            "frame {} has {} atoms, expected {atoms}",
            frame.frame,
            frame.positions.len()
        )));
    }
    let frames = data.frames.len();
    let mut flat = Vec::with_capacity(frames * atoms * 3);
    let mut times = Vec::with_capacity(frames);
    for frame in &data.frames {
        flat.extend(frame.positions.iter().flatten());
        times.push(frame.time.unwrap_or(f64::NAN));
    }
    let positions = PyArray1::from_vec(py, flat).reshape([frames, atoms, 3])?;
    positions.readwrite().make_nonwriteable();
    let times = PyArray1::from_vec(py, times);
    times.readwrite().make_nonwriteable();
    Ok(PyTrajectory {
        format: format_name(data.format),
        positions: positions.unbind(),
        times: times.unbind(),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyTrajectory>()?;
    module.add_function(wrap_pyfunction!(read, module)?)
}
