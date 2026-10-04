//! Mechanical adapters for reading whole trajectories into `NumPy` arrays.

use crate::analysis_result::PyAnalysis;
use crate::policy::PyAnalysisPolicy;
use molframe::trajectory::{
    FrameAlignment, FrameView, TrajectoryReadOptions, analyse_rmsd_to_reference_view,
    read_trajectory_materialized,
};
use numpy::{PyArray1, PyArray3, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
use pyo3::prelude::*;
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
    format: String,
    positions: Py<PyArray3<f32>>,
    times: Py<PyArray1<f64>>,
}

#[pymethods]
impl PyTrajectory {
    /// Container the frames were read from.
    #[getter]
    fn format(&self) -> String {
        self.format.clone()
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

/// Reads every frame of a self-describing trajectory file.
#[pyfunction]
#[pyo3(signature = (path, *, format=None))]
fn read(py: Python<'_>, path: PathBuf, format: Option<&str>) -> PyResult<PyTrajectory> {
    let options = TrajectoryReadOptions {
        format: format
            .map(|name| name.parse().map_err(crate::error::kernel))
            .transpose()?,
        ..TrajectoryReadOptions::default()
    };
    let data = py
        .detach(move || read_trajectory_materialized(&path, &options))
        .map_err(crate::error::failure)?;
    let dense = data.dense_positions().map_err(crate::error::kernel)?;
    let positions =
        PyArray1::from_vec(py, dense.coordinates).reshape([dense.frames, dense.atoms, 3])?;
    positions.readwrite().make_nonwriteable();
    let times = PyArray1::from_vec(py, dense.times);
    times.readwrite().make_nonwriteable();
    Ok(PyTrajectory {
        format: data.format.name().replace('-', "_"),
        positions: positions.unbind(),
        times: times.unbind(),
    })
}

/// RMSD of every frame against one reference frame, with provenance.
///
/// `positions` has shape `(frames, atoms, 3)`. With `align=True` each frame is
/// rigidly fitted onto the reference first, so the series measures shape change
/// rather than drift and rotation.
#[pyfunction]
#[pyo3(signature = (positions, *, reference=0, align=true, policy=None))]
fn rmsd(
    py: Python<'_>,
    positions: &Bound<'_, PyArray3<f32>>,
    reference: usize,
    align: bool,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let array = positions.readonly();
    let shape = array.shape();
    let flat = array.as_slice().map_err(|_| {
        crate::error::value("positions must be C-contiguous; call numpy.ascontiguousarray")
    })?;
    if shape[2] != 3 {
        return Err(crate::error::value(
            "positions must have shape (frames, atoms, 3)",
        ));
    }
    let view = FrameView::new(flat.as_chunks::<3>().0, shape[0], shape[1])
        .map_err(crate::error::failure)?;
    let alignment = if align {
        FrameAlignment::Rigid
    } else {
        FrameAlignment::None
    };
    let policy = crate::policy::policy_of(policy);
    let analysis = py
        .detach(|| analyse_rmsd_to_reference_view(view, reference, alignment, &policy))
        .map_err(crate::error::failure)?;
    let series = match analysis.value() {
        Some(values) => {
            let series = values.to_pyarray(py);
            series.readwrite().make_nonwriteable();
            Some(series.into_any().unbind())
        }
        None => None,
    };
    Ok(PyAnalysis::new(&analysis, series))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyTrajectory>()?;
    module.add_function(wrap_pyfunction!(read, module)?)?;
    module.add_function(wrap_pyfunction!(rmsd, module)?)
}
