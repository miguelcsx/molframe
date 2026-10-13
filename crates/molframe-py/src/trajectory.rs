//! Mechanical adapters for reading whole trajectories into `NumPy` arrays.

use crate::analysis_result::PyAnalysis;
use crate::policy::PyAnalysisPolicy;
use molframe::trajectory::{DenseChunk, FrameAlignment, FrameView, analyse_rmsd_to_reference_view};
use numpy::{PyArray1, PyArray2, PyArray3, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
use pyo3::prelude::*;

/// Frames of one trajectory as read-only arrays, in ångström and picoseconds.
#[derive(Debug)]
#[pyclass(
    name = "Trajectory",
    frozen,
    skip_from_py_object,
    module = "molframe.trajectory"
)]
pub(crate) struct PyTrajectory {
    format: String,
    positions: Py<PyArray3<f32>>,
    times: Py<PyArray1<f64>>,
    cells: Option<Py<PyArray2<f64>>>,
    frame_indices: Py<PyArray1<u64>>,
}

impl PyTrajectory {
    /// Hands a filled chunk's buffers to `NumPy` without copying them.
    pub(crate) fn from_chunk(py: Python<'_>, format: String, chunk: DenseChunk) -> PyResult<Self> {
        let frames = chunk.frames();
        let positions =
            PyArray1::from_vec(py, chunk.coordinates).reshape([frames, chunk.atoms, 3])?;
        positions.readwrite().make_nonwriteable();
        let times = PyArray1::from_vec(py, chunk.times);
        times.readwrite().make_nonwriteable();
        let frame_indices = PyArray1::from_vec(py, chunk.indices);
        frame_indices.readwrite().make_nonwriteable();
        let cells = chunk
            .cells
            .map(|flat| -> PyResult<_> {
                let array = PyArray1::from_vec(py, flat).reshape([frames, 6])?;
                array.readwrite().make_nonwriteable();
                Ok(array.unbind())
            })
            .transpose()?;
        Ok(Self {
            format,
            positions: positions.unbind(),
            times: times.unbind(),
            cells,
            frame_indices: frame_indices.unbind(),
        })
    }
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

    /// Source frame ordinal of each row, as `uint64`; strided and sliced
    /// reads keep the original numbering.
    #[getter]
    fn frame_indices(&self, py: Python<'_>) -> Py<PyArray1<u64>> {
        self.frame_indices.clone_ref(py)
    }

    /// Unit cell of each frame as `(frames, 6)` rows of `a, b, c, alpha, beta,
    /// gamma` (ångström, degrees); `nan` rows where a frame has no cell, and
    /// `None` when no frame has one.
    #[getter]
    fn cells(&self, py: Python<'_>) -> Option<Py<PyArray2<f64>>> {
        self.cells.as_ref().map(|cells| cells.clone_ref(py))
    }

    /// The cell of one frame as a `UnitCell`, or `None` when it has none.
    fn cell(&self, py: Python<'_>, frame: isize) -> PyResult<Option<crate::crystal::PyUnitCell>> {
        let Some(cells) = &self.cells else {
            return Ok(None);
        };
        let frames = self.n_frames(py);
        let index = if frame < 0 {
            frame.checked_add_unsigned(frames)
        } else {
            Some(frame)
        };
        let row = index
            .and_then(|index| usize::try_from(index).ok())
            .filter(|index| *index < frames)
            .ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("frame index out of range"))?;
        let array = cells.bind(py).readonly();
        let values = array.as_array();
        let at = |column: usize| values[[row, column]];
        if (0..6).any(|column| at(column).is_nan()) {
            return Ok(None);
        }
        crate::crystal::PyUnitCell::from_cell(&molframe::UnitCell {
            lengths: [at(0), at(1), at(2)],
            angles: [at(3), at(4), at(5)],
        })
        .map(Some)
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

/// RMSD of every frame against one reference frame, with provenance.
///
/// `positions` is any array-like of shape `(frames, atoms, 3)`; strided views and
/// `float64` are converted to contiguous `float32`. With `align=True` each frame is
/// rigidly fitted onto the reference first, so the series measures shape change
/// rather than drift and rotation.
#[pyfunction]
#[pyo3(signature = (positions, *, reference=0, align=true, policy=None))]
fn rmsd(
    py: Python<'_>,
    positions: &Bound<'_, PyAny>,
    reference: usize,
    align: bool,
    policy: Option<PyRef<'_, PyAnalysisPolicy>>,
) -> PyResult<PyAnalysis> {
    let array = crate::trajectory_ensemble::positions_f32(positions)?;
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
    crate::trajectory_ensemble::register(module)?;
    module.add_class::<PyTrajectory>()?;
    module.add_function(wrap_pyfunction!(crate::trajectory_stream::read, module)?)?;
    crate::trajectory_stream::register(module)?;
    module.add_function(wrap_pyfunction!(rmsd, module)?)
}
