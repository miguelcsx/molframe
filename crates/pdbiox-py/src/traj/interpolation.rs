//! Python trajectory interpolation over topology-aligned coordinate frames.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pyclass(name = "TrajectoryInterpolation", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyTrajectoryInterpolation {
    Linear,
    CentripetalCatmullRom,
}

impl From<PyTrajectoryInterpolation> for pdbiox::traj::TrajectoryInterpolation {
    fn from(value: PyTrajectoryInterpolation) -> Self {
        match value {
            PyTrajectoryInterpolation::Linear => Self::Linear,
            PyTrajectoryInterpolation::CentripetalCatmullRom => Self::CentripetalCatmullRom,
        }
    }
}

#[pyclass(
    name = "TrajectoryInterpolationError",
    frozen,
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyTrajectoryInterpolationError {
    TopologyMismatch,
    NonFiniteInput,
    FractionOutOfRange,
}

#[pyfunction]
fn interpolate_trajectory_frames(
    py: Python<'_>,
    frames: Vec<Vec<[f32; 3]>>,
    fraction: f32,
    interpolation: PyTrajectoryInterpolation,
) -> PyResult<Vec<[f32; 3]>> {
    py.detach(move || -> PyResult<Vec<[f32; 3]>> {
        let [first, second, third, fourth] = frames.as_slice() else {
            return Err(PyValueError::new_err("exactly four frames are required"));
        };
        pdbiox::traj::interpolate_trajectory_frames(
            [first, second, third, fourth],
            fraction,
            interpolation.into(),
        )
        .map_err(|error| PyValueError::new_err(error.to_string()))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyTrajectoryInterpolation>()?;
    module.add_class::<PyTrajectoryInterpolationError>()?;
    module.add_function(wrap_pyfunction!(interpolate_trajectory_frames, module)?)
}
