//! Native bindings for TRZ and GROMOS block trajectories.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;

pyo3::create_exception!(_native, TrzError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, GromosError, pyo3::exceptions::PyValueError);

#[pyclass(name = "TrzTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTrzTrajectory {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) has_forces: bool,
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
}

#[pymethods]
impl PyTrzTrajectory {
    #[new]
    fn new(title: String, has_forces: bool, frames: Vec<PyTimestep>) -> Self {
        Self {
            title,
            has_forces,
            frames,
        }
    }
}

impl TryFrom<pdbiox::traj::TrzTrajectory> for PyTrzTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::TrzTrajectory) -> Result<Self, Self::Error> {
        let frames = value
            .frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            title: value.title.into(),
            has_forces: value.has_forces,
            frames,
        })
    }
}

impl TryFrom<PyTrzTrajectory> for pdbiox::traj::TrzTrajectory {
    type Error = PyErr;

    fn try_from(value: PyTrzTrajectory) -> Result<Self, Self::Error> {
        let frames = value
            .frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            title: value.title.into(),
            has_forces: value.has_forces,
            frames,
        })
    }
}

#[pyfunction]
pub(crate) fn parse_trz(bytes: Vec<u8>) -> PyResult<PyTrzTrajectory> {
    pdbiox::traj::parse_trz(&bytes)
        .map_err(|error| TrzError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
pub(crate) fn write_trz(trajectory: PyTrzTrajectory) -> PyResult<Vec<u8>> {
    pdbiox::traj::write_trz(&trajectory.try_into()?)
        .map_err(|error| TrzError::new_err(error.to_string()))
}

#[pyclass(name = "GromosBoundary", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyGromosBoundary {
    Vacuum,
    Rectangular,
    Triclinic,
    TruncatedOctahedron,
}

impl From<pdbiox::traj::GromosBoundary> for PyGromosBoundary {
    fn from(value: pdbiox::traj::GromosBoundary) -> Self {
        match value {
            pdbiox::traj::GromosBoundary::Vacuum => Self::Vacuum,
            pdbiox::traj::GromosBoundary::Rectangular => Self::Rectangular,
            pdbiox::traj::GromosBoundary::Triclinic => Self::Triclinic,
            pdbiox::traj::GromosBoundary::TruncatedOctahedron => Self::TruncatedOctahedron,
        }
    }
}

#[pyclass(name = "GromosTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGromosTrajectory {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
    #[pyo3(get)]
    pub(crate) boundaries: Vec<Option<PyGromosBoundary>>,
}

impl TryFrom<pdbiox::traj::GromosTrajectory> for PyGromosTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::GromosTrajectory) -> Result<Self, Self::Error> {
        let frames = value
            .frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            title: value.title.into(),
            frames,
            boundaries: value
                .boundaries
                .into_iter()
                .map(|value| value.map(Into::into))
                .collect(),
        })
    }
}

#[pyfunction]
pub(crate) fn parse_gromos11_trc(source: &str) -> PyResult<PyGromosTrajectory> {
    pdbiox::traj::parse_gromos11_trc(source)
        .map_err(|error| GromosError::new_err(error.to_string()))?
        .try_into()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("TrzError", module.py().get_type::<TrzError>())?;
    module.add("GromosError", module.py().get_type::<GromosError>())?;
    module.add_class::<PyTrzTrajectory>()?;
    module.add_class::<PyGromosBoundary>()?;
    module.add_class::<PyGromosTrajectory>()?;
    module.add_function(wrap_pyfunction!(parse_trz, module)?)?;
    module.add_function(wrap_pyfunction!(write_trz, module)?)?;
    module.add_function(wrap_pyfunction!(parse_gromos11_trc, module)?)?;
    Ok(())
}
