//! Native bindings for AMBER `NetCDF` trajectories.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;

pyo3::create_exception!(_native, AmberNetcdfError, pyo3::exceptions::PyValueError);

#[pyclass(name = "AmberNetcdfPrecision", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyAmberNetcdfPrecision {
    Single,
    Double,
}

impl From<PyAmberNetcdfPrecision> for pdbiox::traj::AmberNetcdfPrecision {
    fn from(value: PyAmberNetcdfPrecision) -> Self {
        match value {
            PyAmberNetcdfPrecision::Single => Self::Single,
            PyAmberNetcdfPrecision::Double => Self::Double,
        }
    }
}

impl From<pdbiox::traj::AmberNetcdfPrecision> for PyAmberNetcdfPrecision {
    fn from(value: pdbiox::traj::AmberNetcdfPrecision) -> Self {
        match value {
            pdbiox::traj::AmberNetcdfPrecision::Single => Self::Single,
            pdbiox::traj::AmberNetcdfPrecision::Double => Self::Double,
        }
    }
}

#[pyclass(name = "AmberNetcdfMetadata", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberNetcdfMetadata {
    #[pyo3(get)]
    pub(crate) precision: PyAmberNetcdfPrecision,
    #[pyo3(get)]
    pub(crate) program: Option<String>,
    #[pyo3(get)]
    pub(crate) program_version: Option<String>,
}

impl From<pdbiox::traj::AmberNetcdfMetadata> for PyAmberNetcdfMetadata {
    fn from(value: pdbiox::traj::AmberNetcdfMetadata) -> Self {
        Self {
            precision: value.precision.into(),
            program: value.program.map(Into::into),
            program_version: value.program_version.map(Into::into),
        }
    }
}

#[pyclass(name = "AmberNetcdfTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberNetcdfTrajectory {
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
    #[pyo3(get)]
    pub(crate) metadata: PyAmberNetcdfMetadata,
}

#[pyclass(name = "AmberNetcdfWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAmberNetcdfWriteOptions {
    #[pyo3(get)]
    pub(crate) precision: PyAmberNetcdfPrecision,
    #[pyo3(get)]
    pub(crate) program: String,
    #[pyo3(get)]
    pub(crate) program_version: String,
}

#[pymethods]
impl PyAmberNetcdfWriteOptions {
    #[new]
    #[pyo3(signature = (*, precision=None, program=None, program_version=None))]
    fn new(
        precision: Option<PyAmberNetcdfPrecision>,
        program: Option<String>,
        program_version: Option<String>,
    ) -> Self {
        let defaults = pdbiox::traj::AmberNetcdfWriteOptions::default();
        Self {
            precision: match precision {
                Some(value) => value,
                None => defaults.precision.into(),
            },
            program: match program {
                Some(value) => value,
                None => defaults.program,
            },
            program_version: match program_version {
                Some(value) => value,
                None => defaults.program_version,
            },
        }
    }
}

impl From<PyAmberNetcdfWriteOptions> for pdbiox::traj::AmberNetcdfWriteOptions {
    fn from(value: PyAmberNetcdfWriteOptions) -> Self {
        Self {
            precision: value.precision.into(),
            program: value.program,
            program_version: value.program_version,
        }
    }
}

impl TryFrom<pdbiox::traj::AmberNetcdfTrajectory> for PyAmberNetcdfTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::AmberNetcdfTrajectory) -> Result<Self, Self::Error> {
        let frames = value
            .frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            frames,
            metadata: value.metadata.into(),
        })
    }
}

#[pyfunction]
pub(crate) fn parse_amber_netcdf(py: Python<'_>, bytes: Vec<u8>) -> PyResult<Vec<PyTimestep>> {
    py.detach(move || -> PyResult<Vec<PyTimestep>> {
        pdbiox::traj::parse_amber_netcdf(&bytes)
            .map_err(|error| AmberNetcdfError::new_err(error.to_string()))?
            .into_iter()
            .map(TryInto::try_into)
            .collect()
    })
}

#[pyfunction]
pub(crate) fn parse_amber_netcdf_record(
    py: Python<'_>,
    bytes: Vec<u8>,
) -> PyResult<PyAmberNetcdfTrajectory> {
    py.detach(move || -> PyResult<PyAmberNetcdfTrajectory> {
        pdbiox::traj::parse_amber_netcdf_record(&bytes)
            .map_err(|error| AmberNetcdfError::new_err(error.to_string()))?
            .try_into()
    })
}

#[pyfunction]
pub(crate) fn write_amber_netcdf(
    py: Python<'_>,
    frames: Vec<PyTimestep>,
    options: PyAmberNetcdfWriteOptions,
) -> PyResult<Vec<u8>> {
    py.detach(move || -> PyResult<Vec<u8>> {
        let frames = frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        pdbiox::traj::write_amber_netcdf(&frames, options.into())
            .map_err(|error| AmberNetcdfError::new_err(error.to_string()))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "AmberNetcdfError",
        module.py().get_type::<AmberNetcdfError>(),
    )?;
    module.add_class::<PyAmberNetcdfPrecision>()?;
    module.add_class::<PyAmberNetcdfMetadata>()?;
    module.add_class::<PyAmberNetcdfTrajectory>()?;
    module.add_class::<PyAmberNetcdfWriteOptions>()?;
    module.add_function(wrap_pyfunction!(parse_amber_netcdf, module)?)?;
    module.add_function(wrap_pyfunction!(parse_amber_netcdf_record, module)?)?;
    module.add_function(wrap_pyfunction!(write_amber_netcdf, module)?)?;
    Ok(())
}
