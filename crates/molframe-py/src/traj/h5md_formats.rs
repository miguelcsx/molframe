//! Native bindings for H5MD unit-aware trajectory records.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;

pyo3::create_exception!(_native, H5mdError, pyo3::exceptions::PyValueError);

#[pyclass(name = "H5mdUnitSystem", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyH5mdUnitSystem {
    pub(crate) inner: molframe::traj::H5mdUnitSystem,
}

#[pymethods]
impl PyH5mdUnitSystem {
    #[new]
    #[pyo3(signature = (length_unit, length_to_angstrom, time_unit, time_to_picosecond, velocity_unit, velocity_to_angstrom_per_picosecond, force_unit, force_to_kilojoule_per_mole_angstrom))]
    fn new(
        length_unit: String,
        length_to_angstrom: f64,
        time_unit: String,
        time_to_picosecond: f64,
        velocity_unit: String,
        velocity_to_angstrom_per_picosecond: f64,
        force_unit: String,
        force_to_kilojoule_per_mole_angstrom: f64,
    ) -> PyResult<Self> {
        let inner = molframe::traj::H5mdUnitSystem {
            length_unit,
            length_to_angstrom,
            time_unit,
            time_to_picosecond,
            velocity_unit,
            velocity_to_angstrom_per_picosecond,
            force_unit,
            force_to_kilojoule_per_mole_angstrom,
        };
        if !(length_to_angstrom.is_finite()
            && length_to_angstrom > 0.0
            && time_to_picosecond.is_finite()
            && time_to_picosecond > 0.0
            && velocity_to_angstrom_per_picosecond.is_finite()
            && velocity_to_angstrom_per_picosecond > 0.0
            && force_to_kilojoule_per_mole_angstrom.is_finite()
            && force_to_kilojoule_per_mole_angstrom > 0.0)
        {
            return Err(H5mdError::new_err(
                "H5MD unit scales must be finite and positive",
            ));
        }
        Ok(Self { inner })
    }

    #[staticmethod]
    fn canonical() -> Self {
        Self {
            inner: molframe::traj::H5mdUnitSystem::canonical(),
        }
    }

    #[getter]
    fn length_unit(&self) -> &str {
        &self.inner.length_unit
    }
    #[getter]
    const fn length_to_angstrom(&self) -> f64 {
        self.inner.length_to_angstrom
    }
    #[getter]
    fn time_unit(&self) -> &str {
        &self.inner.time_unit
    }
    #[getter]
    const fn time_to_picosecond(&self) -> f64 {
        self.inner.time_to_picosecond
    }
    #[getter]
    fn velocity_unit(&self) -> &str {
        &self.inner.velocity_unit
    }
    #[getter]
    const fn velocity_to_angstrom_per_picosecond(&self) -> f64 {
        self.inner.velocity_to_angstrom_per_picosecond
    }
    #[getter]
    fn force_unit(&self) -> &str {
        &self.inner.force_unit
    }
    #[getter]
    const fn force_to_kilojoule_per_mole_angstrom(&self) -> f64 {
        self.inner.force_to_kilojoule_per_mole_angstrom
    }
}

#[pyclass(name = "H5mdOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyH5mdOptions {
    pub(crate) inner: molframe::traj::H5mdOptions,
}

#[pymethods]
impl PyH5mdOptions {
    #[new]
    #[pyo3(signature = (*, particle_group="all", units=None))]
    fn new(particle_group: &str, units: Option<PyRef<'_, PyH5mdUnitSystem>>) -> Self {
        let defaults = molframe::traj::H5mdOptions::default();
        Self {
            inner: molframe::traj::H5mdOptions {
                particle_group: particle_group.to_owned(),
                units: match units {
                    Some(value) => value.inner.clone(),
                    None => defaults.units,
                },
            },
        }
    }

    #[getter]
    fn particle_group(&self) -> &str {
        &self.inner.particle_group
    }

    #[getter]
    fn units(&self) -> PyH5mdUnitSystem {
        PyH5mdUnitSystem {
            inner: self.inner.units.clone(),
        }
    }
}

impl From<PyH5mdOptions> for molframe::traj::H5mdOptions {
    fn from(value: PyH5mdOptions) -> Self {
        value.inner
    }
}

#[pyclass(name = "H5mdMetadata", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyH5mdMetadata {
    #[pyo3(get)]
    pub(crate) options: PyH5mdOptions,
    #[pyo3(get)]
    pub(crate) creator_name: String,
    #[pyo3(get)]
    pub(crate) creator_version: String,
}

impl From<molframe::traj::H5mdMetadata> for PyH5mdMetadata {
    fn from(value: molframe::traj::H5mdMetadata) -> Self {
        Self {
            options: PyH5mdOptions {
                inner: value.options,
            },
            creator_name: value.creator_name,
            creator_version: value.creator_version,
        }
    }
}

impl From<PyH5mdMetadata> for molframe::traj::H5mdMetadata {
    fn from(value: PyH5mdMetadata) -> Self {
        Self {
            options: value.options.inner,
            creator_name: value.creator_name,
            creator_version: value.creator_version,
        }
    }
}

#[pyclass(name = "H5mdTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyH5mdTrajectory {
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
    #[pyo3(get)]
    pub(crate) metadata: PyH5mdMetadata,
}

impl TryFrom<molframe::traj::H5mdTrajectory> for PyH5mdTrajectory {
    type Error = PyErr;

    fn try_from(value: molframe::traj::H5mdTrajectory) -> Result<Self, Self::Error> {
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
pub(crate) fn parse_h5md(py: Python<'_>, bytes: Vec<u8>) -> PyResult<Vec<PyTimestep>> {
    py.detach(move || -> PyResult<Vec<PyTimestep>> {
        molframe::traj::parse_h5md(&bytes)
            .map_err(|error| H5mdError::new_err(error.to_string()))?
            .into_iter()
            .map(TryInto::try_into)
            .collect()
    })
}

#[pyfunction]
pub(crate) fn parse_h5md_with_options(
    py: Python<'_>,
    bytes: Vec<u8>,
    options: PyH5mdOptions,
) -> PyResult<Vec<PyTimestep>> {
    py.detach(move || -> PyResult<Vec<PyTimestep>> {
        let options = options.inner;
        molframe::traj::parse_h5md_with_options(&bytes, &options)
            .map_err(|error| H5mdError::new_err(error.to_string()))?
            .into_iter()
            .map(TryInto::try_into)
            .collect()
    })
}

#[pyfunction]
pub(crate) fn parse_h5md_record_with_options(
    py: Python<'_>,
    bytes: Vec<u8>,
    options: PyH5mdOptions,
) -> PyResult<PyH5mdTrajectory> {
    py.detach(move || -> PyResult<PyH5mdTrajectory> {
        let options = options.inner;
        molframe::traj::parse_h5md_record_with_options(&bytes, &options)
            .map_err(|error| H5mdError::new_err(error.to_string()))?
            .try_into()
    })
}

#[pyfunction]
pub(crate) fn write_h5md(
    py: Python<'_>,
    frames: Vec<PyTimestep>,
    options: PyH5mdOptions,
) -> PyResult<Vec<u8>> {
    py.detach(move || -> PyResult<Vec<u8>> {
        let frames = frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        molframe::traj::write_h5md(&frames, &options.inner)
            .map_err(|error| H5mdError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_h5md_with_metadata(
    py: Python<'_>,
    frames: Vec<PyTimestep>,
    metadata: PyH5mdMetadata,
) -> PyResult<Vec<u8>> {
    py.detach(move || -> PyResult<Vec<u8>> {
        let frames = frames
            .into_iter()
            .map(TryInto::try_into)
            .collect::<PyResult<Vec<_>>>()?;
        molframe::traj::write_h5md_with_metadata(&frames, &metadata.into())
            .map_err(|error| H5mdError::new_err(error.to_string()))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("H5mdError", module.py().get_type::<H5mdError>())?;
    module.add_class::<PyH5mdUnitSystem>()?;
    module.add_class::<PyH5mdOptions>()?;
    module.add_class::<PyH5mdMetadata>()?;
    module.add_class::<PyH5mdTrajectory>()?;
    module.add_function(wrap_pyfunction!(parse_h5md, module)?)?;
    module.add_function(wrap_pyfunction!(parse_h5md_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(parse_h5md_record_with_options, module)?)?;
    module.add_function(wrap_pyfunction!(write_h5md, module)?)?;
    module.add_function(wrap_pyfunction!(write_h5md_with_metadata, module)?)?;
    Ok(())
}
