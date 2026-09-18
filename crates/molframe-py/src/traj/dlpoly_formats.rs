//! Native bindings for `DL_POLY` CONFIG and HISTORY records.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;

pyo3::create_exception!(_native, DlPolyError, pyo3::exceptions::PyValueError);

#[pyclass(name = "DlPolyAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDlPolyAtom {
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) index: Option<u32>,
    #[pyo3(get)]
    pub(crate) mass: Option<f64>,
    #[pyo3(get)]
    pub(crate) charge: Option<f64>,
}

#[pymethods]
impl PyDlPolyAtom {
    #[new]
    #[pyo3(signature = (name, index=None, mass=None, charge=None))]
    fn new(name: String, index: Option<u32>, mass: Option<f64>, charge: Option<f64>) -> Self {
        Self {
            name,
            index,
            mass,
            charge,
        }
    }
}

impl From<PyDlPolyAtom> for molframe::traj::DlPolyAtom {
    fn from(value: PyDlPolyAtom) -> Self {
        Self {
            name: value.name.into(),
            index: value.index,
            mass: value.mass,
            charge: value.charge,
        }
    }
}

impl From<molframe::traj::DlPolyAtom> for PyDlPolyAtom {
    fn from(value: molframe::traj::DlPolyAtom) -> Self {
        Self {
            name: value.name.into(),
            index: value.index,
            mass: value.mass,
            charge: value.charge,
        }
    }
}

#[pyclass(name = "DlPolyFrame", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDlPolyFrame {
    #[pyo3(get)]
    pub(crate) step: i64,
    #[pyo3(get)]
    pub(crate) time: f64,
    #[pyo3(get)]
    pub(crate) positions: Vec<[f32; 3]>,
    #[pyo3(get)]
    pub(crate) velocities: Option<Vec<[f32; 3]>>,
    #[pyo3(get)]
    pub(crate) forces: Option<Vec<[f32; 3]>>,
    #[pyo3(get)]
    pub(crate) lattice_vectors: Option<[[f64; 3]; 3]>,
}

#[pymethods]
impl PyDlPolyFrame {
    #[new]
    #[pyo3(signature = (step, time, positions, velocities=None, forces=None, lattice_vectors=None))]
    fn new(
        step: i64,
        time: f64,
        positions: Vec<[f32; 3]>,
        velocities: Option<Vec<[f32; 3]>>,
        forces: Option<Vec<[f32; 3]>>,
        lattice_vectors: Option<[[f64; 3]; 3]>,
    ) -> PyResult<Self> {
        validate_streams(&positions, velocities.as_ref(), forces.as_ref())?;
        Ok(Self {
            step,
            time,
            positions,
            velocities,
            forces,
            lattice_vectors,
        })
    }

    fn to_timestep(&self, frame: usize) -> PyResult<PyTimestep> {
        native_frame(self.clone()).to_timestep(frame).try_into()
    }
}

impl From<PyDlPolyFrame> for molframe::traj::DlPolyFrame {
    fn from(value: PyDlPolyFrame) -> Self {
        Self {
            step: value.step,
            time: value.time,
            positions: value.positions,
            velocities: value.velocities,
            forces: value.forces,
            lattice_vectors: value.lattice_vectors,
        }
    }
}

impl From<molframe::traj::DlPolyFrame> for PyDlPolyFrame {
    fn from(value: molframe::traj::DlPolyFrame) -> Self {
        Self {
            step: value.step,
            time: value.time,
            positions: value.positions,
            velocities: value.velocities,
            forces: value.forces,
            lattice_vectors: value.lattice_vectors,
        }
    }
}

#[pyclass(name = "DlPolyConfig", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDlPolyConfig {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) level: u8,
    #[pyo3(get)]
    pub(crate) boundary: i32,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyDlPolyAtom>,
    #[pyo3(get)]
    pub(crate) frame: PyDlPolyFrame,
}

#[pymethods]
impl PyDlPolyConfig {
    #[new]
    fn new(
        title: String,
        level: u8,
        boundary: i32,
        atoms: Vec<PyDlPolyAtom>,
        frame: PyDlPolyFrame,
    ) -> Self {
        Self {
            title,
            level,
            boundary,
            atoms,
            frame,
        }
    }
}

impl From<PyDlPolyConfig> for molframe::traj::DlPolyConfig {
    fn from(value: PyDlPolyConfig) -> Self {
        Self {
            title: value.title.into(),
            level: value.level,
            boundary: value.boundary,
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            frame: value.frame.into(),
        }
    }
}

impl From<molframe::traj::DlPolyConfig> for PyDlPolyConfig {
    fn from(value: molframe::traj::DlPolyConfig) -> Self {
        Self {
            title: value.title.into(),
            level: value.level,
            boundary: value.boundary,
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            frame: value.frame.into(),
        }
    }
}

#[pyclass(name = "DlPolyHistory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDlPolyHistory {
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) level: u8,
    #[pyo3(get)]
    pub(crate) boundary: i32,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyDlPolyAtom>,
    #[pyo3(get)]
    pub(crate) frames: Vec<PyDlPolyFrame>,
}

#[pymethods]
impl PyDlPolyHistory {
    #[new]
    fn new(
        title: String,
        level: u8,
        boundary: i32,
        atoms: Vec<PyDlPolyAtom>,
        frames: Vec<PyDlPolyFrame>,
    ) -> Self {
        Self {
            title,
            level,
            boundary,
            atoms,
            frames,
        }
    }
}

impl From<PyDlPolyHistory> for molframe::traj::DlPolyHistory {
    fn from(value: PyDlPolyHistory) -> Self {
        Self {
            title: value.title.into(),
            level: value.level,
            boundary: value.boundary,
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            frames: value.frames.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<molframe::traj::DlPolyHistory> for PyDlPolyHistory {
    fn from(value: molframe::traj::DlPolyHistory) -> Self {
        Self {
            title: value.title.into(),
            level: value.level,
            boundary: value.boundary,
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            frames: value.frames.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyfunction]
pub(crate) fn parse_dlpoly_config(py: Python<'_>, text: &str) -> PyResult<PyDlPolyConfig> {
    py.detach(move || -> PyResult<PyDlPolyConfig> {
        molframe::traj::parse_dlpoly_config(text)
            .map(Into::into)
            .map_err(|error| DlPolyError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn parse_dlpoly_history(py: Python<'_>, text: &str) -> PyResult<PyDlPolyHistory> {
    py.detach(move || -> PyResult<PyDlPolyHistory> {
        molframe::traj::parse_dlpoly_history(text)
            .map(Into::into)
            .map_err(|error| DlPolyError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_dlpoly_config(py: Python<'_>, config: PyDlPolyConfig) -> PyResult<String> {
    py.detach(move || -> PyResult<String> {
        molframe::traj::write_dlpoly_config(&config.into())
            .map_err(|error| DlPolyError::new_err(error.to_string()))
    })
}

#[pyfunction]
pub(crate) fn write_dlpoly_history(py: Python<'_>, history: PyDlPolyHistory) -> PyResult<String> {
    py.detach(move || -> PyResult<String> {
        molframe::traj::write_dlpoly_history(&history.into())
            .map_err(|error| DlPolyError::new_err(error.to_string()))
    })
}

fn native_frame(frame: PyDlPolyFrame) -> molframe::traj::DlPolyFrame {
    frame.into()
}

fn validate_streams(
    positions: &[[f32; 3]],
    velocities: Option<&Vec<[f32; 3]>>,
    forces: Option<&Vec<[f32; 3]>>,
) -> PyResult<()> {
    if velocities.is_some_and(|values| values.len() != positions.len())
        || forces.is_some_and(|values| values.len() != positions.len())
    {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "DL_POLY streams must match the position atom count",
        ));
    }
    Ok(())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("DlPolyError", module.py().get_type::<DlPolyError>())?;
    module.add_class::<PyDlPolyAtom>()?;
    module.add_class::<PyDlPolyFrame>()?;
    module.add_class::<PyDlPolyConfig>()?;
    module.add_class::<PyDlPolyHistory>()?;
    module.add_function(wrap_pyfunction!(parse_dlpoly_config, module)?)?;
    module.add_function(wrap_pyfunction!(parse_dlpoly_history, module)?)?;
    module.add_function(wrap_pyfunction!(write_dlpoly_config, module)?)?;
    module.add_function(wrap_pyfunction!(write_dlpoly_history, module)?)?;
    Ok(())
}
