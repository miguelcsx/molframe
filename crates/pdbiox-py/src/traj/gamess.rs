//! Native bindings for GAMESS coordinate-producing output.

use pyo3::prelude::*;

pyo3::create_exception!(_native, GamessError, pyo3::exceptions::PyValueError);

#[pyclass(name = "GamessRunType", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyGamessRunType {
    Optimize,
    Surface,
}

impl From<pdbiox::traj::GamessRunType> for PyGamessRunType {
    fn from(value: pdbiox::traj::GamessRunType) -> Self {
        match value {
            pdbiox::traj::GamessRunType::Optimize => Self::Optimize,
            pdbiox::traj::GamessRunType::Surface => Self::Surface,
        }
    }
}

impl From<PyGamessRunType> for pdbiox::traj::GamessRunType {
    fn from(value: PyGamessRunType) -> Self {
        match value {
            PyGamessRunType::Optimize => Self::Optimize,
            PyGamessRunType::Surface => Self::Surface,
        }
    }
}

#[pyclass(name = "GamessAtom", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGamessAtom {
    #[pyo3(get)]
    pub(crate) name: String,
    #[pyo3(get)]
    pub(crate) nuclear_charge: Option<f64>,
}

#[pymethods]
impl PyGamessAtom {
    #[new]
    #[pyo3(signature = (name, nuclear_charge=None))]
    fn new(name: String, nuclear_charge: Option<f64>) -> Self {
        Self {
            name,
            nuclear_charge,
        }
    }
}

impl From<pdbiox::traj::GamessAtom> for PyGamessAtom {
    fn from(value: pdbiox::traj::GamessAtom) -> Self {
        Self {
            name: value.name.into(),
            nuclear_charge: value.nuclear_charge,
        }
    }
}

impl From<PyGamessAtom> for pdbiox::traj::GamessAtom {
    fn from(value: PyGamessAtom) -> Self {
        Self {
            name: value.name.into(),
            nuclear_charge: value.nuclear_charge,
        }
    }
}

#[pyclass(name = "GamessFrame", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGamessFrame {
    #[pyo3(get)]
    pub(crate) step: i64,
    #[pyo3(get)]
    pub(crate) energy: Option<f64>,
    #[pyo3(get)]
    pub(crate) surface_coordinates: Option<[f64; 2]>,
    #[pyo3(get)]
    pub(crate) positions: Vec<[f32; 3]>,
}

#[pymethods]
impl PyGamessFrame {
    #[new]
    #[pyo3(signature = (step, positions, energy=None, surface_coordinates=None))]
    fn new(
        step: i64,
        positions: Vec<[f32; 3]>,
        energy: Option<f64>,
        surface_coordinates: Option<[f64; 2]>,
    ) -> Self {
        Self {
            step,
            energy,
            surface_coordinates,
            positions,
        }
    }
}

impl From<pdbiox::traj::GamessFrame> for PyGamessFrame {
    fn from(value: pdbiox::traj::GamessFrame) -> Self {
        Self {
            step: value.step,
            energy: value.energy,
            surface_coordinates: value.surface_coordinates,
            positions: value.positions,
        }
    }
}

impl From<PyGamessFrame> for pdbiox::traj::GamessFrame {
    fn from(value: PyGamessFrame) -> Self {
        Self {
            step: value.step,
            energy: value.energy,
            surface_coordinates: value.surface_coordinates,
            positions: value.positions,
        }
    }
}

#[pyclass(name = "GamessTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyGamessTrajectory {
    #[pyo3(get)]
    pub(crate) run_type: PyGamessRunType,
    #[pyo3(get)]
    pub(crate) atoms: Vec<PyGamessAtom>,
    #[pyo3(get)]
    pub(crate) frames: Vec<PyGamessFrame>,
}

#[pymethods]
impl PyGamessTrajectory {
    #[new]
    fn new(
        run_type: PyGamessRunType,
        atoms: Vec<PyGamessAtom>,
        frames: Vec<PyGamessFrame>,
    ) -> Self {
        Self {
            run_type,
            atoms,
            frames,
        }
    }
}

impl From<pdbiox::traj::GamessTrajectory> for PyGamessTrajectory {
    fn from(value: pdbiox::traj::GamessTrajectory) -> Self {
        Self {
            run_type: value.run_type.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            frames: value.frames.into_iter().map(Into::into).collect(),
        }
    }
}

impl From<PyGamessTrajectory> for pdbiox::traj::GamessTrajectory {
    fn from(value: PyGamessTrajectory) -> Self {
        Self {
            run_type: value.run_type.into(),
            atoms: value.atoms.into_iter().map(Into::into).collect(),
            frames: value.frames.into_iter().map(Into::into).collect(),
        }
    }
}

#[pyfunction]
pub(crate) fn parse_gamess_output(text: &str) -> PyResult<PyGamessTrajectory> {
    pdbiox::traj::parse_gamess_output(text)
        .map(Into::into)
        .map_err(|error| GamessError::new_err(error.to_string()))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("GamessError", module.py().get_type::<GamessError>())?;
    module.add_class::<PyGamessRunType>()?;
    module.add_class::<PyGamessAtom>()?;
    module.add_class::<PyGamessFrame>()?;
    module.add_class::<PyGamessTrajectory>()?;
    module.add_function(wrap_pyfunction!(parse_gamess_output, module)?)?;
    Ok(())
}
