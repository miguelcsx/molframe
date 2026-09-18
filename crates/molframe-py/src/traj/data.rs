//! Typed projections for trajectory containers and source metadata.

use super::model::PyTrajectory;
use super::types::PyTrajectoryFormat;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "FormatMetadata", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyFormatMetadata {
    inner: molframe::traj::FormatMetadata,
}

#[pymethods]
impl PyFormatMetadata {
    #[staticmethod]
    fn none() -> Self {
        Self {
            inner: molframe::traj::FormatMetadata::None,
        }
    }

    #[getter]
    fn kind(&self) -> &'static str {
        match &self.inner {
            molframe::traj::FormatMetadata::None => "none",
            molframe::traj::FormatMetadata::Xtc { .. } => "xtc",
            molframe::traj::FormatMetadata::Trr { .. } => "trr",
            molframe::traj::FormatMetadata::Dcd(_) => "dcd",
            molframe::traj::FormatMetadata::AmberNetcdf(_) => "amber_netcdf",
            molframe::traj::FormatMetadata::Tng { .. } => "tng",
            molframe::traj::FormatMetadata::Gsd(_) => "gsd",
            molframe::traj::FormatMetadata::H5md(_) => "h5md",
            molframe::traj::FormatMetadata::Trz { .. } => "trz",
            molframe::traj::FormatMetadata::Namd(_) => "namd",
            molframe::traj::FormatMetadata::AmberRestart { .. } => "amber_restart",
            molframe::traj::FormatMetadata::AmberAscii { .. } => "amber_ascii",
            molframe::traj::FormatMetadata::Gro(_) => "gro",
            molframe::traj::FormatMetadata::Xyz(_) => "xyz",
            molframe::traj::FormatMetadata::Aims(_) => "aims",
            molframe::traj::FormatMetadata::Txyz(_) => "txyz",
            molframe::traj::FormatMetadata::DlPolyConfig(_) => "dlpoly_config",
            molframe::traj::FormatMetadata::DlPolyHistory(_) => "dlpoly_history",
            molframe::traj::FormatMetadata::CharmmCard(_) => "charmm_card",
            molframe::traj::FormatMetadata::Gamess(_) => "gamess",
            molframe::traj::FormatMetadata::Gromos11(_) => "gromos11",
            molframe::traj::FormatMetadata::Dms(_) => "dms",
            _ => "unknown",
        }
    }

    fn __repr__(&self) -> String {
        format!("FormatMetadata(kind={:?})", self.kind())
    }
}

#[pyclass(name = "TrajectoryMetadata", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTrajectoryMetadata {
    steps: Option<Vec<i64>>,
    format: PyFormatMetadata,
}

#[pymethods]
impl PyTrajectoryMetadata {
    #[getter]
    fn steps(&self) -> Option<Vec<i64>> {
        self.steps.clone()
    }

    #[getter]
    fn format(&self) -> PyFormatMetadata {
        self.format.clone()
    }
}

#[pyclass(name = "TrajectoryData", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyTrajectoryData {
    trajectory: PyTrajectory,
    format: Option<PyTrajectoryFormat>,
    metadata: PyTrajectoryMetadata,
}

#[pymethods]
impl PyTrajectoryData {
    #[new]
    fn new(py: Python<'_>, trajectory: PyRef<'_, PyTrajectory>) -> PyResult<Self> {
        Self::from_trajectory_value(py, (*trajectory).clone())
    }

    #[staticmethod]
    fn from_trajectory(py: Python<'_>, trajectory: PyRef<'_, PyTrajectory>) -> PyResult<Self> {
        Self::from_trajectory_value(py, (*trajectory).clone())
    }

    #[getter]
    fn trajectory(&self) -> PyTrajectory {
        self.trajectory.clone()
    }

    #[getter]
    fn format(&self) -> Option<PyTrajectoryFormat> {
        self.format
    }

    #[getter]
    fn metadata(&self) -> PyTrajectoryMetadata {
        self.metadata.clone()
    }

    fn select_frames(&self, py: Python<'_>, indices: Vec<usize>) -> PyResult<Self> {
        let trajectory = self
            .trajectory
            .select_indices(py, indices)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Self::from_trajectory_value(py, trajectory)
    }
}

impl PyTrajectoryData {
    fn from_trajectory_value(py: Python<'_>, trajectory: PyTrajectory) -> PyResult<Self> {
        let format = trajectory.python_format()?;
        Ok(Self {
            metadata: PyTrajectoryMetadata {
                steps: trajectory.source_steps(py)?,
                format: PyFormatMetadata {
                    inner: trajectory.source_format_metadata(),
                },
            },
            trajectory,
            format,
        })
    }
}

#[pyclass(name = "TrzWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTrzWriteOptions {
    #[pyo3(get)]
    pub(crate) title: String,
}

#[pymethods]
impl PyTrzWriteOptions {
    #[new]
    fn new(title: String) -> PyResult<Self> {
        if title.len() > 80 {
            return Err(PyValueError::new_err("TRZ title cannot exceed 80 bytes"));
        }
        Ok(Self { title })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyFormatMetadata>()?;
    module.add_class::<PyTrajectoryMetadata>()?;
    module.add_class::<PyTrajectoryData>()?;
    module.add_class::<PyTrzWriteOptions>()?;
    Ok(())
}
