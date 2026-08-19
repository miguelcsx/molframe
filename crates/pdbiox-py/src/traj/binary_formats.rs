//! Native bindings for DCD, TRR and XTC trajectory containers.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;

pyo3::create_exception!(_native, DcdError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, TrrError, pyo3::exceptions::PyValueError);
pyo3::create_exception!(_native, XtcError, pyo3::exceptions::PyValueError);

#[pyclass(name = "DcdEndian", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyDcdEndian {
    Little,
    Big,
}

impl From<PyDcdEndian> for pdbiox::traj::DcdEndian {
    fn from(value: PyDcdEndian) -> Self {
        match value {
            PyDcdEndian::Little => Self::Little,
            PyDcdEndian::Big => Self::Big,
        }
    }
}

impl From<pdbiox::traj::DcdEndian> for PyDcdEndian {
    fn from(value: pdbiox::traj::DcdEndian) -> Self {
        match value {
            pdbiox::traj::DcdEndian::Little => Self::Little,
            pdbiox::traj::DcdEndian::Big => Self::Big,
        }
    }
}

#[pyclass(name = "DcdHeader", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDcdHeader {
    #[pyo3(get)]
    pub(crate) frame_count: usize,
    #[pyo3(get)]
    pub(crate) atom_count: usize,
    #[pyo3(get)]
    pub(crate) start_step: i32,
    #[pyo3(get)]
    pub(crate) save_interval: i32,
    #[pyo3(get)]
    pub(crate) delta_akma: f64,
    #[pyo3(get)]
    pub(crate) fixed_atom_count: usize,
    #[pyo3(get)]
    pub(crate) titles: Vec<String>,
    #[pyo3(get)]
    pub(crate) endian: PyDcdEndian,
    #[pyo3(get)]
    pub(crate) charmm: bool,
    #[pyo3(get)]
    pub(crate) has_unit_cell: bool,
    #[pyo3(get)]
    pub(crate) has_fourth_dimension: bool,
}

#[pymethods]
impl PyDcdHeader {
    #[new]
    #[pyo3(signature = (frame_count, atom_count, start_step, save_interval, delta_akma, fixed_atom_count, titles, endian, charmm, has_unit_cell, has_fourth_dimension))]
    fn new(
        frame_count: usize,
        atom_count: usize,
        start_step: i32,
        save_interval: i32,
        delta_akma: f64,
        fixed_atom_count: usize,
        titles: Vec<String>,
        endian: PyDcdEndian,
        charmm: bool,
        has_unit_cell: bool,
        has_fourth_dimension: bool,
    ) -> Self {
        Self {
            frame_count,
            atom_count,
            start_step,
            save_interval,
            delta_akma,
            fixed_atom_count,
            titles,
            endian,
            charmm,
            has_unit_cell,
            has_fourth_dimension,
        }
    }
}

impl From<pdbiox::traj::DcdHeader> for PyDcdHeader {
    fn from(value: pdbiox::traj::DcdHeader) -> Self {
        Self {
            frame_count: value.frame_count,
            atom_count: value.atom_count,
            start_step: value.start_step,
            save_interval: value.save_interval,
            delta_akma: value.delta_akma,
            fixed_atom_count: value.fixed_atom_count,
            titles: value.titles.into_iter().map(Into::into).collect(),
            endian: value.endian.into(),
            charmm: value.charmm,
            has_unit_cell: value.has_unit_cell,
            has_fourth_dimension: value.has_fourth_dimension,
        }
    }
}

impl From<PyDcdHeader> for pdbiox::traj::DcdHeader {
    fn from(value: PyDcdHeader) -> Self {
        Self {
            frame_count: value.frame_count,
            atom_count: value.atom_count,
            start_step: value.start_step,
            save_interval: value.save_interval,
            delta_akma: value.delta_akma,
            fixed_atom_count: value.fixed_atom_count,
            titles: value.titles.into_iter().map(Into::into).collect(),
            endian: value.endian.into(),
            charmm: value.charmm,
            has_unit_cell: value.has_unit_cell,
            has_fourth_dimension: value.has_fourth_dimension,
        }
    }
}

#[pyclass(name = "DcdTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDcdTrajectory {
    #[pyo3(get)]
    pub(crate) header: PyDcdHeader,
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
}

#[pymethods]
impl PyDcdTrajectory {
    #[new]
    fn new(header: PyDcdHeader, frames: Vec<PyTimestep>) -> Self {
        Self { header, frames }
    }
}

#[pyclass(name = "DcdWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDcdWriteOptions {
    #[pyo3(get)]
    pub(crate) endian: PyDcdEndian,
    #[pyo3(get)]
    pub(crate) title: String,
    #[pyo3(get)]
    pub(crate) start_step: i32,
    #[pyo3(get)]
    pub(crate) save_interval: i32,
    #[pyo3(get)]
    pub(crate) delta_akma: f32,
}

#[pymethods]
impl PyDcdWriteOptions {
    #[new]
    #[pyo3(signature = (*, endian=None, title=None, start_step=None, save_interval=None, delta_akma=None))]
    fn new(
        endian: Option<PyDcdEndian>,
        title: Option<String>,
        start_step: Option<i32>,
        save_interval: Option<i32>,
        delta_akma: Option<f32>,
    ) -> Self {
        let defaults = pdbiox::traj::DcdWriteOptions::default();
        Self {
            endian: match endian {
                Some(value) => value,
                None => defaults.endian.into(),
            },
            title: match title {
                Some(value) => value,
                None => defaults.title.into(),
            },
            start_step: match start_step {
                Some(value) => value,
                None => defaults.start_step,
            },
            save_interval: match save_interval {
                Some(value) => value,
                None => defaults.save_interval,
            },
            delta_akma: match delta_akma {
                Some(value) => value,
                None => defaults.delta_akma,
            },
        }
    }
}

impl From<PyDcdWriteOptions> for pdbiox::traj::DcdWriteOptions {
    fn from(value: PyDcdWriteOptions) -> Self {
        Self {
            endian: value.endian.into(),
            title: value.title.into(),
            start_step: value.start_step,
            save_interval: value.save_interval,
            delta_akma: value.delta_akma,
        }
    }
}

#[pyclass(name = "TrrPrecision", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyTrrPrecision {
    Single,
    Double,
}

impl From<PyTrrPrecision> for pdbiox::traj::TrrPrecision {
    fn from(value: PyTrrPrecision) -> Self {
        match value {
            PyTrrPrecision::Single => Self::Single,
            PyTrrPrecision::Double => Self::Double,
        }
    }
}

impl From<pdbiox::traj::TrrPrecision> for PyTrrPrecision {
    fn from(value: pdbiox::traj::TrrPrecision) -> Self {
        match value {
            pdbiox::traj::TrrPrecision::Single => Self::Single,
            pdbiox::traj::TrrPrecision::Double => Self::Double,
        }
    }
}

#[pyclass(name = "TrrTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTrrTrajectory {
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
    #[pyo3(get)]
    pub(crate) steps: Vec<i32>,
    #[pyo3(get)]
    pub(crate) precision: Vec<PyTrrPrecision>,
}

#[pyclass(name = "TrrWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyTrrWriteOptions {
    #[pyo3(get)]
    pub(crate) precision: PyTrrPrecision,
}

#[pymethods]
impl PyTrrWriteOptions {
    #[new]
    #[pyo3(signature = (precision=None))]
    fn new(precision: Option<PyTrrPrecision>) -> Self {
        Self {
            precision: match precision {
                Some(value) => value,
                None => PyTrrPrecision::Single,
            },
        }
    }
}

impl From<PyTrrWriteOptions> for pdbiox::traj::TrrWriteOptions {
    fn from(value: PyTrrWriteOptions) -> Self {
        Self {
            precision: value.precision.into(),
        }
    }
}

#[pyclass(name = "XtcTrajectory", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyXtcTrajectory {
    #[pyo3(get)]
    pub(crate) frames: Vec<PyTimestep>,
    #[pyo3(get)]
    pub(crate) steps: Vec<u32>,
    #[pyo3(get)]
    pub(crate) precision: Vec<f32>,
}

#[pyclass(name = "XtcWriteOptions", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyXtcWriteOptions {
    #[pyo3(get)]
    pub(crate) precision: f32,
}

#[pymethods]
impl PyXtcWriteOptions {
    #[new]
    #[pyo3(signature = (precision=None))]
    fn new(precision: Option<f32>) -> Self {
        Self {
            precision: match precision {
                Some(value) => value,
                None => 1_000.0,
            },
        }
    }
}

impl TryFrom<pdbiox::traj::DcdTrajectory> for PyDcdTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::DcdTrajectory) -> Result<Self, Self::Error> {
        Ok(Self {
            header: value.header.into(),
            frames: convert_frames(value.frames)?,
        })
    }
}

impl TryFrom<pdbiox::traj::TrrTrajectory> for PyTrrTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::TrrTrajectory) -> Result<Self, Self::Error> {
        Ok(Self {
            frames: convert_frames(value.frames)?,
            steps: value.steps,
            precision: value.precision.into_iter().map(Into::into).collect(),
        })
    }
}

impl TryFrom<pdbiox::traj::XtcTrajectory> for PyXtcTrajectory {
    type Error = PyErr;

    fn try_from(value: pdbiox::traj::XtcTrajectory) -> Result<Self, Self::Error> {
        Ok(Self {
            frames: convert_frames(value.frames)?,
            steps: value.steps,
            precision: value.precision,
        })
    }
}

#[pyfunction]
pub(crate) fn parse_dcd(bytes: Vec<u8>) -> PyResult<PyDcdTrajectory> {
    pdbiox::traj::parse_dcd(&bytes)
        .map_err(|error| DcdError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
pub(crate) fn write_dcd(frames: Vec<PyTimestep>, options: PyDcdWriteOptions) -> PyResult<Vec<u8>> {
    let frames = native_frames(frames)?;
    pdbiox::traj::write_dcd(&frames, &options.into())
        .map_err(|error| DcdError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn parse_trr(bytes: Vec<u8>) -> PyResult<PyTrrTrajectory> {
    pdbiox::traj::parse_trr(&bytes)
        .map_err(|error| TrrError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
pub(crate) fn write_trr(frames: Vec<PyTimestep>, options: PyTrrWriteOptions) -> PyResult<Vec<u8>> {
    let frames = native_frames(frames)?;
    pdbiox::traj::write_trr(&frames, options.into())
        .map_err(|error| TrrError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_trr_with_precisions(
    frames: Vec<PyTimestep>,
    precisions: Vec<PyTrrPrecision>,
) -> PyResult<Vec<u8>> {
    let frames = native_frames(frames)?;
    let precisions = precisions.into_iter().map(Into::into).collect::<Vec<_>>();
    pdbiox::traj::write_trr_with_precisions(&frames, &precisions)
        .map_err(|error| TrrError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn parse_xtc(bytes: Vec<u8>) -> PyResult<PyXtcTrajectory> {
    pdbiox::traj::parse_xtc(&bytes)
        .map_err(|error| XtcError::new_err(error.to_string()))?
        .try_into()
}

#[pyfunction]
pub(crate) fn write_xtc(frames: Vec<PyTimestep>, options: PyXtcWriteOptions) -> PyResult<Vec<u8>> {
    let frames = native_frames(frames)?;
    pdbiox::traj::write_xtc(
        &frames,
        pdbiox::traj::XtcWriteOptions {
            precision: options.precision,
        },
    )
    .map_err(|error| XtcError::new_err(error.to_string()))
}

#[pyfunction]
pub(crate) fn write_xtc_with_precisions(
    frames: Vec<PyTimestep>,
    precisions: Vec<f32>,
) -> PyResult<Vec<u8>> {
    let frames = native_frames(frames)?;
    pdbiox::traj::write_xtc_with_precisions(&frames, &precisions)
        .map_err(|error| XtcError::new_err(error.to_string()))
}

fn native_frames(frames: Vec<PyTimestep>) -> PyResult<Vec<pdbiox::traj::Timestep>> {
    frames.into_iter().map(TryInto::try_into).collect()
}

fn convert_frames(frames: Vec<pdbiox::traj::Timestep>) -> PyResult<Vec<PyTimestep>> {
    frames.into_iter().map(TryInto::try_into).collect()
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("DcdError", module.py().get_type::<DcdError>())?;
    module.add("TrrError", module.py().get_type::<TrrError>())?;
    module.add("XtcError", module.py().get_type::<XtcError>())?;
    module.add_class::<PyDcdEndian>()?;
    module.add_class::<PyDcdHeader>()?;
    module.add_class::<PyDcdTrajectory>()?;
    module.add_class::<PyDcdWriteOptions>()?;
    module.add_class::<PyTrrPrecision>()?;
    module.add_class::<PyTrrTrajectory>()?;
    module.add_class::<PyTrrWriteOptions>()?;
    module.add_class::<PyXtcTrajectory>()?;
    module.add_class::<PyXtcWriteOptions>()?;
    module.add_function(wrap_pyfunction!(parse_dcd, module)?)?;
    module.add_function(wrap_pyfunction!(write_dcd, module)?)?;
    module.add_function(wrap_pyfunction!(parse_trr, module)?)?;
    module.add_function(wrap_pyfunction!(write_trr, module)?)?;
    module.add_function(wrap_pyfunction!(write_trr_with_precisions, module)?)?;
    module.add_function(wrap_pyfunction!(parse_xtc, module)?)?;
    module.add_function(wrap_pyfunction!(write_xtc, module)?)?;
    module.add_function(wrap_pyfunction!(write_xtc_with_precisions, module)?)?;
    Ok(())
}
