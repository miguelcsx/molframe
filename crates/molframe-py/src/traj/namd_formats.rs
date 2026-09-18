//! Native bindings for NAMD binary coordinate snapshots.

use super::reader_types::PyTimestep;
use pyo3::prelude::*;

pyo3::create_exception!(_native, NamdError, pyo3::exceptions::PyValueError);

#[pyclass(name = "NamdEndian", frozen, eq, eq_int, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum PyNamdEndian {
    Little,
    Big,
}

impl From<PyNamdEndian> for molframe::traj::NamdEndian {
    fn from(value: PyNamdEndian) -> Self {
        match value {
            PyNamdEndian::Little => Self::Little,
            PyNamdEndian::Big => Self::Big,
        }
    }
}

impl From<molframe::traj::NamdEndian> for PyNamdEndian {
    fn from(value: molframe::traj::NamdEndian) -> Self {
        match value {
            molframe::traj::NamdEndian::Little => Self::Little,
            molframe::traj::NamdEndian::Big => Self::Big,
        }
    }
}

#[pyclass(name = "NamdBinary", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyNamdBinary {
    #[pyo3(get)]
    pub(crate) frame: PyTimestep,
    #[pyo3(get)]
    pub(crate) endian: PyNamdEndian,
}

impl TryFrom<molframe::traj::NamdBinary> for PyNamdBinary {
    type Error = PyErr;

    fn try_from(value: molframe::traj::NamdBinary) -> Result<Self, Self::Error> {
        Ok(Self {
            frame: value.frame.try_into()?,
            endian: value.endian.into(),
        })
    }
}

#[pyfunction]
pub(crate) fn parse_namd_binary(py: Python<'_>, bytes: Vec<u8>) -> PyResult<PyNamdBinary> {
    py.detach(move || -> PyResult<PyNamdBinary> {
        molframe::traj::parse_namd_binary(&bytes)
            .map_err(|error| NamdError::new_err(error.to_string()))?
            .try_into()
    })
}

#[pyfunction]
pub(crate) fn write_namd_binary(
    py: Python<'_>,
    frame: PyTimestep,
    endian: PyNamdEndian,
) -> PyResult<Vec<u8>> {
    py.detach(move || -> PyResult<Vec<u8>> {
        let frame = frame.try_into()?;
        molframe::traj::write_namd_binary(&frame, endian.into())
            .map_err(|error| NamdError::new_err(error.to_string()))
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("NamdError", module.py().get_type::<NamdError>())?;
    module.add_class::<PyNamdEndian>()?;
    module.add_class::<PyNamdBinary>()?;
    module.add_function(wrap_pyfunction!(parse_namd_binary, module)?)?;
    module.add_function(wrap_pyfunction!(write_namd_binary, module)?)?;
    Ok(())
}
