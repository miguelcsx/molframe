//! Mechanical Arrow `PyCapsule` exports over the Rust table adapters.

use crate::atom::PyAtoms;
use crate::bonds::PyBonds;
use crate::hierarchy::{PyChains, PyResidues};
use pyo3::exceptions::{PyNotImplementedError, PyRuntimeError};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;

#[pymethods]
impl PyAtoms {
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        capsule(
            py,
            pdbiox::AtomArrowTable::new(self.structure()).arrow_stream(),
        )
    }
}

#[pymethods]
impl PyResidues {
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        if !self.is_full_table() {
            return Err(PyNotImplementedError::new_err(
                "Arrow export of a filtered residue collection is not supported",
            ));
        }
        capsule(
            py,
            pdbiox::ResidueArrowTable::new(self.structure()).arrow_stream(),
        )
    }
}

#[pymethods]
impl PyChains {
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        if !self.is_full_table() {
            return Err(PyNotImplementedError::new_err(
                "Arrow export of a filtered chain collection is not supported",
            ));
        }
        capsule(
            py,
            pdbiox::ChainArrowTable::new(self.structure()).arrow_stream(),
        )
    }
}

#[pymethods]
impl PyBonds {
    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        capsule(
            py,
            pdbiox::BondArrowTable::new(self.structure()).arrow_stream(),
        )
    }
}

pub(crate) fn capsule<E: std::fmt::Display>(
    py: Python<'_>,
    stream: Result<pdbiox::ArrowStream, E>,
) -> PyResult<Bound<'_, PyCapsule>> {
    let stream = stream.map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    PyCapsule::new_with_value(py, stream.into_ffi(), c"arrow_array_stream")
}

#[cfg(test)]
#[path = "arrow_tests.rs"]
mod tests;
