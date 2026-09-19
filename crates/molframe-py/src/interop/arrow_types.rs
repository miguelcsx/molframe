//! Python handles for the public Rust Arrow table adapters.

use crate::arrow::capsule;
use crate::structure::PyStructure;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyCapsule};
use std::sync::{Mutex, MutexGuard};

#[pyclass(name = "ArrowStream", skip_from_py_object)]
pub(crate) struct PyArrowStream {
    inner: Mutex<Option<molframe::ArrowStream>>,
}

#[pymethods]
impl PyArrowStream {
    #[getter]
    fn consumed(&self) -> PyResult<bool> {
        Ok(self.lock()?.is_none())
    }

    #[pyo3(signature = (_requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        _requested_schema: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let stream = self
            .lock()?
            .take()
            .ok_or_else(|| PyRuntimeError::new_err("Arrow stream has already been consumed"))?;
        capsule(py, Ok::<_, std::convert::Infallible>(stream))
    }

    fn __repr__(&self) -> PyResult<&'static str> {
        if self.lock()?.is_some() {
            Ok("ArrowStream(consumed=False)")
        } else {
            Ok("ArrowStream(consumed=True)")
        }
    }
}

impl PyArrowStream {
    pub(crate) fn from_native(inner: molframe::ArrowStream) -> Self {
        Self {
            inner: Mutex::new(Some(inner)),
        }
    }

    fn lock(&self) -> PyResult<MutexGuard<'_, Option<molframe::ArrowStream>>> {
        self.inner
            .lock()
            .map_err(|_| PyRuntimeError::new_err("Arrow stream state is poisoned"))
    }
}

macro_rules! table_binding {
    ($rust:ident, $python:literal, $native:ty) => {
        #[pyclass(name = $python, skip_from_py_object)]
        pub(crate) struct $rust {
            inner: $native,
        }

        #[pymethods]
        impl $rust {
            #[new]
            fn new(structure: &PyStructure) -> Self {
                Self {
                    inner: <$native>::new(structure.structure()),
                }
            }

            fn schema(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                import_reader(py, native_stream(self.inner.arrow_stream())?)
                    .map(|reader| reader.getattr("schema").map(Bound::unbind))?
            }

            fn record_batches(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
                let reader = import_reader(py, native_stream(self.inner.arrow_stream())?)?;
                py.import("builtins")?
                    .getattr("list")?
                    .call1((reader,))
                    .map(Bound::unbind)
            }

            fn arrow_stream(&self) -> PyResult<PyArrowStream> {
                self.inner
                    .arrow_stream()
                    .map(PyArrowStream::from_native)
                    .map_err(|error| PyRuntimeError::new_err(error.to_string()))
            }

            fn __repr__(&self) -> &'static str {
                concat!($python, "(snapshot='retained')")
            }
        }
    };
}

table_binding!(PyAtomTable, "AtomArrowTable", molframe::AtomTable);
table_binding!(PyResidueTable, "ResidueArrowTable", molframe::ResidueTable);
table_binding!(PyChainTable, "ChainArrowTable", molframe::ChainTable);
table_binding!(PyBondTable, "BondArrowTable", molframe::interop::BondTable);

fn import_reader(py: Python<'_>, stream: molframe::ArrowStream) -> PyResult<Bound<'_, PyAny>> {
    let capsule = capsule(py, Ok::<_, std::convert::Infallible>(stream))?;
    py.import("pyarrow")?
        .getattr("RecordBatchReader")?
        .getattr("_import_from_c")?
        .call1((capsule,))
}

fn native_stream<E: std::fmt::Display>(
    stream: Result<molframe::ArrowStream, E>,
) -> PyResult<molframe::ArrowStream> {
    stream.map_err(|error| PyRuntimeError::new_err(error.to_string()))
}

pub(crate) fn register_classes(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyArrowStream>()?;
    module.add_class::<PyAtomTable>()?;
    module.add_class::<PyResidueTable>()?;
    module.add_class::<PyChainTable>()?;
    module.add_class::<PyBondTable>()?;
    Ok(())
}
