//! Exact provenance replay with one operation-level Python callback.

use crate::contract::PyProvenance;
use crate::errors::ReexecutionError;
use crate::query::PyAnalysisPolicy;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[pyclass(name = "ReexecutionEnvironment", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyReexecutionEnvironment {
    molframe: String,
    schema: Option<String>,
    component: Option<String>,
}

#[pymethods]
impl PyReexecutionEnvironment {
    #[new]
    #[pyo3(signature = (*, molframe_version=None, schema_version=None, component_version=None))]
    fn new(
        molframe_version: Option<String>,
        schema_version: Option<String>,
        component_version: Option<String>,
    ) -> Self {
        let current = molframe::core::contract::ReexecutionEnvironment::current();
        let molframe_version = match molframe_version {
            Some(value) => value,
            None => current.molframe_version.to_owned(),
        };
        Self {
            molframe: molframe_version,
            schema: schema_version,
            component: component_version,
        }
    }

    #[staticmethod]
    fn current() -> Self {
        Self::new(None, None, None)
    }

    #[getter]
    fn molframe_version(&self) -> &str {
        &self.molframe
    }

    #[getter]
    fn schema_version(&self) -> Option<&str> {
        self.schema.as_deref()
    }

    #[getter]
    fn component_version(&self) -> Option<&str> {
        self.component.as_deref()
    }
}

impl PyReexecutionEnvironment {
    fn native(&self) -> molframe::core::contract::ReexecutionEnvironment<'_> {
        molframe::core::contract::ReexecutionEnvironment {
            molframe_version: &self.molframe,
            schema_version: self.schema.as_deref(),
            component_version: self.component.as_deref(),
        }
    }
}

#[pyclass(name = "Reexecution", frozen, skip_from_py_object)]
pub(crate) struct PyReexecution {
    value: Py<PyAny>,
    provenance: PyProvenance,
}

#[pymethods]
impl PyReexecution {
    #[getter]
    fn value(&self, py: Python<'_>) -> Py<PyAny> {
        self.value.clone_ref(py)
    }

    #[getter]
    fn provenance(&self) -> PyProvenance {
        self.provenance.clone()
    }
}

#[pyfunction]
#[pyo3(signature = (provenance, input, environment, run))]
pub(crate) fn reexecute_from_provenance(
    py: Python<'_>,
    provenance: &PyProvenance,
    input: &Bound<'_, PyBytes>,
    environment: &PyReexecutionEnvironment,
    run: Py<PyAny>,
) -> PyResult<PyReexecution> {
    let bytes = input.as_bytes();
    let native = molframe::core::contract::reexecute_from_provenance(
        &provenance.inner,
        bytes,
        environment.native(),
        |bytes, policy| -> PyResult<Py<PyAny>> {
            let input = PyBytes::new(py, bytes);
            let policy = Py::new(py, PyAnalysisPolicy::from(policy.clone()))?;
            run.call1(py, (input, policy))
        },
    )
    .map_err(|error| ReexecutionError::new_err(error.to_string()))?;
    Ok(PyReexecution {
        value: native.value?,
        provenance: native.provenance.into(),
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyReexecutionEnvironment>()?;
    module.add_class::<PyReexecution>()?;
    module.add_function(wrap_pyfunction!(reexecute_from_provenance, module)?)?;
    Ok(())
}
