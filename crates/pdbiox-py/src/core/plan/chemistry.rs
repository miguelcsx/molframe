//! Declarative chemistry operations backed by the facade plan.

use super::dispatch::require_operation_tag;
use super::{Operation, PyPlan, execute_native};
use crate::facade::{PyBondInference, PyBondInferenceReport};
use crate::graph::PySpatialBackend;
use crate::structure::PyStructure;
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

/// Infers covalent bonds from element radii and coordinates.
#[pyclass(name = "InferBonds", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyInferBonds {
    pub(crate) options: PyBondInference,
}

#[pymethods]
impl PyInferBonds {
    #[new]
    #[pyo3(signature = (*, options=None))]
    fn new(options: Option<PyBondInference>) -> Self {
        Self {
            options: match options {
                Some(options) => options,
                None => PyBondInference(pdbiox::BondInference::default()),
            },
        }
    }

    #[getter]
    fn options(&self) -> PyBondInference {
        self.options
    }

    fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyBondInferenceReport> {
        let plan = PyPlan {
            operations: vec![("result".to_owned(), Operation::BondInference(self.clone()))],
        };
        let result = execute_native(&plan, py, Some(structure), None)?;
        let Some(entry) = result.entries.into_iter().next() else {
            return Err(PyValueError::new_err(
                "native bond-inference plan returned no result",
            ));
        };
        match entry.value {
            pdbiox::PlanValue::BondInference(report) => {
                Ok(PyBondInferenceReport::from_native(*report))
            }
            _ => Err(PyValueError::new_err(
                "native bond-inference plan returned an incompatible result",
            )),
        }
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.explain_native(py)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.to_dict_native(py)
    }
}

impl PyInferBonds {
    pub(super) fn explain_native<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("operation", "infer_bonds")?;
        result.set_item("backend", format!("{:?}", self.options.0.backend))?;
        result.set_item("complexity", "O(n + candidate_pairs)")?;
        result.set_item("execution", "native")?;
        result.set_item("materializes", true)?;
        result.set_item("output", "structure_with_bonds")?;
        result.set_item("scale", self.options.0.scale)?;
        result.set_item("lower_bound", self.options.0.lower_bound)?;
        result.set_item(
            "exclude_across_chains",
            self.options.0.exclude_across_chains,
        )?;
        result.set_item("respect_existing", self.options.0.respect_existing)?;
        Ok(result)
    }

    pub(super) fn to_dict_native<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("operation", "infer_bonds")?;
        result.set_item("scale", self.options.0.scale)?;
        result.set_item("lower_bound", self.options.0.lower_bound)?;
        result.set_item(
            "exclude_across_chains",
            self.options.0.exclude_across_chains,
        )?;
        result.set_item("respect_existing", self.options.0.respect_existing)?;
        result.set_item("backend", PySpatialBackend::from(self.options.0.backend))?;
        Ok(result)
    }
}

#[pymethods]
impl PyInferBonds {
    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, "infer_bonds")?;
        let scale = match optional(config, "scale")? {
            Some(value) => value.extract()?,
            None => 1.15,
        };
        let lower_bound = match optional(config, "lower_bound")? {
            Some(value) => value.extract()?,
            None => 0.4,
        };
        let exclude_across_chains = match optional(config, "exclude_across_chains")? {
            Some(value) => value.extract()?,
            None => false,
        };
        let respect_existing = match optional(config, "respect_existing")? {
            Some(value) => value.extract()?,
            None => true,
        };
        let backend: PySpatialBackend = match optional(config, "backend")? {
            Some(value) => value.extract()?,
            None => PySpatialBackend::Auto,
        };
        Ok(Self {
            options: PyBondInference(pdbiox::BondInference {
                scale,
                lower_bound,
                exclude_across_chains,
                respect_existing,
                backend: backend.into(),
            }),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "InferBonds(scale={}, lower_bound={}, backend={:?})",
            self.options.0.scale, self.options.0.lower_bound, self.options.0.backend
        )
    }
}

pub(super) fn serialized_request(config: &Bound<'_, PyDict>) -> PyResult<Option<PyInferBonds>> {
    let Some(operation) = config.get_item("operation")? else {
        return Ok(None);
    };
    let operation: &str = operation.extract()?;
    if operation != "infer_bonds" {
        return Ok(None);
    }
    PyInferBonds::from_dict(config).map(Some)
}

fn optional<'py>(
    config: &'py Bound<'py, PyDict>,
    name: &str,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    config
        .get_item(name)
        .map_err(|error| PyKeyError::new_err(error.to_string()))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyInferBonds>()
}
