//! Declarative query selection executed by the native facade plan.

use super::dispatch::require_operation_tag;
use super::{Operation, PyPlan, execute_native};
use crate::query::{PyAnalysisPolicy, PyEvaluation, PyQuery};
use crate::structure::PyStructure;
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

#[pyclass(name = "SelectQuery", frozen)]
pub(crate) struct PySelectQuery {
    pub(crate) query: Py<PyAny>,
    policy: Option<Py<PyAny>>,
    pub(crate) request: pdbiox::SelectionRequest,
}

impl PySelectQuery {
    pub(crate) fn clone_ref(&self, py: Python<'_>) -> Self {
        Self {
            query: self.query.clone_ref(py),
            policy: self.policy.as_ref().map(|value| value.clone_ref(py)),
            request: self.request.clone(),
        }
    }
}

pub(super) fn explain<'py>(
    py: Python<'py>,
    operation: &PySelectQuery,
) -> PyResult<Bound<'py, PyDict>> {
    operation.explain(py)
}

pub(super) fn to_dict<'py>(
    py: Python<'py>,
    operation: &PySelectQuery,
) -> PyResult<Bound<'py, PyDict>> {
    operation.to_dict(py)
}

#[pymethods]
impl PySelectQuery {
    #[new]
    #[pyo3(signature = (*, query, policy=None))]
    fn new(py: Python<'_>, query: Py<PyAny>, policy: Option<Py<PyAny>>) -> PyResult<Self> {
        let query_ref = query.bind(py).extract::<PyRef<'_, PyQuery>>()?;
        let policy_inner = policy
            .as_ref()
            .map(|value| {
                value
                    .bind(py)
                    .extract::<PyRef<'_, PyAnalysisPolicy>>()
                    .map(|policy| policy.inner.clone())
            })
            .transpose()?;
        let request_policy = match policy_inner {
            Some(policy) => policy,
            None => pdbiox::AnalysisPolicy::default(),
        };
        Ok(Self {
            query,
            policy,
            request: pdbiox::SelectionRequest::from_query(query_ref.inner.clone(), request_policy),
        })
    }

    fn execute(&self, py: Python<'_>, structure: &PyStructure) -> PyResult<PyEvaluation> {
        let plan = PyPlan {
            operations: vec![(
                "result".to_owned(),
                Operation::Selection(self.clone_ref(py)),
            )],
        };
        let result = execute_native(&plan, py, Some(structure))?;
        let Some(entry) = result.entries.into_iter().next() else {
            return Err(PyValueError::new_err(
                "native selection plan returned no result",
            ));
        };
        match entry.value {
            pdbiox::PlanValue::Selection(value) => Ok(PyEvaluation::from(*value)),
            _ => Err(PyValueError::new_err(
                "native selection plan returned an incompatible result",
            )),
        }
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("operation", "select")?;
        result.set_item("execution", "native")?;
        result.set_item("complexity", "planner-dependent")?;
        result.set_item("materializes", "selected atom indices")?;
        Ok(result)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let result = PyDict::new(py);
        result.set_item("operation", "select")?;
        result.set_item("query", self.query.bind(py))?;
        result.set_item(
            "policy",
            self.policy
                .as_ref()
                .map_or_else(|| py.None(), |value| value.clone_ref(py)),
        )?;
        Ok(result)
    }

    #[staticmethod]
    fn from_dict(py: Python<'_>, config: &Bound<'_, PyDict>) -> PyResult<Self> {
        require_operation_tag(config, "select")?;
        let query = required(config, "query")?.unbind();
        let policy = match config.get_item("policy")? {
            Some(value) if !value.is_none() => Some(value.unbind()),
            _ => None,
        };
        Self::new(py, query, policy)
    }

    fn __repr__(&self) -> String {
        "SelectQuery(query=<Query>, policy=<AnalysisPolicy-or-none>)".to_owned()
    }
}

pub(super) fn typed_request(value: &Bound<'_, PyAny>) -> Option<PySelectQuery> {
    if let Ok(value) = value.extract::<PyRef<'_, PySelectQuery>>() {
        return Some(value.clone_ref(value.py()));
    }
    None
}

pub(super) fn serialized_request(
    py: Python<'_>,
    config: &Bound<'_, PyDict>,
) -> PyResult<Option<PySelectQuery>> {
    let Some(value) = config.get_item("operation")? else {
        return Ok(None);
    };
    if value.extract::<&str>()? != "select" {
        return Ok(None);
    }
    PySelectQuery::from_dict(py, config).map(Some)
}

fn required<'py>(config: &'py Bound<'py, PyDict>, name: &str) -> PyResult<Bound<'py, PyAny>> {
    config
        .get_item(name)?
        .ok_or_else(|| PyKeyError::new_err(format!("missing operation field: {name}")))
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PySelectQuery>()
}
