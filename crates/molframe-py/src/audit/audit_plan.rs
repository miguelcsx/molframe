use super::{PyPolicyDimension, PyPolicyField};
use crate::query::PyAnalysisPolicy;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(name = "PolicySpace", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPolicySpace(pub(crate) molframe::PolicySpace);

#[pymethods]
impl PyPolicySpace {
    #[new]
    #[pyo3(signature = (policy=None))]
    pub(crate) fn new(policy: Option<&PyAnalysisPolicy>) -> Self {
        Self(molframe::PolicySpace::new(
            policy.map_or_else(molframe::AnalysisPolicy::default, |value| {
                value.inner.clone()
            }),
        ))
    }
    pub(crate) fn vary(&self, dimension: &PyPolicyDimension) -> Self {
        Self(self.0.clone().vary(dimension.0.clone()))
    }
    fn with_max_runs(&self, max_runs: usize) -> Self {
        Self(self.0.clone().with_max_runs(max_runs))
    }
    pub(crate) fn cost(&self) -> PyResult<usize> {
        self.0
            .cost()
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }
    pub(crate) fn plan(&self) -> PyResult<PyAuditPlan> {
        self.0
            .clone()
            .plan()
            .map(PyAuditPlan)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }
}

#[pyclass(name = "AuditPlan", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyAuditPlan(pub(crate) molframe::AuditPlan);

#[pymethods]
impl PyAuditPlan {
    #[getter]
    pub(crate) fn cost(&self) -> usize {
        self.0.cost()
    }
    #[getter]
    pub(crate) fn fields(&self) -> Vec<PyPolicyField> {
        self.0.fields().iter().copied().map(Into::into).collect()
    }
    #[getter]
    fn policies(&self) -> Vec<PyAnalysisPolicy> {
        self.0
            .policies()
            .iter()
            .cloned()
            .map(|inner| PyAnalysisPolicy { inner })
            .collect()
    }
}
