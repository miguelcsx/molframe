//! Reusable Rust-compiled queries and compact selected-index results.

use super::PyAnalysisPolicy;
use crate::contract::PyDiagnostic;
use crate::contract::{PyAnalysis, analysis_with_value};
use crate::errors::read_error;
use crate::structure::PyStructure;
use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pdbiox::query::{LogicalPlan, PhysicalQuery};
use pdbiox::{AnalysisPolicy, AtomSelection, Groups, Query, QueryStructure};
use pyo3::exceptions::PyOverflowError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(name = "Selection", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySelection {
    pub(crate) inner: AtomSelection,
}

#[pyclass(name = "Evaluation", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyEvaluation {
    #[pyo3(get)]
    pub(crate) selection: PySelection,
    #[pyo3(get)]
    pub(crate) warnings: Vec<PyDiagnostic>,
}

#[pymethods]
impl PySelection {
    #[new]
    fn new(indices: Vec<u32>) -> Self {
        Self {
            inner: indices.into_iter().collect(),
        }
    }
    fn __len__(&self) -> PyResult<usize> {
        usize::try_from(self.inner.len())
            .map_err(|_| PyOverflowError::new_err("selection length exceeds Python limits"))
    }
    fn __contains__(&self, index: u32) -> bool {
        self.inner.contains(index)
    }
    #[getter]
    fn indices<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        Array1::from_iter(self.inner.iter()).into_pyarray(py)
    }

    fn intersect(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.intersect(&other.inner),
        }
    }

    fn union(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.union(&other.inner),
        }
    }

    fn difference(&self, other: &Self) -> Self {
        Self {
            inner: self.inner.difference(&other.inner),
        }
    }

    fn complement(&self, atom_count: u32) -> Self {
        Self {
            inner: self.inner.complement(atom_count),
        }
    }
}

#[pyclass(name = "Query", frozen, from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyQuery {
    pub(crate) inner: Query,
}

#[pyclass(name = "LogicalPlan", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyLogicalPlan {
    inner: LogicalPlan,
}

#[pyclass(name = "PhysicalQuery", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyPhysicalQuery {
    inner: PhysicalQuery,
    structure: pdbiox::Structure,
    policy: AnalysisPolicy,
}

#[pymethods]
impl PyQuery {
    #[new]
    fn new(py: Python<'_>, source: &str) -> PyResult<Self> {
        Query::compile(source)
            .map(|inner| Self { inner })
            .map_err(|findings| read_error(py, &findings))
    }

    #[staticmethod]
    fn from_builder(builder: &super::PyQueryBuilder) -> Self {
        Self {
            inner: Query::from_builder(builder.0.clone()),
        }
    }

    fn logical_plan(&self) -> PyLogicalPlan {
        PyLogicalPlan {
            inner: self.inner.logical_plan(),
        }
    }

    #[pyo3(signature = (structure, *, policy=None))]
    fn plan(&self, structure: &PyStructure, policy: Option<&PyAnalysisPolicy>) -> PyPhysicalQuery {
        let policy = policy.map_or_else(AnalysisPolicy::default, |value| value.inner.clone());
        PyPhysicalQuery {
            inner: self.inner.plan(structure.structure(), &policy),
            structure: structure.structure().clone(),
            policy,
        }
    }
}

#[pymethods]
impl PyLogicalPlan {
    fn __repr__(&self) -> String {
        format!("LogicalPlan({:?})", self.inner)
    }
}

#[pymethods]
impl PyPhysicalQuery {
    #[pyo3(signature = (*, groups=None))]
    fn execute(
        &self,
        py: Python<'_>,
        groups: Option<BTreeMap<String, PySelection>>,
    ) -> PyResult<PyEvaluation> {
        let groups = groups_from_python(groups);
        let query = self.inner.clone();
        let structure = self.structure.clone();
        let policy = self.policy.clone();
        py.detach(move || query.evaluate(&structure, &policy, &groups, None))
            .map(PyEvaluation::from)
            .map_err(|findings| read_error(py, &findings))
    }
}

#[pymethods]
impl PyStructure {
    fn resolve_altlocs(&self, py: Python<'_>, policy: &PyAnalysisPolicy) -> PyResult<PyAnalysis> {
        let analysis = self.structure().resolve_altlocs(&policy.inner);
        analysis_with_value(py, analysis, |py, selection| {
            Py::new(py, PySelection { inner: selection }).map(Py::into_any)
        })
    }

    #[pyo3(signature = (query, *, policy=None, groups=None))]
    fn select(
        &self,
        py: Python<'_>,
        query: &PyQuery,
        policy: Option<&PyAnalysisPolicy>,
        groups: Option<BTreeMap<String, PySelection>>,
    ) -> PyResult<PySelection> {
        let structure = self.structure().clone();
        let query = query.inner.clone();
        let policy = policy.map_or_else(AnalysisPolicy::default, |value| value.inner.clone());
        let groups = groups_from_python(groups);
        py.detach(move || structure.select(&query, &policy, &groups))
            .map(|evaluation| PySelection {
                inner: evaluation.selection,
            })
            .map_err(|findings| read_error(py, &findings))
    }

    #[pyo3(signature = (query, *, policy=None, groups=None))]
    fn evaluate(
        &self,
        py: Python<'_>,
        query: &PyQuery,
        policy: Option<&PyAnalysisPolicy>,
        groups: Option<BTreeMap<String, PySelection>>,
    ) -> PyResult<PyEvaluation> {
        let structure = self.structure().clone();
        let query = query.inner.clone();
        let policy = policy.map_or_else(AnalysisPolicy::default, |value| value.inner.clone());
        let groups = groups_from_python(groups);
        py.detach(move || query.evaluate(&structure, &policy, &groups, None))
            .map(PyEvaluation::from)
            .map_err(|findings| read_error(py, &findings))
    }
}

impl From<pdbiox::Evaluation> for PyEvaluation {
    fn from(value: pdbiox::Evaluation) -> Self {
        Self {
            selection: PySelection {
                inner: value.selection,
            },
            warnings: value.warnings.into_iter().map(Into::into).collect(),
        }
    }
}

pub(crate) fn groups_from_python(groups: Option<BTreeMap<String, PySelection>>) -> Groups {
    let groups = match groups {
        Some(groups) => groups,
        None => BTreeMap::new(),
    };
    groups
        .into_iter()
        .map(|(name, selection)| (name.into(), selection.inner))
        .collect()
}

#[cfg(test)]
#[path = "selection_tests.rs"]
mod tests;
