//! Python-owned audit results and callback execution.

use crate::audit::{PyAuditPlan, PyPolicyField};
use crate::query::PyAnalysisPolicy;
use pyo3::prelude::*;
use std::collections::{BTreeMap, BTreeSet};

#[pyclass(name = "AuditRun", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyAuditRun {
    #[pyo3(get)]
    policy: PyAnalysisPolicy,
    result: Py<PyAny>,
}

#[pymethods]
impl PyAuditRun {
    #[getter]
    fn result(&self, py: Python<'_>) -> Py<PyAny> {
        self.result.clone_ref(py)
    }
}

#[pyclass(name = "SensitiveItem", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PySensitiveItem {
    #[pyo3(get)]
    item: String,
    #[pyo3(get)]
    present_in: Vec<usize>,
}

#[pyclass(name = "DimensionSensitivity", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyDimensionSensitivity {
    #[pyo3(get)]
    field: PyPolicyField,
    #[pyo3(get)]
    sensitive_items: Vec<String>,
    #[pyo3(get)]
    mean_change: f64,
}

#[pyclass(name = "AuditReport", frozen, skip_from_py_object)]
#[derive(Debug)]
pub(crate) struct PyAuditReport {
    runs: Vec<PyAuditRun>,
    #[pyo3(get)]
    stability: f64,
    #[pyo3(get)]
    sensitive_items: Vec<PySensitiveItem>,
    #[pyo3(get)]
    dimensions: Vec<PyDimensionSensitivity>,
}

/// The pooled sensitivity of one policy field across a batch audit.
#[pyclass(name = "BatchDimension", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBatchDimension {
    #[pyo3(get)]
    field: PyPolicyField,
    #[pyo3(get)]
    mean_change: f64,
}

/// Aggregate stability statistics for a set of independently audited subjects.
#[pyclass(name = "BatchAudit", frozen, skip_from_py_object)]
#[derive(Clone, Debug)]
pub(crate) struct PyBatchAudit {
    #[pyo3(get)]
    subjects: usize,
    #[pyo3(get)]
    mean_stability: f64,
    #[pyo3(get)]
    min_stability: f64,
    #[pyo3(get)]
    max_stability: f64,
    #[pyo3(get)]
    dimensions: Vec<PyBatchDimension>,
}

#[pymethods]
impl PyAuditReport {
    #[getter]
    fn runs(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        self.runs
            .iter()
            .map(|run| {
                Py::new(
                    py,
                    PyAuditRun {
                        policy: run.policy.clone(),
                        result: run.result.clone_ref(py),
                    },
                )
                .map(Py::into_any)
            })
            .collect()
    }
}

#[pyfunction(name = "run_audit")]
pub(crate) fn audit(
    py: Python<'_>,
    plan: &PyAuditPlan,
    analyse: &Bound<'_, PyAny>,
    items: &Bound<'_, PyAny>,
) -> PyResult<PyAuditReport> {
    let mut runs = Vec::with_capacity(plan.0.cost());
    let mut sets = Vec::with_capacity(plan.0.cost());
    for policy in plan.0.policies() {
        let policy_object: Py<PyAnalysisPolicy> = Py::new(
            py,
            PyAnalysisPolicy {
                inner: policy.clone(),
            },
        )?;
        let result = analyse.call1((&policy_object,))?;
        let projected = items.call1((&result,))?;
        let mut set = BTreeSet::new();
        for value in projected.try_iter()? {
            set.insert(value?.str()?.to_str()?.to_owned());
        }
        sets.push(set);
        runs.push(PyAuditRun {
            policy: PyAnalysisPolicy {
                inner: policy.clone(),
            },
            result: result.unbind(),
        });
    }
    let frequencies = frequencies(&sets);
    let all_count = sets.len();
    let union_size = frequencies.len();
    let invariant = frequencies
        .values()
        .filter(|count| **count == all_count)
        .count();
    let stability = if union_size == 0 {
        1.0
    } else {
        count_to_f64(invariant) / count_to_f64(union_size)
    };
    let sensitive_items = frequencies
        .into_iter()
        .filter_map(|(item, count)| {
            (count != all_count).then(|| PySensitiveItem {
                present_in: sets
                    .iter()
                    .enumerate()
                    .filter_map(|(index, set)| set.contains(&item).then_some(index))
                    .collect(),
                item,
            })
        })
        .collect();
    let dimensions = plan
        .0
        .fields()
        .iter()
        .enumerate()
        .map(|(axis, field)| {
            let mut changed = BTreeSet::new();
            let mut loss = 0.0;
            let mut comparisons = 0usize;
            for first in 0..sets.len() {
                for second in first + 1..sets.len() {
                    let coordinates = plan.0.coordinates();
                    if coordinates[first].get(axis) == coordinates[second].get(axis)
                        || coordinates[first]
                            .iter()
                            .zip(&coordinates[second])
                            .enumerate()
                            .any(|(index, (left, right))| index != axis && left != right)
                    {
                        continue;
                    }
                    changed.extend(sets[first].symmetric_difference(&sets[second]).cloned());
                    let union = sets[first].union(&sets[second]).count();
                    if union > 0 {
                        loss += 1.0
                            - count_to_f64(sets[first].intersection(&sets[second]).count())
                                / count_to_f64(union);
                    }
                    comparisons += 1;
                }
            }
            PyDimensionSensitivity {
                field: (*field).into(),
                sensitive_items: changed.into_iter().collect(),
                mean_change: if comparisons == 0 {
                    0.0
                } else {
                    loss / count_to_f64(comparisons)
                },
            }
        })
        .collect();
    Ok(PyAuditReport {
        runs,
        stability,
        sensitive_items,
        dimensions,
    })
}

/// Executes one native batch audit, invoking caller callbacks once per policy
/// and subject.  Scientific data stays inside each callback's native kernel.
#[pyfunction]
pub(crate) fn audit_batch(
    plan: &PyAuditPlan,
    subjects: &Bound<'_, PyAny>,
    analyse: &Bound<'_, PyAny>,
    items: &Bound<'_, PyAny>,
) -> PyResult<PyBatchAudit> {
    let subjects = subjects
        .try_iter()?
        .map(|subject| subject.map(Bound::unbind))
        .collect::<PyResult<Vec<_>>>()?;
    let analyse = analyse.clone().unbind();
    let items = items.clone().unbind();
    let report = pdbiox::audit::audit_batch(
        &plan.0,
        &subjects,
        |subject, policy| {
            Python::attach(|py| {
                let policy = Py::new(
                    py,
                    PyAnalysisPolicy {
                        inner: policy.clone(),
                    },
                )?;
                let result = analyse.bind(py).call1((subject.bind(py), policy))?;
                let projected = items.bind(py).call1((&result,))?;
                let mut values = BTreeSet::new();
                for value in projected.try_iter()? {
                    values.insert(value?.str()?.to_str()?.to_owned());
                }
                Ok::<BTreeSet<String>, PyErr>(values)
            })
        },
        Clone::clone,
    )?;
    Ok(PyBatchAudit {
        subjects: report.subjects,
        mean_stability: report.mean_stability,
        min_stability: report.min_stability,
        max_stability: report.max_stability,
        dimensions: report
            .dimensions
            .into_iter()
            .map(|dimension| PyBatchDimension {
                field: dimension.field.into(),
                mean_change: dimension.mean_change,
            })
            .collect(),
    })
}

fn frequencies(sets: &[BTreeSet<String>]) -> BTreeMap<String, usize> {
    let mut values = BTreeMap::new();
    for set in sets {
        for item in set {
            *values.entry(item.clone()).or_insert(0) += 1;
        }
    }
    values
}

fn count_to_f64(value: usize) -> f64 {
    match value.to_string().parse() {
        Ok(value) => value,
        Err(_) => f64::INFINITY,
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyAuditRun>()?;
    module.add_class::<PySensitiveItem>()?;
    module.add_class::<PyDimensionSensitivity>()?;
    module.add_class::<PyAuditReport>()?;
    module.add_class::<PyBatchDimension>()?;
    module.add_class::<PyBatchAudit>()?;
    module.add_function(wrap_pyfunction!(audit, module)?)?;
    module.add_function(wrap_pyfunction!(audit_batch, module)?)?;
    Ok(())
}
