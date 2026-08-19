//! `PyO3` operation classes for declarative trajectory analyses.

use super::model::PyTrajectoryOperation;
use super::{execute_operation, explain_operation, native, operation_to_dict};
use crate::contract::PyAnalysis;
use crate::intrinsic::PyFrameAlignment;
use crate::query::PyAnalysisPolicy;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict};

#[pyclass(name = "RmsdToReference", frozen)]
pub(crate) struct PyRmsdToReference {
    pub(crate) operation: PyTrajectoryOperation,
}

#[pymethods]
impl PyRmsdToReference {
    #[new]
    #[pyo3(signature = (*, frames, reference, alignment=None, policy=None))]
    fn new(
        py: Python<'_>,
        frames: Py<PyAny>,
        reference: usize,
        alignment: Option<PyFrameAlignment>,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        let operation = PyTrajectoryOperation::RmsdToReference {
            frames,
            reference,
            alignment: match alignment {
                Some(alignment) => alignment,
                None => PyFrameAlignment::Unaligned,
            },
            policy: native::policy_or_default(policy),
        };
        native::validate_operation(py, &operation)?;
        Ok(Self { operation })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<PyAnalysis> {
        execute_operation(py, &self.operation)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = native::from_dict(config)?;
        if matches!(operation, PyTrajectoryOperation::RmsdToReference { .. }) {
            Ok(Self { operation })
        } else {
            Err(expected("trajectory_rmsd_to_reference"))
        }
    }

    fn __repr__(&self) -> String {
        "RmsdToReference(frames=<array>, reference=<index>, alignment=<FrameAlignment>)".to_owned()
    }
}

#[pyclass(name = "MeanSquaredDisplacementOp", frozen)]
pub(crate) struct PyMeanSquaredDisplacementOperation {
    pub(crate) operation: PyTrajectoryOperation,
}

#[pymethods]
impl PyMeanSquaredDisplacementOperation {
    #[new]
    #[pyo3(signature = (*, frames, maximum_lag, atoms=None, policy=None))]
    fn new(
        py: Python<'_>,
        frames: Py<PyAny>,
        maximum_lag: usize,
        atoms: Option<Py<PyAny>>,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        let operation = PyTrajectoryOperation::MeanSquaredDisplacement {
            frames,
            atoms,
            maximum_lag,
            policy: native::policy_or_default(policy),
        };
        native::validate_operation(py, &operation)?;
        Ok(Self { operation })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<PyAnalysis> {
        execute_operation(py, &self.operation)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = native::from_dict(config)?;
        if matches!(
            operation,
            PyTrajectoryOperation::MeanSquaredDisplacement { .. }
        ) {
            Ok(Self { operation })
        } else {
            Err(expected("mean_squared_displacement"))
        }
    }

    fn __repr__(&self) -> String {
        "MeanSquaredDisplacementOp(frames=<array>, maximum_lag=<int>, atoms=<array-or-none>)"
            .to_owned()
    }
}

#[pyclass(name = "PairwiseFittedRmsd", frozen)]
pub(crate) struct PyPairwiseFittedRmsd {
    pub(crate) operation: PyTrajectoryOperation,
}

#[pymethods]
impl PyPairwiseFittedRmsd {
    #[new]
    #[pyo3(signature = (*, frames, memory_limit=None, policy=None))]
    fn new(
        py: Python<'_>,
        frames: Py<PyAny>,
        memory_limit: Option<usize>,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        let operation = PyTrajectoryOperation::PairwiseFittedRmsd {
            frames,
            memory_limit: match memory_limit {
                Some(memory_limit) => memory_limit,
                None => pdbiox::traj::DEFAULT_PAIRWISE_MEMORY_LIMIT,
            },
            policy: native::policy_or_default(policy),
        };
        native::validate_operation(py, &operation)?;
        Ok(Self { operation })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<PyAnalysis> {
        execute_operation(py, &self.operation)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = native::from_dict(config)?;
        if matches!(operation, PyTrajectoryOperation::PairwiseFittedRmsd { .. }) {
            Ok(Self { operation })
        } else {
            Err(expected("pairwise_fitted_rmsd"))
        }
    }

    fn __repr__(&self) -> String {
        "PairwiseFittedRmsd(frames=<array>, memory_limit=<bytes>)".to_owned()
    }
}

#[pyclass(name = "GeneralizedProcrustesMean", frozen)]
pub(crate) struct PyGeneralizedProcrustesMean {
    pub(crate) operation: PyTrajectoryOperation,
}

#[pymethods]
impl PyGeneralizedProcrustesMean {
    #[new]
    #[pyo3(signature = (*, frames, tolerance=None, maximum_iterations=None, policy=None))]
    fn new(
        py: Python<'_>,
        frames: Py<PyAny>,
        tolerance: Option<f64>,
        maximum_iterations: Option<usize>,
        policy: Option<PyAnalysisPolicy>,
    ) -> PyResult<Self> {
        let operation = PyTrajectoryOperation::GeneralizedProcrustesMean {
            frames,
            tolerance: match tolerance {
                Some(tolerance) => tolerance,
                None => pdbiox::traj::DEFAULT_PROCRUSTES_TOLERANCE,
            },
            maximum_iterations: match maximum_iterations {
                Some(maximum_iterations) => maximum_iterations,
                None => pdbiox::traj::DEFAULT_PROCRUSTES_MAXIMUM_ITERATIONS,
            },
            policy: native::policy_or_default(policy),
        };
        native::validate_operation(py, &operation)?;
        Ok(Self { operation })
    }

    fn execute(&self, py: Python<'_>) -> PyResult<PyAnalysis> {
        execute_operation(py, &self.operation)
    }

    fn explain<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        explain_operation(py, &self.operation)
    }

    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        operation_to_dict(py, &self.operation)
    }

    #[staticmethod]
    fn from_dict(config: &Bound<'_, PyDict>) -> PyResult<Self> {
        let operation = native::from_dict(config)?;
        if matches!(
            operation,
            PyTrajectoryOperation::GeneralizedProcrustesMean { .. }
        ) {
            Ok(Self { operation })
        } else {
            Err(expected("generalized_procrustes_mean"))
        }
    }

    fn __repr__(&self) -> String {
        "GeneralizedProcrustesMean(frames=<array>, tolerance=<float>, maximum_iterations=<int>)"
            .to_owned()
    }
}

pub(super) fn typed_classes(value: &Bound<'_, PyAny>) -> PyResult<Option<PyTrajectoryOperation>> {
    macro_rules! extract {
        ($type:ty) => {
            if let Ok(value) = value.extract::<PyRef<'_, $type>>() {
                return Ok(Some(value.operation.clone_ref(value.py())));
            }
        };
    }
    extract!(PyRmsdToReference);
    extract!(PyMeanSquaredDisplacementOperation);
    extract!(PyPairwiseFittedRmsd);
    extract!(PyGeneralizedProcrustesMean);
    Ok(None)
}

pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyRmsdToReference>()?;
    module.add_class::<PyMeanSquaredDisplacementOperation>()?;
    module.add_class::<PyPairwiseFittedRmsd>()?;
    module.add_class::<PyGeneralizedProcrustesMean>()?;
    Ok(())
}

fn expected(tag: &str) -> PyErr {
    PyValueError::new_err(format!("expected serialized operation {tag:?}"))
}

#[cfg(test)]
#[path = "classes_tests.rs"]
mod tests;
