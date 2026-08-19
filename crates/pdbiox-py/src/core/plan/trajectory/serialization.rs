//! Reproducible configuration and execution explanations for trajectory nodes.

use super::model::PyTrajectoryOperation;
use crate::query::PyAnalysisPolicy;
use pyo3::prelude::*;
use pyo3::types::PyDict;

pub(crate) fn explain_operation<'py>(
    py: Python<'py>,
    operation: &PyTrajectoryOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation.tag())?;
    result.set_item("execution", "rust")?;
    result.set_item("requires_c_contiguous", true)?;
    result.set_item("input_materializes", false)?;
    match operation {
        PyTrajectoryOperation::RmsdToReference {
            reference,
            alignment,
            policy,
            ..
        } => {
            result.set_item("complexity", "O(frames × atoms)")?;
            result.set_item("reference", reference)?;
            result.set_item("alignment", Py::new(py, *alignment)?)?;
            result.set_item("result", "NumPy float64 vector")?;
            policy_fields(py, &result, policy)?;
        }
        PyTrajectoryOperation::MeanSquaredDisplacement {
            atoms,
            maximum_lag,
            policy,
            ..
        } => {
            result.set_item("complexity", "O(maximum_lag × frames × selected_atoms)")?;
            result.set_item("maximum_lag", maximum_lag)?;
            result.set_item("atom_indices", atoms.is_some())?;
            result.set_item("result", "typed NumPy MSD series")?;
            policy_fields(py, &result, policy)?;
        }
        PyTrajectoryOperation::PairwiseFittedRmsd {
            memory_limit,
            policy,
            ..
        } => {
            result.set_item("complexity", "O(frames² × atoms) time and output")?;
            result.set_item("memory_limit", memory_limit)?;
            result.set_item("result", "bounded NumPy float64 matrix")?;
            policy_fields(py, &result, policy)?;
        }
        PyTrajectoryOperation::GeneralizedProcrustesMean {
            tolerance,
            maximum_iterations,
            policy,
            ..
        } => {
            result.set_item("complexity", "O(iterations × frames × atoms)")?;
            result.set_item("tolerance", tolerance)?;
            result.set_item("maximum_iterations", maximum_iterations)?;
            result.set_item("working_buffers", "mean and per-atom accumulator")?;
            result.set_item("result", "NumPy float32 coordinates")?;
            policy_fields(py, &result, policy)?;
        }
    }
    Ok(result)
}

pub(crate) fn operation_to_dict<'py>(
    py: Python<'py>,
    operation: &PyTrajectoryOperation,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("operation", operation.tag())?;
    match operation {
        PyTrajectoryOperation::RmsdToReference {
            frames,
            reference,
            alignment,
            policy,
        } => {
            result.set_item("frames", frames.bind(py))?;
            result.set_item("reference", reference)?;
            result.set_item("alignment", Py::new(py, *alignment)?)?;
            policy_fields(py, &result, policy)?;
        }
        PyTrajectoryOperation::MeanSquaredDisplacement {
            frames,
            atoms,
            maximum_lag,
            policy,
        } => {
            result.set_item("frames", frames.bind(py))?;
            result.set_item(
                "atoms",
                atoms
                    .as_ref()
                    .map_or_else(|| py.None(), |value| value.clone_ref(py)),
            )?;
            result.set_item("maximum_lag", maximum_lag)?;
            policy_fields(py, &result, policy)?;
        }
        PyTrajectoryOperation::PairwiseFittedRmsd {
            frames,
            memory_limit,
            policy,
        } => {
            result.set_item("frames", frames.bind(py))?;
            result.set_item("memory_limit", memory_limit)?;
            policy_fields(py, &result, policy)?;
        }
        PyTrajectoryOperation::GeneralizedProcrustesMean {
            frames,
            tolerance,
            maximum_iterations,
            policy,
        } => {
            result.set_item("frames", frames.bind(py))?;
            result.set_item("tolerance", tolerance)?;
            result.set_item("maximum_iterations", maximum_iterations)?;
            policy_fields(py, &result, policy)?;
        }
    }
    Ok(result)
}

fn policy_fields(
    py: Python<'_>,
    result: &Bound<'_, PyDict>,
    policy: &pdbiox::AnalysisPolicy,
) -> PyResult<()> {
    result.set_item(
        "policy",
        Py::new(py, PyAnalysisPolicy::from(policy.clone()))?,
    )?;
    result.set_item("policy_fingerprint", policy.fingerprint().to_string())
}
