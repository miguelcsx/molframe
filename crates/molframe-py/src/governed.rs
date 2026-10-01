//! Running a structure kernel under an explicit policy and wrapping its result.

use crate::analysis_result::PyAnalysis;
use crate::bindings::PyStructure;
use crate::policy::PyAnalysisPolicy;
use molframe::analysis::{StructureKernel, analyse_structure};
use pyo3::{exceptions::PyValueError, prelude::*};
use std::fmt::Debug;

/// The policy a call names, or the default profile.
pub(crate) fn policy_of(policy: Option<PyRef<'_, PyAnalysisPolicy>>) -> molframe::AnalysisPolicy {
    policy.map_or_else(molframe::AnalysisPolicy::default, |policy| policy.0.clone())
}

/// Runs `kernel` on the first model under `policy` and converts the value.
///
/// Alternate conformations are resolved once under the policy, so the result's
/// coverage, status and assumptions describe what was actually analysed.
pub(crate) fn run<K>(
    py: Python<'_>,
    structure: &PyStructure,
    policy: &molframe::AnalysisPolicy,
    kernel: &K,
    convert: impl FnOnce(Python<'_>, K::Output) -> PyResult<Py<PyAny>>,
) -> PyResult<PyAnalysis>
where
    K: StructureKernel,
    K::Error: Debug,
{
    let structure = structure.inner.clone();
    let analysis = py
        .detach(|| {
            analyse_structure(
                structure.engine(),
                policy,
                kernel,
                &molframe::ExecutionContext::default(),
            )
        })
        .map_err(|error| PyValueError::new_err(format!("{error:?}")))?;
    let envelope = PyAnalysis::new(&analysis, py.None());
    let value = convert(py, analysis.value)?;
    Ok(envelope.with_value(value))
}
