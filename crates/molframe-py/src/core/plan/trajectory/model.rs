//! Owned Python references and native options for trajectory plan nodes.

use crate::intrinsic::PyFrameAlignment;
use pyo3::prelude::*;
use pyo3::types::PyAny;

/// One typed trajectory operation retained by a Python plan.
#[derive(Debug)]
pub(crate) enum PyTrajectoryOperation {
    RmsdToReference {
        frames: Py<PyAny>,
        reference: usize,
        alignment: PyFrameAlignment,
        policy: molframe::AnalysisPolicy,
    },
    MeanSquaredDisplacement {
        frames: Py<PyAny>,
        atoms: Option<Py<PyAny>>,
        maximum_lag: usize,
        policy: molframe::AnalysisPolicy,
    },
    PairwiseFittedRmsd {
        frames: Py<PyAny>,
        memory_limit: usize,
        policy: molframe::AnalysisPolicy,
    },
    GeneralizedProcrustesMean {
        frames: Py<PyAny>,
        tolerance: f64,
        maximum_iterations: usize,
        policy: molframe::AnalysisPolicy,
    },
}

impl PyTrajectoryOperation {
    pub(crate) fn clone_ref(&self, py: Python<'_>) -> Self {
        match self {
            Self::RmsdToReference {
                frames,
                reference,
                alignment,
                policy,
            } => Self::RmsdToReference {
                frames: frames.clone_ref(py),
                reference: *reference,
                alignment: *alignment,
                policy: policy.clone(),
            },
            Self::MeanSquaredDisplacement {
                frames,
                atoms,
                maximum_lag,
                policy,
            } => Self::MeanSquaredDisplacement {
                frames: frames.clone_ref(py),
                atoms: atoms.as_ref().map(|value| value.clone_ref(py)),
                maximum_lag: *maximum_lag,
                policy: policy.clone(),
            },
            Self::PairwiseFittedRmsd {
                frames,
                memory_limit,
                policy,
            } => Self::PairwiseFittedRmsd {
                frames: frames.clone_ref(py),
                memory_limit: *memory_limit,
                policy: policy.clone(),
            },
            Self::GeneralizedProcrustesMean {
                frames,
                tolerance,
                maximum_iterations,
                policy,
            } => Self::GeneralizedProcrustesMean {
                frames: frames.clone_ref(py),
                tolerance: *tolerance,
                maximum_iterations: *maximum_iterations,
                policy: policy.clone(),
            },
        }
    }

    pub(crate) const fn tag(&self) -> &'static str {
        match self {
            Self::RmsdToReference { .. } => "trajectory_rmsd_to_reference",
            Self::MeanSquaredDisplacement { .. } => "mean_squared_displacement",
            Self::PairwiseFittedRmsd { .. } => "pairwise_fitted_rmsd",
            Self::GeneralizedProcrustesMean { .. } => "generalized_procrustes_mean",
        }
    }
}
