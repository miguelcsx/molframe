//! Typed trajectory requests executed over borrowed contiguous frame buffers.

use super::plan::inputs::FrameInput;
use super::plan::value::{ExecutionPlanError, PlanOperation};
use molframe_core::contract::{Analysis, AnalysisPolicy};
use molframe_traj::{
    EnsembleDistanceMatrix, FrameAlignment, FrameView, MeanSquaredDisplacement,
    analyse_generalized_procrustes_mean_view, analyse_mean_squared_displacement_view,
    analyse_pairwise_fitted_rmsd_view, analyse_rmsd_to_reference_view,
};

/// A reusable trajectory analysis over one borrowed frame-array slot.
///
/// Every variant borrows a C-contiguous `(frames, atoms, 3)` input. It never
/// repacks frames into per-frame vectors; only the documented result and
/// algorithmic working buffers are allocated.
#[derive(Clone, Debug)]
pub enum TrajectoryRequest {
    /// RMSD from each frame to one reference frame in `O(frames * atoms)`.
    RmsdToReference {
        /// Borrowed frame-array slot.
        frames: usize,
        /// Reference-frame index.
        reference: usize,
        /// Coordinate-frame alignment policy.
        alignment: FrameAlignment,
        /// Analysis contract policy.
        policy: AnalysisPolicy,
    },
    /// Window-averaged MSD in `O(maximum_lag * frames * selected_atoms)`.
    MeanSquaredDisplacement {
        /// Borrowed frame-array slot.
        frames: usize,
        /// Optional borrowed atom-index slot; `None` means all atoms.
        atoms: Option<usize>,
        /// Largest included frame lag.
        maximum_lag: usize,
        /// Analysis contract policy.
        policy: AnalysisPolicy,
    },
    /// Dense fitted RMSD matrix in `O(frames^2 * atoms)`.
    PairwiseFittedRmsd {
        /// Borrowed frame-array slot.
        frames: usize,
        /// Explicit upper bound for the dense result allocation.
        memory_limit: usize,
        /// Analysis contract policy.
        policy: AnalysisPolicy,
    },
    /// Generalized Procrustes mean in `O(iterations * frames * atoms)`.
    GeneralizedProcrustesMean {
        /// Borrowed frame-array slot.
        frames: usize,
        /// Positive convergence tolerance.
        tolerance: f64,
        /// Positive iteration limit.
        maximum_iterations: usize,
        /// Analysis contract policy.
        policy: AnalysisPolicy,
    },
}

impl TrajectoryRequest {
    /// Borrowed frame-array slot consumed by this request.
    #[must_use]
    pub const fn frame_slot(&self) -> usize {
        match self {
            Self::RmsdToReference { frames, .. }
            | Self::MeanSquaredDisplacement { frames, .. }
            | Self::PairwiseFittedRmsd { frames, .. }
            | Self::GeneralizedProcrustesMean { frames, .. } => *frames,
        }
    }

    /// Optional borrowed atom-index slot consumed by this request.
    #[must_use]
    pub const fn atom_slot(&self) -> Option<usize> {
        match self {
            Self::MeanSquaredDisplacement { atoms, .. } => *atoms,
            Self::RmsdToReference { .. }
            | Self::PairwiseFittedRmsd { .. }
            | Self::GeneralizedProcrustesMean { .. } => None,
        }
    }

    /// Analysis contract policy retained by this request.
    #[must_use]
    pub const fn policy(&self) -> &AnalysisPolicy {
        match self {
            Self::RmsdToReference { policy, .. }
            | Self::MeanSquaredDisplacement { policy, .. }
            | Self::PairwiseFittedRmsd { policy, .. }
            | Self::GeneralizedProcrustesMean { policy, .. } => policy,
        }
    }
}

impl From<TrajectoryRequest> for PlanOperation {
    fn from(value: TrajectoryRequest) -> Self {
        Self::Trajectory(Box::new(value))
    }
}

/// Typed trajectory values emitted by [`TrajectoryRequest`].
#[derive(Clone, Debug)]
pub enum TrajectoryValue {
    /// One RMSD value per input frame.
    RmsdToReference(Analysis<Vec<f64>>),
    /// One window-averaged displacement result per lag.
    MeanSquaredDisplacement(Analysis<Vec<MeanSquaredDisplacement>>),
    /// Dense fitted RMSD matrix under an explicit allocation ceiling.
    PairwiseFittedRmsd(Analysis<EnsembleDistanceMatrix>),
    /// Cartesian generalized Procrustes mean coordinates.
    GeneralizedProcrustesMean(Analysis<Vec<[f32; 3]>>),
}

/// Executes one trajectory request without materialising the borrowed frames.
pub(super) fn execute(
    request: &TrajectoryRequest,
    input: &FrameInput<'_>,
    atoms: &[usize],
) -> Result<TrajectoryValue, ExecutionPlanError> {
    let frames = FrameView::new(input.positions, input.frame_count, input.atom_count)?;
    match request {
        TrajectoryRequest::RmsdToReference {
            reference,
            alignment,
            policy,
            ..
        } => analyse_rmsd_to_reference_view(frames, *reference, *alignment, policy)
            .map(TrajectoryValue::RmsdToReference)
            .map_err(Into::into),
        TrajectoryRequest::MeanSquaredDisplacement {
            maximum_lag,
            policy,
            ..
        } => analyse_mean_squared_displacement_view(frames, atoms, *maximum_lag, policy)
            .map(TrajectoryValue::MeanSquaredDisplacement)
            .map_err(Into::into),
        TrajectoryRequest::PairwiseFittedRmsd {
            memory_limit,
            policy,
            ..
        } => analyse_pairwise_fitted_rmsd_view(frames, *memory_limit, policy)
            .map(TrajectoryValue::PairwiseFittedRmsd)
            .map_err(Into::into),
        TrajectoryRequest::GeneralizedProcrustesMean {
            tolerance,
            maximum_iterations,
            policy,
            ..
        } => analyse_generalized_procrustes_mean_view(
            frames,
            *tolerance,
            *maximum_iterations,
            policy,
        )
        .map(TrajectoryValue::GeneralizedProcrustesMean)
        .map_err(Into::into),
    }
}
