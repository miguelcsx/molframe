//! One structure kernel over the models of a dense ensemble, without a caller loop.
//!
//! The models of a dense ensemble share one topology, so each is a frame of the
//! same system. The frames are gathered once, charged against the execution
//! budget before they are allocated, and handed to the trajectory executor, which
//! is what already gives every kernel its coverage, status and provenance.

use super::execute::analyse_trajectory;
use super::{GovernedAnalysisError, StructureKernel};
use molframe_core::ExecutionContext;
use molframe_core::contract::{Analysis, AnalysisPolicy, ModelChoice, Outcome};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::index::ModelIndex;
use molframe_core::structure::Structure;
use molframe_traj::{Frame, Trajectory, TrajectoryError};

/// The finding for a model position the structure does not have.
pub(super) fn missing_model(structure: &Structure, model: u32) -> Diagnostic {
    Diagnostic::new(Code::E6003)
        .with_message(format!(
            "model {model} does not exist: the structure has {} models, addressed by position from 0",
            structure.model_count()
        ))
        .with_context("model", model.to_string())
        .with_context("models", structure.model_count().to_string())
}

/// Runs one selected model through the same adapter used for trajectories.
///
/// # Errors
///
/// Returns typed setup, kernel, missing-data or execution errors. Policies that
/// select all models must use [`analyse_trajectory`] so their output cardinality
/// remains explicit.
pub fn analyse_structure<K: StructureKernel>(
    structure: &Structure,
    policy: &AnalysisPolicy,
    kernel: &K,
    context: &ExecutionContext,
) -> Result<Analysis<K::Output>, GovernedAnalysisError<K::Error>> {
    structure
        .require_resolved_model()
        .map_err(|finding| GovernedAnalysisError::InvalidStructure(vec![finding]))?;
    let model = match policy.model {
        ModelChoice::First => 0,
        ModelChoice::Index(index) => index,
        ModelChoice::All | ModelChoice::Ensemble => {
            return Err(GovernedAnalysisError::MultipleModelsRequested);
        }
        _ => return Err(GovernedAnalysisError::UnsupportedPolicyValue("model")),
    };
    let Some(positions) = structure.model_positions(ModelIndex::new(model)) else {
        return Err(GovernedAnalysisError::InvalidStructure(vec![
            missing_model(structure, model),
        ]));
    };
    let trajectory = Trajectory::from_frames(vec![Frame {
        positions: positions.to_vec(),
    }])
    .map_err(TrajectoryError::from)?;
    let result = analyse_trajectory(structure, &trajectory, policy, kernel, context)?;
    let origin = result.atom_origin().map(<[u32]>::to_vec);
    let (outcome, quality, coverage, warnings, assumptions, provenance) = result.into_parts();
    let outcome = match outcome {
        Outcome::Determinate(mut values) => {
            if values.len() != 1 {
                return Err(GovernedAnalysisError::MissingFrameOutput);
            }
            Outcome::Determinate(values.remove(0))
        }
        Outcome::Indeterminate(reason) => Outcome::Indeterminate(reason),
    };
    let analysis = Analysis::from_parts(
        outcome,
        quality,
        coverage,
        warnings,
        assumptions,
        provenance,
    );
    Ok(match origin {
        Some(origin) => analysis.with_atom_origin(origin),
        None => analysis,
    })
}

/// Runs `kernel` on each requested model (every model when `models` is `None`).
///
/// The value is one result per model, in the order asked. Peak memory is one
/// copy of the selected coordinates, reserved from `context` up front: when the
/// selection does not fit, the error names the bytes needed and the caller
/// narrows `models` instead of the process running out of memory.
///
/// # Errors
///
/// A ragged ensemble, an absent model and an over-budget selection are refused
/// with typed findings; kernel and executor failures pass through as for
/// [`analyse_trajectory`].
pub fn analyse_models<K: StructureKernel>(
    structure: &Structure,
    models: Option<&[ModelIndex]>,
    policy: &AnalysisPolicy,
    kernel: &K,
    context: &ExecutionContext,
) -> Result<Analysis<Vec<K::Output>>, GovernedAnalysisError<K::Error>> {
    structure
        .require_resolved_model()
        .map_err(|finding| GovernedAnalysisError::InvalidStructure(vec![finding]))?;
    let count = structure.model_count();
    let chosen: Vec<ModelIndex> = match models {
        Some(chosen) => chosen.to_vec(),
        None => (0..count)
            .filter_map(|index| u32::try_from(index).ok().map(ModelIndex::new))
            .collect(),
    };
    let atoms = structure.positions().len();
    let bytes = chosen
        .len()
        .saturating_mul(atoms)
        .saturating_mul(std::mem::size_of::<[f32; 3]>());
    let _reservation = context
        .try_reserve(bytes)
        .map_err(|_| TrajectoryError::MemoryLimit {
            required: bytes,
            limit: context
                .memory_budget()
                .bytes()
                .saturating_sub(context.reserved_bytes()),
        })?;
    let mut frames = Vec::with_capacity(chosen.len());
    for model in &chosen {
        let Some(positions) = structure.model_positions(*model) else {
            return Err(GovernedAnalysisError::InvalidStructure(vec![
                missing_model(structure, model.get()),
            ]));
        };
        frames.push(Frame {
            positions: positions.to_vec(),
        });
    }
    let trajectory = Trajectory::from_frames(frames).map_err(TrajectoryError::from)?;
    analyse_trajectory(structure, &trajectory, policy, kernel, context)
}
