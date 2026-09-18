//! Single-structure and trajectory execution through one frame adapter.

use super::{FrameKernelResult, GovernedAnalysisError, StructureKernel};
use molframe_core::contract::{
    Analysis, AnalysisPolicy, AssemblyChoice, Coverage, MissingPolicy, ModelChoice, ParameterValue,
    Provenance, SourceRef, Status, SymmetryPolicy,
};
use molframe_core::index::ModelIndex;
use molframe_core::structure::Structure;
use molframe_core::{ExecutionContext, MemoryBudgetError};
use molframe_traj::{Frame, FrameAnalysis, Timestep, Trajectory, TrajectoryError, run_analysis};

/// An existing structure kernel bound to topology, policy and altloc selection.
#[derive(Debug)]
pub struct GovernedStructureAnalysis<'a, K> {
    template: Structure,
    source_atom_count: usize,
    selected_atoms: Vec<usize>,
    policy: AnalysisPolicy,
    kernel: &'a K,
    input_status: Status,
    input_warnings: Vec<molframe_core::diagnostic::Diagnostic>,
    input_assumptions: Vec<molframe_core::contract::Assumption>,
}

impl<'a, K: StructureKernel> GovernedStructureAnalysis<'a, K> {
    /// Resolves alternate conformations once and retains their source atom map.
    ///
    /// # Errors
    ///
    /// Returns every structural diagnostic if the selected topology cannot be
    /// materialised as a valid immutable structure.
    pub fn new(
        template: &Structure,
        policy: &AnalysisPolicy,
        kernel: &'a K,
    ) -> Result<Self, GovernedAnalysisError<K::Error>> {
        if !matches!(policy.assembly, AssemblyChoice::AsymmetricUnit) {
            return Err(GovernedAnalysisError::UnsupportedPolicyValue("assembly"));
        }
        if !matches!(policy.symmetry, SymmetryPolicy::None) {
            return Err(GovernedAnalysisError::UnsupportedPolicyValue("symmetry"));
        }
        let resolution = template.resolve_altlocs(policy);
        let selected_atoms = resolution.value.iter().map(|atom| atom as usize).collect();
        let selected = if resolution.value.len() == u64::from(template.atom_count()) {
            template.clone()
        } else {
            template
                .materialize(&resolution.value)
                .map_err(GovernedAnalysisError::InvalidStructure)?
        };
        let input_status = if resolution.coverage.ambiguous > 0 {
            combine_status(resolution.status, Status::Ambiguous)
        } else {
            resolution.status
        };
        Ok(Self {
            template: selected,
            source_atom_count: template.atom_count() as usize,
            selected_atoms,
            policy: policy.clone(),
            kernel,
            input_status,
            input_warnings: resolution.warnings,
            input_assumptions: resolution.assumptions,
        })
    }
}

impl<K: StructureKernel> FrameAnalysis for GovernedStructureAnalysis<'_, K> {
    type Output = Vec<K::Output>;
    type Partial = Vec<(usize, FrameKernelResult<K::Output>)>;
    type Error = GovernedAnalysisError<K::Error>;

    const NAME: &'static str = "governed-structure-kernel";
    const PARALLELIZABLE: bool = true;

    fn prepare(&self, _trajectory: &Trajectory) -> Result<Self::Partial, Self::Error> {
        Ok(Vec::new())
    }

    fn single_frame(
        &self,
        timestep: &Timestep,
        partial: &mut Self::Partial,
        context: &ExecutionContext,
    ) -> Result<(), Self::Error> {
        if timestep.positions.len() != self.source_atom_count {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: self.source_atom_count,
                found: timestep.positions.len(),
            }
            .into());
        }
        let mut editor = self
            .template
            .edit_coordinates(context)
            .map_err(|error| GovernedAnalysisError::Trajectory(memory_limit(error, context)))?;
        let Some(positions) = editor.positions_mut(ModelIndex::new(0)) else {
            return Err(TrajectoryError::AtomCountMismatch {
                expected: self.selected_atoms.len(),
                found: 0,
            }
            .into());
        };
        for (target, source) in positions.iter_mut().zip(&self.selected_atoms) {
            let Some(position) = timestep.positions.get(*source) else {
                return Err(TrajectoryError::SelectionOutOfRange {
                    index: *source,
                    atoms: timestep.positions.len(),
                }
                .into());
            };
            *target = *position;
        }
        let frame = editor
            .commit()
            .map_err(GovernedAnalysisError::InvalidStructure)?;
        let mut result = self
            .kernel
            .analyse_mapped(&frame, &self.policy, &self.selected_atoms, context)
            .map_err(GovernedAnalysisError::Kernel)?;
        enforce_missing_policy(&mut result, self.policy.missing_atoms)?;
        partial.push((timestep.frame, result));
        Ok(())
    }

    fn merge(&self, partials: Vec<Self::Partial>) -> Result<Self::Partial, Self::Error> {
        let mut merged: Self::Partial = partials.into_iter().flatten().collect();
        merged.sort_by_key(|(frame, _)| *frame);
        Ok(merged)
    }

    fn conclude(&self, partial: Self::Partial) -> Result<Analysis<Self::Output>, Self::Error> {
        combine_frames(self, partial)
    }
}

fn memory_limit(error: MemoryBudgetError, context: &ExecutionContext) -> TrajectoryError {
    let required = match error {
        MemoryBudgetError::Zero => 1,
        MemoryBudgetError::Exhausted { requested, .. } => requested,
    };
    TrajectoryError::MemoryLimit {
        required,
        limit: context.memory_budget().bytes(),
    }
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
    let positions = match policy.model {
        ModelChoice::First => structure.model_positions(ModelIndex::new(0)),
        ModelChoice::Index(index) => structure.model_positions(ModelIndex::new(index)),
        ModelChoice::All | ModelChoice::Ensemble => {
            return Err(GovernedAnalysisError::MultipleModelsRequested);
        }
        _ => return Err(GovernedAnalysisError::UnsupportedPolicyValue("model")),
    };
    let Some(positions) = positions else {
        return Err(TrajectoryError::AtomCountMismatch {
            expected: structure.atom_count() as usize,
            found: 0,
        }
        .into());
    };
    let trajectory = Trajectory::from_frames(vec![Frame {
        positions: positions.to_vec(),
    }])
    .map_err(TrajectoryError::from)?;
    let result = analyse_trajectory(structure, &trajectory, policy, kernel, context)?;
    let mut values = result.value;
    if values.len() != 1 {
        return Err(GovernedAnalysisError::MissingFrameOutput);
    }
    let value = values.remove(0);
    Ok(Analysis {
        value,
        status: result.status,
        coverage: result.coverage,
        warnings: result.warnings,
        assumptions: result.assumptions,
        provenance: result.provenance,
    })
}

/// Runs a structure kernel over every trajectory frame with fixed block order.
///
/// # Errors
///
/// Returns setup, kernel, missing-data or trajectory execution errors without
/// publishing a partial frame series.
pub fn analyse_trajectory<K: StructureKernel>(
    topology: &Structure,
    trajectory: &Trajectory,
    policy: &AnalysisPolicy,
    kernel: &K,
    context: &ExecutionContext,
) -> Result<Analysis<Vec<K::Output>>, GovernedAnalysisError<K::Error>> {
    let analysis = GovernedStructureAnalysis::new(topology, policy, kernel)?;
    run_analysis(trajectory, &analysis, context)
}

fn enforce_missing_policy<T, E>(
    result: &mut FrameKernelResult<T>,
    policy: MissingPolicy,
) -> Result<(), GovernedAnalysisError<E>> {
    validate_coverage(result.coverage)?;
    let incomplete = result.coverage.missing > 0 || result.coverage.ambiguous > 0;
    if !incomplete {
        return Ok(());
    }
    match policy {
        MissingPolicy::Ignore => Ok(()),
        MissingPolicy::Report => {
            result.status = combine_status(result.status, Status::Partial);
            Ok(())
        }
        MissingPolicy::Indeterminate => {
            result.status = Status::Indeterminate;
            Ok(())
        }
        MissingPolicy::Fail => Err(GovernedAnalysisError::MissingData {
            missing: result.coverage.missing,
            ambiguous: result.coverage.ambiguous,
        }),
        _ => Err(GovernedAnalysisError::UnsupportedPolicyValue(
            "missing_atoms",
        )),
    }
}

fn validate_coverage<E>(coverage: Coverage) -> Result<(), GovernedAnalysisError<E>> {
    let classified = coverage
        .used
        .checked_add(coverage.missing)
        .and_then(|count| count.checked_add(coverage.ambiguous));
    if classified == Some(coverage.intended) {
        Ok(())
    } else {
        Err(GovernedAnalysisError::InvalidCoverage {
            intended: coverage.intended,
            used: coverage.used,
            missing: coverage.missing,
            ambiguous: coverage.ambiguous,
        })
    }
}

fn combine_frames<K: StructureKernel>(
    analysis: &GovernedStructureAnalysis<'_, K>,
    partial: Vec<(usize, FrameKernelResult<K::Output>)>,
) -> Result<Analysis<Vec<K::Output>>, GovernedAnalysisError<K::Error>> {
    let frame_count =
        i64::try_from(partial.len()).map_err(|_| GovernedAnalysisError::CoverageOverflow)?;
    let mut result = Analysis {
        value: Vec::with_capacity(partial.len()),
        status: analysis.input_status,
        coverage: Coverage::default(),
        warnings: analysis.input_warnings.clone(),
        assumptions: analysis.input_assumptions.clone(),
        provenance: analysis
            .kernel
            .descriptor()
            .apply(Provenance::new(&analysis.policy).with_source(SourceRef::Memory))
            .with_parameter("frame_count", ParameterValue::Integer(frame_count)),
    };
    for (_, frame) in partial {
        result.status = combine_status(result.status, frame.status);
        result.coverage = add_coverage(result.coverage, frame.coverage)?;
        result.warnings.extend(frame.warnings);
        result.assumptions.extend(frame.assumptions);
        result.value.push(frame.value);
    }
    Ok(result)
}

fn add_coverage<E>(left: Coverage, right: Coverage) -> Result<Coverage, GovernedAnalysisError<E>> {
    Ok(Coverage {
        intended: left
            .intended
            .checked_add(right.intended)
            .ok_or(GovernedAnalysisError::CoverageOverflow)?,
        used: left
            .used
            .checked_add(right.used)
            .ok_or(GovernedAnalysisError::CoverageOverflow)?,
        missing: left
            .missing
            .checked_add(right.missing)
            .ok_or(GovernedAnalysisError::CoverageOverflow)?,
        ambiguous: left
            .ambiguous
            .checked_add(right.ambiguous)
            .ok_or(GovernedAnalysisError::CoverageOverflow)?,
    })
}

fn combine_status(left: Status, right: Status) -> Status {
    if matches!(left, Status::Indeterminate) || matches!(right, Status::Indeterminate) {
        Status::Indeterminate
    } else if matches!(left, Status::Ambiguous) || matches!(right, Status::Ambiguous) {
        Status::Ambiguous
    } else if matches!(left, Status::Partial) || matches!(right, Status::Partial) {
        Status::Partial
    } else {
        Status::Complete
    }
}

#[cfg(test)]
#[path = "execute_tests.rs"]
mod tests;
