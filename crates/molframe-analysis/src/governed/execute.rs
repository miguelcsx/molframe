//! Single-structure and trajectory execution through one frame adapter.

use super::requirements::without_hydrogens;
use super::system::AnalysisSystem;
use super::{FrameKernelResult, GovernedAnalysisError, StructureKernel};
use molframe_core::contract::{
    Analysis, AnalysisPolicy, Assumption, AssumptionSource, Coverage, HydrogenPolicy, Impact,
    Indeterminacy, MissingPolicy, MissingPolicyError, ModelChoice, Outcome, ParameterValue,
    PolicyField, Provenance, Quality, SourceRef, resolve_missing,
};
use molframe_core::index::ModelIndex;
use molframe_core::structure::Structure;
use molframe_core::{ExecutionContext, MemoryBudgetError};
use molframe_geom::Rigid;
use molframe_traj::{Frame, FrameAnalysis, Timestep, Trajectory, TrajectoryError, run_analysis};

/// The policy fields the executor itself applies to every kernel, before the kernel
/// runs: which system, which model, which conformations, which hydrogens, and what
/// to do about missing atoms.
const EXECUTOR_READS: [PolicyField; 6] = [
    PolicyField::Assembly,
    PolicyField::Symmetry,
    PolicyField::Model,
    PolicyField::Altloc,
    PolicyField::Hydrogens,
    PolicyField::MissingAtoms,
];

/// An existing structure kernel bound to topology, policy and altloc selection.
#[derive(Debug)]
pub struct GovernedStructureAnalysis<'a, K> {
    template: Structure,
    source_atom_count: usize,
    selected_atoms: Vec<usize>,
    motions: Vec<Option<Rigid>>,
    policy: AnalysisPolicy,
    kernel: &'a K,
    input_quality: Quality,
    input_indeterminacy: Option<Indeterminacy>,
    input_warnings: Vec<molframe_core::diagnostic::Diagnostic>,
    input_assumptions: Vec<Assumption>,
}

impl<'a, K: StructureKernel> GovernedStructureAnalysis<'a, K> {
    /// Builds the analysed system, resolves its alternate conformations and hydrogens
    /// once, and retains the source atom and motion of every atom it keeps.
    ///
    /// Equivalent to [`GovernedStructureAnalysis::new_in`] under a default context.
    ///
    /// # Errors
    ///
    /// See [`GovernedStructureAnalysis::new_in`].
    pub fn new(
        template: &Structure,
        policy: &AnalysisPolicy,
        kernel: &'a K,
    ) -> Result<Self, GovernedAnalysisError<K::Error>> {
        Self::new_in(template, policy, kernel, &ExecutionContext::default())
    }

    /// Builds the analysed system under an execution context.
    ///
    /// The policy decides what the system is. The asymmetric unit is the input;
    /// a biological assembly is its generated chain instances; crystal contacts
    /// are the unit plus every chain that touches it from a neighbouring cell.
    /// Alternate conformations and hydrogens are then resolved over that system,
    /// so a rule that picks one conformer picks it in every copy.
    ///
    /// # Errors
    ///
    /// Returns the findings that explain a contradictory policy or a system that
    /// cannot be built, an unsupported hydrogen rule, a kernel that cannot run over
    /// replicated atoms, or the diagnostics of a topology that cannot be materialised.
    pub fn new_in(
        template: &Structure,
        policy: &AnalysisPolicy,
        kernel: &'a K,
        context: &ExecutionContext,
    ) -> Result<Self, GovernedAnalysisError<K::Error>> {
        if matches!(policy.hydrogens, HydrogenPolicy::IncludeInferred) {
            // Adding hydrogens that were not modelled is not something this library does,
            // and silently using the modelled ones would be a different decision.
            return Err(GovernedAnalysisError::UnsupportedPolicyValue("hydrogens"));
        }
        if let Some(refused) = kernel
            .descriptor()
            .forbidden_resolutions()
            .iter()
            .find(|resolution| resolution.forbids(policy))
        {
            return Err(GovernedAnalysisError::ForbiddenResolution {
                analysis: kernel.descriptor().name().into(),
                field: refused.field().name(),
                reason: refused.reason(),
            });
        }
        let system = AnalysisSystem::build(template, policy, context)
            .map_err(GovernedAnalysisError::System)?;
        if system.is_replicated() && !kernel.descriptor().allows_replicated_systems() {
            return Err(GovernedAnalysisError::ReplicatedSystemUnsupported(
                kernel.descriptor().name().into(),
            ));
        }
        let mut assumptions: Vec<Assumption> = system.assumption().into_iter().collect();
        if matches!(policy.hydrogens, HydrogenPolicy::Exclude) {
            assumptions.push(Assumption::new(
                PolicyField::Hydrogens,
                "hydrogens excluded from the analysed atoms",
                AssumptionSource::Explicit,
                Impact::Unmeasured,
            ));
        }
        let resolution = system.structure.resolve_altlocs(policy);
        let (selected, selected_atoms, motions, input_quality, input_indeterminacy) =
            match resolution.value() {
                None => (
                    template.clone(),
                    Vec::new(),
                    Vec::new(),
                    Quality::Complete,
                    resolution.indeterminacy().cloned(),
                ),
                Some(chosen) => {
                    let chosen = match policy.hydrogens {
                        HydrogenPolicy::Exclude => without_hydrogens(&system.structure, chosen),
                        _ => chosen.clone(),
                    };
                    let kept: Vec<usize> = chosen.iter().map(|atom| atom as usize).collect();
                    let selected_atoms: Vec<usize> = kept
                        .iter()
                        .map(|&atom| system.placements[atom].source)
                        .collect();
                    let motions: Vec<Option<Rigid>> = kept
                        .iter()
                        .map(|&atom| system.placements[atom].motion)
                        .collect();
                    let structure = if chosen.len() == u64::from(system.structure.atom_count()) {
                        system.structure.clone()
                    } else {
                        system
                            .structure
                            .materialize(&chosen)
                            .map_err(GovernedAnalysisError::InvalidStructure)?
                    };
                    let resolved = match resolution.quality() {
                        Some(quality) => quality,
                        None => Quality::Complete,
                    };
                    let quality = if resolution.coverage.ambiguous > 0 {
                        resolved.worst(Quality::Ambiguous)
                    } else {
                        resolved
                    };
                    (structure, selected_atoms, motions, quality, None)
                }
            };
        assumptions.extend(resolution.assumptions);
        let input_indeterminacy = input_indeterminacy.or_else(|| {
            kernel
                .descriptor()
                .required_information()
                .iter()
                .find(|requirement| !requirement.satisfied_by(&selected))
                .map(|requirement| Indeterminacy::Other(requirement.reason().into()))
        });
        Ok(Self {
            template: selected,
            source_atom_count: template.atom_count() as usize,
            selected_atoms,
            motions,
            policy: policy.clone(),
            kernel,
            input_quality,
            input_indeterminacy,
            input_warnings: resolution.warnings,
            input_assumptions: assumptions,
        })
    }
}

impl<K: StructureKernel> FrameAnalysis for GovernedStructureAnalysis<'_, K> {
    type Output = Vec<K::Output>;
    type Partial = Vec<(usize, FrameRecord<K::Output>)>;
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
        if self.input_indeterminacy.is_some() {
            // No atoms could be chosen under the policy: nothing is run, and the
            // conclusion says why.
            return Ok(());
        }
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
        for ((target, source), motion) in positions
            .iter_mut()
            .zip(&self.selected_atoms)
            .zip(&self.motions)
        {
            let Some(position) = timestep.positions.get(*source) else {
                return Err(TrajectoryError::SelectionOutOfRange {
                    index: *source,
                    atoms: timestep.positions.len(),
                }
                .into());
            };
            *target = match motion {
                Some(motion) => motion.apply(*position),
                None => *position,
            };
        }
        let frame = editor
            .commit()
            .map_err(GovernedAnalysisError::InvalidStructure)?;
        let result = self
            .kernel
            .analyse_mapped(&frame, &self.policy, &self.selected_atoms, context)
            .map_err(GovernedAnalysisError::Kernel)?;
        let record = apply_missing_policy(result, self.policy.missing_atoms)?;
        partial.push((timestep.frame, record));
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
    let analysis = GovernedStructureAnalysis::new_in(topology, policy, kernel, context)?;
    run_analysis(trajectory, &analysis, context)
}

/// One frame's answer after the missing-data policy has had its say.
#[derive(Debug)]
pub struct FrameRecord<T> {
    outcome: Outcome<T>,
    quality: Quality,
    coverage: Coverage,
    warnings: Vec<molframe_core::diagnostic::Diagnostic>,
    assumptions: Vec<molframe_core::contract::Assumption>,
}

fn apply_missing_policy<T, E>(
    result: FrameKernelResult<T>,
    policy: MissingPolicy,
) -> Result<FrameRecord<T>, GovernedAnalysisError<E>> {
    validate_coverage(result.coverage)?;
    let (outcome, quality) = resolve_missing(result.value, result.quality, result.coverage, policy)
        .map_err(|error| match error {
            MissingPolicyError::Fail { missing, ambiguous } => {
                GovernedAnalysisError::MissingData { missing, ambiguous }
            }
            _ => GovernedAnalysisError::UnsupportedPolicyValue("missing_atoms"),
        })?;
    Ok(FrameRecord {
        outcome,
        quality,
        coverage: result.coverage,
        warnings: result.warnings,
        assumptions: result.assumptions,
    })
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

fn executor_reads(declared: &[PolicyField]) -> Vec<PolicyField> {
    let mut reads = EXECUTOR_READS.to_vec();
    reads.extend_from_slice(declared);
    reads
}

fn combine_frames<K: StructureKernel>(
    analysis: &GovernedStructureAnalysis<'_, K>,
    partial: Vec<(usize, FrameRecord<K::Output>)>,
) -> Result<Analysis<Vec<K::Output>>, GovernedAnalysisError<K::Error>> {
    let frame_count =
        i64::try_from(partial.len()).map_err(|_| GovernedAnalysisError::CoverageOverflow)?;
    let provenance = analysis
        .kernel
        .descriptor()
        .apply(Provenance::new(&analysis.policy).with_source(SourceRef::Memory))
        .with_parameter("frame_count", ParameterValue::Integer(frame_count))
        .with_policy_reads(&executor_reads(analysis.kernel.descriptor().policy_reads()));
    let mut quality = analysis.input_quality;
    let mut coverage = Coverage::default();
    let mut warnings = analysis.input_warnings.clone();
    let mut assumptions = analysis.input_assumptions.clone();
    let mut values = Vec::with_capacity(partial.len());
    let mut refusal = analysis.input_indeterminacy.clone();
    for (frame, record) in partial {
        quality = quality.worst(record.quality);
        coverage = add_coverage(coverage, record.coverage)?;
        warnings.extend(record.warnings);
        assumptions.extend(record.assumptions);
        match record.outcome {
            Outcome::Determinate(value) => values.push(value),
            Outcome::Indeterminate(reason) => {
                // One frame with no defensible answer leaves the series without one.
                refusal.get_or_insert(Indeterminacy::Frame {
                    frame,
                    reason: Box::new(reason),
                });
            }
        }
    }
    let outcome = match refusal {
        Some(reason) => Outcome::Indeterminate(reason),
        None => Outcome::Determinate(values),
    };
    let origin = analysis
        .selected_atoms
        .iter()
        .map(|&atom| u32::try_from(atom).map_err(|_| GovernedAnalysisError::CoverageOverflow))
        .collect::<Result<Vec<u32>, _>>()?;
    Ok(Analysis::from_parts(
        outcome,
        quality,
        coverage,
        warnings,
        assumptions,
        provenance,
    )
    .with_atom_origin(origin))
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

#[cfg(test)]
#[path = "execute_tests.rs"]
mod tests;
