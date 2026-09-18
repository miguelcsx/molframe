//! Governed native-contact comparison with explicit atom correspondence.

use super::common::{backend, complete, descriptor, float};
use super::{StructureKernel, mapped_structure_kernel};
use crate::{NativeContacts, NativeError, native_contact_fraction};
use molframe_core::contract::{AnalysisPolicy, ModelChoice, PeriodicPolicy};
use molframe_core::index::ModelIndex;
use molframe_core::{AtomSelection, Diagnostic, ExecutionContext, Structure};
use molframe_spatial::SpatialBackend;

/// Setup or scientific failure for a governed native-contact comparison.
#[derive(Debug, thiserror::Error)]
pub enum GovernedNativeError {
    /// Coordinate copy-on-write could not fit the shared execution account.
    #[error("target coordinate copy exceeds the execution memory budget: {0}")]
    Memory(#[from] molframe_core::MemoryBudgetError),
    /// The native-contact kernel refused the prepared structures.
    #[error(transparent)]
    Native(#[from] NativeError),
    /// Target topology materialisation failed.
    #[error("target atom correspondence could not be materialised")]
    InvalidTarget(Vec<Diagnostic>),
    /// A source atom index cannot be represented by the public topology.
    #[error("source atom index {0} exceeds the public u32 topology range")]
    SourceIndexOverflow(usize),
    /// A corresponding atom is absent from the target topology.
    #[error("source atom index {index} is absent from target with {atoms} atoms")]
    TargetAtomOutOfRange {
        /// Requested source atom position.
        index: usize,
        /// Available target atoms.
        atoms: u32,
    },
    /// The requested model is absent from the target structure.
    #[error("target does not contain requested model {0}")]
    MissingTargetModel(u32),
    /// Multi-model policy requires the trajectory API rather than one target.
    #[error("native-contact comparison requires one target model")]
    MultipleModelsRequested,
    /// Native-contact distance evaluation does not support periodic images.
    #[error("native-contact comparison currently requires non-periodic coordinates")]
    PeriodicUnsupported,
    /// The installed policy has a variant this crate does not understand.
    #[error("unsupported analysis policy field: {0}")]
    UnsupportedPolicy(&'static str),
}

fn target_model(policy: &AnalysisPolicy) -> Result<u32, GovernedNativeError> {
    match policy.model {
        ModelChoice::First => Ok(0),
        ModelChoice::Index(index) => Ok(index),
        ModelChoice::All | ModelChoice::Ensemble => {
            Err(GovernedNativeError::MultipleModelsRequested)
        }
        _ => Err(GovernedNativeError::UnsupportedPolicy("model")),
    }
}

fn prepare_target(
    target: &Structure,
    source_atoms: &[usize],
    policy: &AnalysisPolicy,
    context: &ExecutionContext,
) -> Result<Structure, GovernedNativeError> {
    let model = target_model(policy)?;
    let positions = target
        .model_positions(ModelIndex::new(model))
        .ok_or(GovernedNativeError::MissingTargetModel(model))?;
    let selection = source_atoms
        .iter()
        .map(|atom| {
            u32::try_from(*atom).map_err(|_| GovernedNativeError::SourceIndexOverflow(*atom))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let selected_positions = source_atoms
        .iter()
        .map(|atom| {
            positions
                .get(*atom)
                .copied()
                .ok_or(GovernedNativeError::TargetAtomOutOfRange {
                    index: *atom,
                    atoms: target.atom_count(),
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let selected = target
        .materialize(&AtomSelection::from_sorted(selection))
        .map_err(GovernedNativeError::InvalidTarget)?;
    let mut editor = selected.edit_coordinates(context)?;
    let output = editor
        .positions_mut(ModelIndex::new(0))
        .ok_or(GovernedNativeError::MissingTargetModel(0))?;
    output.copy_from_slice(&selected_positions);
    editor.commit().map_err(GovernedNativeError::InvalidTarget)
}

/// Governed native-contact kernel sharing reference atom correspondence with a target.
///
/// Alternate conformations are resolved once on the reference by the common
/// executor. The resulting source atom map is then projected onto the target,
/// which preserves the explicit same-index contract of native-contact Q.
#[must_use]
pub fn native_contact_fraction_kernel(
    target: &Structure,
    cutoff: f32,
    tolerance: f32,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = NativeContacts, Error = GovernedNativeError> + '_ {
    mapped_structure_kernel(
        descriptor("native-contact-fraction")
            .with_parameter("cutoff", float(cutoff))
            .with_parameter("tolerance", float(tolerance))
            .with_parameter("spatial_backend", backend(spatial)),
        move |reference: &Structure,
              policy: &AnalysisPolicy,
              source_atoms: &[usize],
              context: &ExecutionContext| {
            if !matches!(policy.periodic, PeriodicPolicy::None) {
                return Err(GovernedNativeError::PeriodicUnsupported);
            }
            let target = prepare_target(target, source_atoms, policy, context)?;
            let value =
                native_contact_fraction(reference, &target, cutoff, tolerance, spatial, context)?;
            Ok(complete(reference, value))
        },
    )
}
