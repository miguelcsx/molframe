//! Governed native-contact comparison with explicit atom correspondence.

use super::common::{backend, complete, descriptor, float, integer};
use super::{StructureKernel, mapped_structure_kernel};
use crate::{
    NativeContacts, NativeError, match_atoms_by_name, native_contact_fraction,
    native_contact_fraction_matched,
};
use molframe_core::contract::{
    AnalysisPolicy, ModelChoice, ParameterValue, PeriodicPolicy, PolicyField,
};
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
    #[error(
        "source atom index {index} is absent from target with {atoms} atoms; \
         structures that number atoms differently need a correspondence by name or an explicit one"
    )]
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
    /// An explicit correspondence names a reference atom the reference lacks.
    #[error("correspondence names reference atom {0}, which is not in the reference")]
    ReferenceAtomAbsent(u32),
    /// An explicit correspondence names a target atom the target lacks.
    #[error("correspondence names target atom {index}, but the target has {atoms} atoms")]
    CorrespondenceTargetAbsent {
        /// The offending target atom index.
        index: u32,
        /// Available target atoms.
        atoms: usize,
    },
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

/// How the target's atoms are paired with the reference's.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NativeCorrespondence {
    /// Atom `i` of the reference is atom `i` of the target.
    Index,
    /// Atoms are paired by chain, residue number, insertion code, residue name
    /// and atom name.
    Names,
    /// Explicit `(reference atom, target atom)` pairs, in the original atom
    /// numbering of each structure. Reference atoms not listed take no part.
    Explicit(Vec<(u32, u32)>),
}

impl NativeCorrespondence {
    fn label(&self) -> &'static str {
        match self {
            Self::Index => "index",
            Self::Names => "names",
            Self::Explicit(_) => "explicit",
        }
    }
}

/// The target atom standing for each atom of the resolved reference.
fn explicit_map(
    pairs: &[(u32, u32)],
    source_atoms: &[usize],
    target_atoms: usize,
) -> Result<Vec<Option<u32>>, GovernedNativeError> {
    let mut local_of = std::collections::HashMap::with_capacity(source_atoms.len());
    for (local, source) in source_atoms.iter().enumerate() {
        local_of.insert(*source, local);
    }
    let mut map = vec![None; source_atoms.len()];
    for &(reference, target) in pairs {
        if target as usize >= target_atoms {
            return Err(GovernedNativeError::CorrespondenceTargetAbsent {
                index: target,
                atoms: target_atoms,
            });
        }
        let slot = local_of
            .get(&(reference as usize))
            .and_then(|local| map.get_mut(*local))
            .ok_or(GovernedNativeError::ReferenceAtomAbsent(reference))?;
        *slot = Some(target);
    }
    Ok(map)
}

fn matched(
    reference: &Structure,
    target: &Structure,
    source_atoms: &[usize],
    correspondence: &NativeCorrespondence,
    policy: &AnalysisPolicy,
) -> Result<(Vec<[f32; 3]>, Vec<Option<u32>>), GovernedNativeError> {
    let model = target_model(policy)?;
    let positions = target
        .model_positions(ModelIndex::new(model))
        .ok_or(GovernedNativeError::MissingTargetModel(model))?;
    let map = match correspondence {
        NativeCorrespondence::Names => match_atoms_by_name(reference, target),
        NativeCorrespondence::Explicit(pairs) => {
            explicit_map(pairs, source_atoms, positions.len())?
        }
        NativeCorrespondence::Index => Vec::new(),
    };
    Ok((positions.to_vec(), map))
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

/// Governed native-contact kernel pairing atoms by index.
///
/// `tolerance` is an additive distance in ångström: a reference contact within
/// `cutoff` is kept when the target pair is within `cutoff + tolerance`.
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
    native_contact_fraction_kernel_with(
        target,
        cutoff,
        tolerance,
        NativeCorrespondence::Index,
        spatial,
    )
}

/// Governed native-contact kernel with an explicit atom correspondence.
///
/// Reference atoms without a counterpart in the target take no part in any
/// contact. The result reports how many reference atoms were matched.
#[must_use]
pub fn native_contact_fraction_kernel_with(
    target: &Structure,
    cutoff: f32,
    tolerance: f32,
    correspondence: NativeCorrespondence,
    spatial: SpatialBackend,
) -> impl StructureKernel<Output = NativeContacts, Error = GovernedNativeError> + '_ {
    let pairs = match &correspondence {
        NativeCorrespondence::Explicit(pairs) => pairs.len(),
        _ => 0,
    };
    mapped_structure_kernel(
        descriptor("native-contact-fraction")
            .estimating("the fraction of the reference structure's contacts that the frame keeps")
            .reading(&[PolicyField::Periodic])
            .without_replicated_systems()
            .with_parameter("cutoff", float(cutoff))
            .with_parameter("tolerance", float(tolerance))
            .with_parameter(
                "correspondence",
                ParameterValue::Text(correspondence.label().into()),
            )
            .with_parameter("correspondence_pairs", integer(pairs))
            .with_parameter("spatial_backend", backend(spatial)),
        move |reference: &Structure,
              policy: &AnalysisPolicy,
              source_atoms: &[usize],
              context: &ExecutionContext| {
            if !matches!(policy.periodic, PeriodicPolicy::None) {
                return Err(GovernedNativeError::PeriodicUnsupported);
            }
            let value = if matches!(correspondence, NativeCorrespondence::Index) {
                let target = prepare_target(target, source_atoms, policy, context)?;
                native_contact_fraction(reference, &target, cutoff, tolerance, spatial, context)?
            } else {
                let (positions, map) =
                    matched(reference, target, source_atoms, &correspondence, policy)?;
                native_contact_fraction_matched(
                    reference, &positions, &map, cutoff, tolerance, spatial, context,
                )?
            };
            Ok(complete(reference, value))
        },
    )
}
