//! Governed B-factor distribution and TLS consistency.

use crate::{
    BFactorDistribution, BFactorError, TlsBFactorReport, TlsGroup, b_factor_distribution,
    tls_b_factor_consistency,
};
use pdbiox_analysis::{
    AnalysisDescriptor, FrameKernelResult, StructureKernel, mapped_structure_kernel,
};
use pdbiox_core::contract::{AnalysisPolicy, Coverage, ParameterValue, Status};
use pdbiox_core::{AtomSelection, ExecutionContext, Structure};

/// B-factor kernel or coverage failure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, thiserror::Error)]
pub enum BFactorKernelError {
    /// Scientific input failure.
    #[error(transparent)]
    Validation(#[from] BFactorError),
    /// Coverage count overflow.
    #[error("B-factor count exceeds u32 coverage limits")]
    CoverageOverflow,
}

/// Governed B-factor distribution over a source-topology selection.
#[must_use]
pub fn b_factor_distribution_kernel(
    selection: &AtomSelection,
    outlier_standard_deviations: f64,
) -> impl StructureKernel<Output = BFactorDistribution, Error = BFactorKernelError> + '_ {
    mapped_structure_kernel(
        AnalysisDescriptor::new("b-factor-distribution", "1").with_parameter(
            "outlier_standard_deviations",
            ParameterValue::Float(outlier_standard_deviations.to_bits()),
        ),
        move |structure: &Structure,
              _policy: &AnalysisPolicy,
              source_atoms: &[usize],
              _context: &ExecutionContext| {
            let value = b_factor_distribution(
                structure,
                &project(selection, source_atoms),
                outlier_standard_deviations,
            )?;
            report(value, |value| (value.intended, value.assessed))
        },
    )
}

/// Governed consistency against explicit TLS declarations.
#[must_use]
pub fn tls_b_factor_consistency_kernel(
    groups: &[TlsGroup],
    maximum_absolute_deviation: f64,
    symmetry_tolerance: f64,
) -> impl StructureKernel<Output = TlsBFactorReport, Error = BFactorKernelError> + '_ {
    mapped_structure_kernel(
        AnalysisDescriptor::new("tls-b-factor-consistency", "1")
            .with_parameter(
                "maximum_absolute_deviation",
                ParameterValue::Float(maximum_absolute_deviation.to_bits()),
            )
            .with_parameter(
                "symmetry_tolerance",
                ParameterValue::Float(symmetry_tolerance.to_bits()),
            )
            .with_parameter("groups", ParameterValue::Text(format!("{groups:?}").into())),
        move |structure: &Structure,
              _policy: &AnalysisPolicy,
              source_atoms: &[usize],
              _context: &ExecutionContext| {
            let projected: Vec<_> = groups
                .iter()
                .map(|group| TlsGroup {
                    id: group.id.clone(),
                    atoms: project(&group.atoms, source_atoms),
                    model: group.model,
                })
                .collect();
            let value = tls_b_factor_consistency(
                structure,
                &projected,
                maximum_absolute_deviation,
                symmetry_tolerance,
            )?;
            report(value, |value| (value.intended, value.assessed))
        },
    )
}

fn project(selection: &AtomSelection, source_atoms: &[usize]) -> AtomSelection {
    AtomSelection::from_sorted(
        source_atoms
            .iter()
            .enumerate()
            .filter_map(|(target, source)| {
                let source = u32::try_from(*source).ok()?;
                selection
                    .contains(source)
                    .then(|| u32::try_from(target).ok())?
            })
            .collect(),
    )
}

fn report<T>(
    value: T,
    counts: impl FnOnce(&T) -> (usize, usize),
) -> Result<FrameKernelResult<T>, BFactorKernelError> {
    let (intended, used) = counts(&value);
    let intended = u32::try_from(intended).map_err(|_| BFactorKernelError::CoverageOverflow)?;
    let used = u32::try_from(used).map_err(|_| BFactorKernelError::CoverageOverflow)?;
    let missing = intended
        .checked_sub(used)
        .ok_or(BFactorKernelError::CoverageOverflow)?;
    Ok(FrameKernelResult::governed(
        value,
        if missing == 0 {
            Status::Complete
        } else {
            Status::Partial
        },
        Coverage {
            intended,
            used,
            missing,
            ambiguous: 0,
        },
    ))
}
