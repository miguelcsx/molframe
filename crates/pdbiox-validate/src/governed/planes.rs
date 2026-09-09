//! Governed explicit plane restraints with atom source-map projection.

use crate::{
    PlanarityError, PlanarityOptions, PlaneRestraint, PlaneRestraintReport,
    plane_restraint_outliers,
};
use pdbiox_analysis::{
    AnalysisDescriptor, FrameKernelResult, StructureKernel, mapped_structure_kernel,
};
use pdbiox_core::contract::{AnalysisPolicy, Coverage, ParameterValue, Status};
use pdbiox_core::{AtomSelection, ExecutionContext, Structure};

/// Plane fitting or coverage failure.
#[derive(Debug, thiserror::Error)]
pub enum PlaneRestraintKernelError {
    /// Plane fitting failure.
    #[error(transparent)]
    Validation(#[from] PlanarityError),
    /// Coverage count overflow.
    #[error("plane restraint count exceeds u32 coverage limits")]
    CoverageOverflow,
}

/// Governed non-aromatic explicit plane restraints.
#[must_use]
pub fn plane_restraint_outliers_kernel(
    restraints: &[PlaneRestraint],
    options: PlanarityOptions,
) -> impl StructureKernel<Output = PlaneRestraintReport, Error = PlaneRestraintKernelError> + '_ {
    mapped_structure_kernel(
        AnalysisDescriptor::new("explicit-plane-restraints", "1").with_parameter(
            "maximum_deviation",
            ParameterValue::Float(options.maximum_deviation.to_bits()),
        ),
        move |structure: &Structure,
              _policy: &AnalysisPolicy,
              source_atoms: &[usize],
              _context: &ExecutionContext| {
            let projected: Vec<_> = restraints
                .iter()
                .map(|restraint| project(restraint, source_atoms))
                .collect();
            let value = plane_restraint_outliers(structure, &projected, options)?;
            let intended = u32::try_from(value.intended)
                .map_err(|_| PlaneRestraintKernelError::CoverageOverflow)?;
            let used = u32::try_from(value.assessed)
                .map_err(|_| PlaneRestraintKernelError::CoverageOverflow)?;
            Ok(FrameKernelResult::governed(
                value,
                if intended == used {
                    Status::Complete
                } else {
                    Status::Partial
                },
                Coverage {
                    intended,
                    used,
                    missing: intended - used,
                    ambiguous: 0,
                },
            ))
        },
    )
}

fn project(restraint: &PlaneRestraint, source_atoms: &[usize]) -> PlaneRestraint {
    let atoms = source_atoms
        .iter()
        .enumerate()
        .filter_map(|(target, source)| {
            let source = u32::try_from(*source).ok()?;
            restraint
                .atoms
                .contains(source)
                .then(|| u32::try_from(target).ok())?
        })
        .collect();
    PlaneRestraint {
        id: restraint.id.clone(),
        atoms: AtomSelection::from_sorted(atoms),
    }
}
