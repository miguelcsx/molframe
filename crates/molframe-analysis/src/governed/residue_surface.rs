//! Per-residue solvent-accessible surface area under the policy's radii.
//!
//! The same decisions as [`super::sasa_kernel`] (one option set, the radii rule of
//! [`RadiusSelection`], the applied table recorded as an assumption), aggregated by
//! residue by [`crate::residue_areas`]: each residue's area is the sum of its atoms'
//! areas from the one surface calculation.

use super::FrameKernelResult;
use super::common::{descriptor, float};
use super::surface_area::radii_assumption;
use super::{StructureKernel, structure_kernel};
use crate::residue_sasa::{ResidueAreas, ResidueSasaError, residue_areas};
use molframe_chem::{RadiusSelection, TableError};
use molframe_core::contract::{AnalysisPolicy, Coverage, ParameterValue, PolicyField, Quality};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{ExecutionContext, Structure, diagnostic_from};
use molframe_surface::{SasaMethod, SasaOptions};

/// Why governed per-residue areas could not be computed.
#[derive(Debug, thiserror::Error)]
pub enum ResidueKernelError {
    /// The radius table conflicts with the policy or could not be supplied.
    #[error(transparent)]
    Radii(#[from] TableError),
    /// The per-residue calculation was refused.
    #[error(transparent)]
    Residues(#[from] ResidueSasaError),
}

diagnostic_from!(ResidueKernelError, |error| match error {
    ResidueKernelError::Radii(inner) => Diagnostic::from(inner).code(),
    ResidueKernelError::Residues(ResidueSasaError::NoAtoms) => Code::E5103,
    ResidueKernelError::Residues(ResidueSasaError::Surface(inner)) =>
        Diagnostic::from(inner).code(),
});

/// Governed per-residue solvent-accessible area, in square ångström.
///
/// `options` and `selection` mean what they do for [`super::sasa_kernel`]. Included
/// atoms whose element has no radius in the table are reported as missing coverage, so
/// the policy's `missing_atoms` decides whether the answer is partial, and a warning
/// names the policy field (`vdw_radii`).
#[must_use]
pub fn residue_sasa_kernel(
    options: SasaOptions,
    selection: RadiusSelection,
) -> impl StructureKernel<Output = ResidueAreas, Error = ResidueKernelError> {
    let method = match options.method {
        SasaMethod::ShrakeRupley { points } => format!("shrake-rupley:{points}"),
        SasaMethod::LeeRichards { slices } => format!("lee-richards:{slices}"),
    };
    structure_kernel(
        descriptor("residue-solvent-accessible-surface")
            .estimating(
                "the area of each residue's van der Waals surface a solvent sphere can touch",
            )
            .reading(&[PolicyField::VdwRadii])
            .with_parameter("probe", float(options.probe))
            .with_parameter("method", ParameterValue::Text(method.into()))
            .with_parameter(
                "hetero",
                ParameterValue::Text(options.hetero.to_string().into()),
            )
            .with_parameter(
                "metals",
                ParameterValue::Text(options.metals.to_string().into()),
            )
            .with_parameter(
                "radii_table",
                ParameterValue::Text(options.radii.identity().to_string().into()),
            ),
        move |structure: &Structure, policy: &AnalysisPolicy, context: &ExecutionContext| {
            let applied = radii_assumption(&options, &selection, policy)?;
            let areas = residue_areas(structure, &options, context)?;
            let used: usize = areas.residues.iter().map(|residue| residue.atoms).sum();
            let missing = areas.missing_radius;
            let table = options.radii.identity().name.to_string();
            let intended = count(used + missing);
            if missing == 0 {
                let used = count(used);
                return Ok(FrameKernelResult::complete(areas, used).with_assumption(applied));
            }
            let coverage = Coverage {
                intended,
                used: count(used),
                missing: count(missing),
                ambiguous: 0,
            };
            let warning = Diagnostic::new(Code::W6004)
                .about_field(PolicyField::VdwRadii.name())
                .with_message(format!(
                    "the {table} radius set has no radius for {missing} included atoms; they are left out of their residues' areas and do not occlude (choose a set that covers them, such as alvarez)"
                ));
            Ok(
                FrameKernelResult::governed(areas, Quality::Partial, coverage)
                    .with_assumption(applied)
                    .with_warning(warning),
            )
        },
    )
}

fn count(value: usize) -> u32 {
    match u32::try_from(value) {
        Ok(count) => count,
        Err(_) => u32::MAX,
    }
}

#[cfg(test)]
#[path = "residue_surface_tests.rs"]
mod tests;
