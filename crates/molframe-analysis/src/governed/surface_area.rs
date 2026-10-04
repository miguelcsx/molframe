//! Solvent-accessible surface area under the policy's radii.
//!
//! The area of an atom is a function of the radii it and its neighbours are given, of
//! whether the hydrogens the file carries occlude, and of which system it sits in (the
//! unit alone, its assembly, or its crystal). Each is a decision, so each is applied
//! here: the radii are read from `vdw_radii`, and the executor resolves hydrogens and
//! the system before the kernel sees the atoms.

use super::common::{complete, descriptor, float};
use super::definition::DefinitionError;
use super::{StructureKernel, structure_kernel};
use molframe_chem::{RadiusSet, atom_radii};
use molframe_core::contract::{AnalysisPolicy, PolicyField};
use molframe_core::{ExecutionContext, Structure};
use molframe_surface::shrake_rupley;

/// Governed per-atom solvent-accessible surface area, in square ångström.
///
/// Shrake–Rupley with `points` samples per atom and a solvent sphere of radius
/// `probe`. The radii are the policy's `vdw_radii`; an atom whose element has no radius in
/// the set has no area, and the analysis refuses rather than guess one.
#[must_use]
pub fn sasa_kernel(
    probe: f32,
    points: u16,
) -> impl StructureKernel<Output = Vec<f64>, Error = DefinitionError> {
    structure_kernel(
        descriptor("solvent-accessible-surface")
            .estimating("the area of each atom's van der Waals surface a solvent sphere can touch")
            .reading(&[PolicyField::VdwRadii])
            .with_parameter("probe", float(probe))
            .with_parameter("points", float(f32::from(points))),
        move |structure: &Structure, policy: &AnalysisPolicy, context: &ExecutionContext| {
            let set: RadiusSet =
                policy.vdw_radii.name().parse().map_err(|_| {
                    DefinitionError::UnknownRadiusSet(policy.vdw_radii.name().into())
                })?;
            let radii = atom_radii(structure, set);
            let unknown = radii.iter().filter(|radius| !radius.is_finite()).count();
            if unknown > 0 {
                return Err(DefinitionError::UnknownRadius {
                    atoms: unknown,
                    set: set.name(),
                });
            }
            let areas = shrake_rupley(structure.positions(), &radii, probe, points, context)?;
            Ok(complete(structure, areas))
        },
    )
}
