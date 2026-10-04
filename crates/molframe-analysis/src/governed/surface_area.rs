//! Solvent-accessible surface area under the policy's radii.
//!
//! The area of an atom is a function of the radii it and its neighbours are given, of
//! whether the hydrogens the file carries occlude, and of which system it sits in (the
//! unit alone, its assembly, or its crystal). Each is a decision, so each is applied
//! here: the radii are read from `vdw_radii`, and the executor resolves hydrogens and
//! the system before the kernel sees the atoms.

use super::common::{complete, descriptor, float};
use super::definition::DefinitionError;
use super::{FrameKernelResult, StructureKernel, structure_kernel};
use molframe_chem::{RadiusSet, atom_radii};
use molframe_core::contract::{
    AnalysisPolicy, Assumption, AssumptionSource, Coverage, Impact, PolicyField, Quality,
};
use molframe_core::{ExecutionContext, Structure};
use molframe_surface::shrake_rupley;

/// Governed per-atom solvent-accessible surface area, in square ångström.
///
/// Shrake–Rupley with `points` samples per atom and a solvent sphere of radius
/// `probe`. The radii are the policy's `vdw_radii`. An atom whose element has no radius
/// in the set has no area (`NaN`) and does not occlude its neighbours; it is reported as
/// missing coverage, so the policy's `missing_atoms` decides whether that is a partial
/// answer, no answer, or a failure.
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
            let known: Vec<usize> = (0..radii.len())
                .filter(|&atom| radii[atom].is_finite())
                .collect();
            let missing = radii.len() - known.len();
            if missing == 0 {
                let areas = shrake_rupley(structure.positions(), &radii, probe, points, context)?;
                return Ok(complete(structure, areas));
            }
            let positions: Vec<[f32; 3]> = known
                .iter()
                .map(|&atom| structure.positions()[atom])
                .collect();
            let kept: Vec<f32> = known.iter().map(|&atom| radii[atom]).collect();
            let measured = shrake_rupley(&positions, &kept, probe, points, context)?;
            let mut areas = vec![f64::NAN; radii.len()];
            for (&atom, area) in known.iter().zip(measured) {
                areas[atom] = area;
            }
            let atoms = structure.atom_count();
            let absent = u32::try_from(missing).map_or(atoms, |count| count.min(atoms));
            let coverage = Coverage {
                intended: atoms,
                used: atoms - absent,
                missing: absent,
                ambiguous: 0,
            };
            Ok(FrameKernelResult::governed(areas, Quality::Partial, coverage).with_assumption(
                Assumption::new(
                    PolicyField::VdwRadii,
                    format!(
                        "{absent} atoms have an element with no radius in the {} set; they have no area and do not occlude",
                        set.name()
                    ),
                    AssumptionSource::Inferred,
                    Impact::Unmeasured,
                ),
            ))
        },
    )
}
