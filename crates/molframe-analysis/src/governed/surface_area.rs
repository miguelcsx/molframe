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
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{ExecutionContext, Structure};
use molframe_surface::shrake_rupley;
use std::collections::BTreeMap;

/// Governed per-atom solvent-accessible surface area, in square ångström.
///
/// Shrake–Rupley with `points` samples per atom and a solvent sphere of radius
/// `probe`. The radii are the policy's `vdw_radii`. An atom whose element has no radius
/// in the set has no area (`NaN`) and does not occlude its neighbours; it is reported as
/// missing coverage, so the policy's `missing_atoms` decides whether that is a partial
/// answer, no answer, or a failure, and a warning names the elements and the policy field
/// (`vdw_radii`) so another radius set can be chosen deliberately.
///
/// Every atom of the analysed system occludes and gets an area: crystallographic waters
/// and heteroatoms included, which buries the protein atoms they touch. The result
/// records how many of each were counted; the protein-only convention of DSSP and
/// FreeSASA needs an input without them.
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
            let scope = scope_assumption(structure);
            let known: Vec<usize> = (0..radii.len())
                .filter(|&atom| radii[atom].is_finite())
                .collect();
            let missing = radii.len() - known.len();
            if missing == 0 {
                let areas = shrake_rupley(structure.positions(), &radii, probe, points, context)?;
                return Ok(add_scope(complete(structure, areas), scope));
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
            let result = FrameKernelResult::governed(areas, Quality::Partial, coverage)
                .with_assumption(Assumption::new(
                    PolicyField::VdwRadii,
                    format!(
                        "{absent} atoms have an element with no radius in the {} set; they have no area and do not occlude",
                        set.name()
                    ),
                    AssumptionSource::Inferred,
                    Impact::Unmeasured,
                ))
                .with_warning(
                    Diagnostic::new(Code::W6004)
                        .about_field(PolicyField::VdwRadii.name())
                        .with_message(format!(
                            "the {} radius set has no radius for {}; those atoms are excluded, so the total is partial (choose a set that covers them, such as alvarez, or remove them from the selection)",
                            set.name(),
                            missing_elements(structure, &radii),
                        )),
                );
            Ok(add_scope(result, scope))
        },
    )
}

fn add_scope<T>(result: FrameKernelResult<T>, scope: Option<Assumption>) -> FrameKernelResult<T> {
    match scope {
        Some(scope) => result.with_assumption(scope),
        None => result,
    }
}

/// The elements without a radius and how many atoms each has, in symbol order.
fn missing_elements(structure: &Structure, radii: &[f32]) -> String {
    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for atom in structure.data().atoms() {
        let missing = radii
            .get(atom.index().as_usize())
            .is_some_and(|radius| !radius.is_finite());
        if missing {
            let symbol = atom
                .element()
                .map_or("unknown element", |element| element.symbol());
            *counts.entry(symbol).or_default() += 1;
        }
    }
    counts
        .iter()
        .map(|(symbol, count)| format!("{symbol} ({count})"))
        .collect::<Vec<_>>()
        .join(", ")
}

/// A note when the analysed atoms include waters or other heteroatoms, which
/// occlude and receive areas like any polymer atom.
fn scope_assumption(structure: &Structure) -> Option<Assumption> {
    let (mut water, mut hetero) = (0_u32, 0_u32);
    for atom in structure.data().atoms() {
        let Some(residue) = atom.residue() else {
            continue;
        };
        if matches!(atom.component_name(), Some("HOH" | "WAT" | "DOD" | "H2O")) {
            water += 1;
        } else if residue.is_het() {
            hetero += 1;
        }
    }
    (water + hetero > 0).then(|| {
        Assumption::new(
            PolicyField::ContactDef,
            format!(
                "{water} water and {hetero} other heteroatom atoms are included as occluders and receive areas; analyse an input without them for a protein-only area"
            ),
            AssumptionSource::ProfileDefault,
            Impact::Unmeasured,
        )
    })
}
