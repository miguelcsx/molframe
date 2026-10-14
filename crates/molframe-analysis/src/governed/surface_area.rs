//! Solvent-accessible surface area under the policy's radii.
//!
//! The area of an atom is a function of the radii it and its neighbours are given, of
//! whether the hydrogens the file carries occlude, and of which system it sits in (the
//! unit alone, its assembly, or its crystal). Each is a decision, so each is applied
//! here: the radii are read from `vdw_radii`, and the executor resolves hydrogens and
//! the system before the kernel sees the atoms.

use super::common::{descriptor, float};
use super::definition::DefinitionError;
use super::{FrameKernelResult, StructureKernel, structure_kernel};
use molframe_chem::included_atoms;
use molframe_core::contract::{
    AnalysisPolicy, Assumption, AssumptionSource, Coverage, Impact, ParameterValue, PolicyField,
    Quality,
};
use molframe_core::diagnostic::{Code, Diagnostic};
use molframe_core::{ExecutionContext, Structure};
use molframe_surface::{SasaMethod, SasaOptions, measure_atoms};
use std::collections::BTreeMap;

/// Governed per-atom solvent-accessible surface area, in square ångström.
///
/// `options` is the one option set every front end shares: the sampling method and its
/// density, the probe radius, which heteroatom classes take part (`hetero`, the polymer
/// alone by default) and where a metal the radius set lacks takes its radius (`metals`).
/// The radii are `options.radii`, a table the caller resolved. The result records the
/// table that was applied (name, version, source, licence) as a `vdw_radii` assumption,
/// which is a profile default when it carries the name the policy's `vdw_radii` gives and
/// explicit otherwise. An atom the options
/// leave out, or whose element has no radius in the set, has no area (`NaN`) and does not
/// occlude its neighbours. A left-out heteroatom is a choice and not a gap; an included
/// atom with no radius is reported as missing coverage, so the policy's `missing_atoms`
/// decides whether that is a partial answer, no answer, or a failure, and a warning names
/// the elements and the policy field (`vdw_radii`) so another radius set can be chosen
/// deliberately.
#[must_use]
pub fn sasa_kernel(
    options: SasaOptions,
) -> impl StructureKernel<Output = Vec<f64>, Error = DefinitionError> {
    structure_kernel(
        descriptor("solvent-accessible-surface")
            .estimating("the area of each atom's van der Waals surface a solvent sphere can touch")
            .reading(&[PolicyField::VdwRadii])
            .with_parameter("probe", float(options.probe))
            .with_parameter("method", method_parameter(options.method))
            .with_parameter("hetero", text(&options.hetero.to_string()))
            .with_parameter("metals", text(&options.metals.to_string()))
            .with_parameter("radii_table", radii_parameter(&options)),
        move |structure: &Structure, policy: &AnalysisPolicy, context: &ExecutionContext| {
            let origin = if options.radii.identity().name == policy.vdw_radii.name() {
                AssumptionSource::ProfileDefault
            } else {
                AssumptionSource::Explicit
            };
            let applied = applied_assumption(&options, origin);
            let table = options.radii.identity().name.to_string();
            let radii = options.atom_radii(structure);
            let included = included_atoms(structure, options.hetero);
            let areas = measure_atoms(structure.positions(), &radii, &options, context)?;
            let intended = included.iter().filter(|&&atom| atom).count();
            let used = areas.iter().filter(|area| area.is_finite()).count();
            let scope = scope_assumption(structure, &areas);
            let missing = intended.saturating_sub(used);
            if missing == 0 {
                let result = FrameKernelResult::complete(areas, count(used));
                return Ok(add_scope(result.with_assumption(applied), scope));
            }
            let coverage = Coverage {
                intended: count(intended),
                used: count(used),
                missing: count(missing),
                ambiguous: 0,
            };
            let result = FrameKernelResult::governed(areas, Quality::Partial, coverage)
                .with_assumption(Assumption::new(
                    PolicyField::VdwRadii,
                    format!(
                        "{missing} atoms have an element with no radius in the {table} set; they have no area and do not occlude"
                    ),
                    AssumptionSource::Inferred,
                    Impact::Unmeasured,
                ))
                .with_warning(
                    Diagnostic::new(Code::W6004)
                        .about_field(PolicyField::VdwRadii.name())
                        .with_message(format!(
                            "the {} radius set has no radius for {}; those atoms are excluded, so the total is partial (choose a set that covers them, such as alvarez, or remove them from the selection)",
                            table,
                            missing_elements(structure, &radii, &included),
                        )),
                );
            Ok(add_scope(result.with_assumption(applied), scope))
        },
    )
}

fn count(value: usize) -> u32 {
    match u32::try_from(value) {
        Ok(count) => count,
        Err(_) => u32::MAX,
    }
}

fn text(value: &str) -> ParameterValue {
    ParameterValue::Text(value.into())
}

fn radii_parameter(options: &SasaOptions) -> ParameterValue {
    text(&options.radii.identity().to_string())
}

/// The radius policy as applied: every table behind the radii, by name, version, source
/// and licence.
fn applied_assumption(options: &SasaOptions, origin: AssumptionSource) -> Assumption {
    Assumption::new(
        PolicyField::VdwRadii,
        options.applied_policy().to_string(),
        origin,
        Impact::Unmeasured,
    )
}

fn method_parameter(method: SasaMethod) -> ParameterValue {
    match method {
        SasaMethod::ShrakeRupley { points } => text(&format!("shrake-rupley:{points}")),
        SasaMethod::LeeRichards { slices } => text(&format!("lee-richards:{slices}")),
    }
}

fn add_scope<T>(result: FrameKernelResult<T>, scope: Option<Assumption>) -> FrameKernelResult<T> {
    match scope {
        Some(scope) => result.with_assumption(scope),
        None => result,
    }
}

/// The elements without a radius and how many atoms each has, in symbol order.
fn missing_elements(structure: &Structure, radii: &[f32], included: &[bool]) -> String {
    let mut counts: BTreeMap<&'static str, usize> = BTreeMap::new();
    for atom in structure.data().atoms() {
        let index = atom.index().as_usize();
        let missing = matches!(included.get(index), Some(true))
            && radii
                .get(index)
                .is_some_and(|radius| !(radius.is_finite() && *radius > 0.0));
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

/// A note when waters or other heteroatoms were measured, so they occlude and receive
/// areas like any polymer atom.
fn scope_assumption(structure: &Structure, areas: &[f64]) -> Option<Assumption> {
    let (mut water, mut hetero) = (0_u32, 0_u32);
    for atom in structure.data().atoms() {
        let measured = areas
            .get(atom.index().as_usize())
            .is_some_and(|area| area.is_finite());
        let Some(residue) = atom.residue().filter(|_| measured) else {
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
                "{water} water and {hetero} other heteroatom atoms are measured as occluders and receive areas because the options include them"
            ),
            AssumptionSource::Explicit,
            Impact::Unmeasured,
        )
    })
}
