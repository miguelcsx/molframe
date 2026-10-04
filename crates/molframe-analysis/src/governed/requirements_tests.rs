use super::{ForbiddenResolution, GovernedAnalysisError, Requirement, analyse_structure};
use crate::policy_execution::{AnalysisDescriptor, FrameKernelResult, structure_kernel};
use molframe_core::ExecutionContext;
use molframe_core::contract::{AnalysisPolicy, HydrogenPolicy, PolicyField};
use molframe_core::io::{InputBuffer, ReadOptions};
use molframe_core::structure::Structure;
use std::convert::Infallible;

const HEADER: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n";

fn structure(with_hydrogen: bool) -> Structure {
    let mut source = format!("{HEADER}ATOM 1 N N1 LIG A 1 0 0 0\n");
    if with_hydrogen {
        source.push_str("ATOM 2 H H1 LIG A 1 1 0 0\n");
    }
    let input = InputBuffer::from_bytes(source.into_bytes());
    match molframe_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn atoms_kernel(
    descriptor: AnalysisDescriptor,
) -> impl crate::policy_execution::StructureKernel<Output = u32, Error = Infallible> {
    structure_kernel(
        descriptor,
        |structure: &Structure, _policy: &AnalysisPolicy, _context: &ExecutionContext| {
            Ok(FrameKernelResult::complete(
                structure.atom_count(),
                structure.atom_count(),
            ))
        },
    )
}

fn hydrogen_geometry() -> AnalysisDescriptor {
    AnalysisDescriptor::new("hydrogen-geometry", "1")
        .forbidding(ForbiddenResolution::new(
            PolicyField::Hydrogens,
            "the geometry is the hydrogen",
            |policy| matches!(policy.hydrogens, HydrogenPolicy::Exclude),
        ))
        .requiring(Requirement::ExplicitHydrogens)
        .estimating("the geometry of a modelled hydrogen")
}

#[test]
fn a_resolution_the_analysis_forbids_is_refused_not_answered_empty() {
    let policy = AnalysisPolicy {
        hydrogens: HydrogenPolicy::Exclude,
        ..AnalysisPolicy::default()
    };
    let result = analyse_structure(
        &structure(true),
        &policy,
        &atoms_kernel(hydrogen_geometry()),
        &ExecutionContext::default(),
    );
    let Err(GovernedAnalysisError::ForbiddenResolution {
        analysis,
        field,
        reason,
    }) = result
    else {
        panic!("excluding hydrogens must be refused for a hydrogen geometry");
    };
    assert_eq!(&*analysis, "hydrogen-geometry");
    assert_eq!(field, "hydrogens");
    assert!(reason.contains("hydrogen"));
}

#[test]
fn an_input_without_the_required_information_has_no_answer() {
    let Ok(result) = analyse_structure(
        &structure(false),
        &AnalysisPolicy::default(),
        &atoms_kernel(hydrogen_geometry()),
        &ExecutionContext::default(),
    ) else {
        panic!("a missing requirement is an outcome, not an error");
    };
    assert!(result.value().is_none());
    let reason = result
        .indeterminacy()
        .map(ToString::to_string)
        .unwrap_or_default();
    assert!(reason.contains("no hydrogen"), "{reason}");
}

#[test]
fn an_input_with_the_required_information_is_answered_and_names_its_estimand() {
    let Ok(result) = analyse_structure(
        &structure(true),
        &AnalysisPolicy::default(),
        &atoms_kernel(hydrogen_geometry()),
        &ExecutionContext::default(),
    ) else {
        panic!("a structure with hydrogens satisfies the requirement");
    };
    assert_eq!(result.value().copied(), Some(2));
    assert_eq!(
        result.provenance.estimand(),
        Some("the geometry of a modelled hydrogen")
    );
}

#[test]
fn an_analysis_that_states_nothing_forbids_nothing() {
    let policy = AnalysisPolicy {
        hydrogens: HydrogenPolicy::Exclude,
        ..AnalysisPolicy::default()
    };
    let Ok(result) = analyse_structure(
        &structure(true),
        &policy,
        &atoms_kernel(AnalysisDescriptor::new("atoms", "1")),
        &ExecutionContext::default(),
    ) else {
        panic!("an unconstrained analysis runs under any hydrogen rule");
    };
    assert_eq!(result.value().copied(), Some(1));
    assert_eq!(result.provenance.estimand(), None);
}
