use super::{analyse_structure, analyse_trajectory};
use crate::numeric::usize_to_f32;
use crate::policy_execution::{AnalysisDescriptor, FrameKernelResult, structure_kernel};
use pdbiox_core::contract::{
    AnalysisPolicy, AssemblyChoice, Coverage, MissingPolicy, ParameterValue, Status,
};
use pdbiox_core::io::{InputBuffer, ReadOptions};
use pdbiox_core::structure::Structure;
use pdbiox_traj::{Frame, Trajectory};
use std::convert::Infallible;

const SOURCE: &str = "data_s\n\
loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n\
_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n\
_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n\
ATOM 1 C C1 LIG A 1 0 0 0\n\
ATOM 2 C C2 LIG A 1 1 0 0\n";

fn structure() -> Structure {
    let input = InputBuffer::from_bytes(SOURCE.as_bytes().to_vec());
    match pdbiox_cif::read(&input, &ReadOptions::new()) {
        Ok((structure, _)) => structure,
        Err(findings) => panic!("fixture failed: {findings:?}"),
    }
}

fn x_sum_kernel() -> impl crate::policy_execution::StructureKernel<Output = f64, Error = Infallible>
{
    structure_kernel(
        AnalysisDescriptor::new("coordinate-x-sum", "1")
            .with_parameter("axis", ParameterValue::Text("x".into())),
        |structure: &Structure, _policy: &AnalysisPolicy| {
            let sum = structure
                .positions()
                .iter()
                .map(|position| f64::from(position[0]))
                .sum();
            Ok(FrameKernelResult::complete(sum, structure.atom_count()))
        },
    )
}

#[test]
fn unmaterialized_assembly_policy_is_rejected() {
    let policy = AnalysisPolicy {
        assembly: AssemblyChoice::Biological("1".into()),
        ..AnalysisPolicy::default()
    };
    let result = analyse_structure(&structure(), &policy, &x_sum_kernel());
    assert!(matches!(
        result,
        Err(super::GovernedAnalysisError::UnsupportedPolicyValue(
            "assembly"
        ))
    ));
}

#[test]
fn a_structure_is_the_one_frame_case_of_the_same_adapter() {
    let structure = structure();
    let policy = AnalysisPolicy::default();
    let kernel = x_sum_kernel();
    let Ok(single) = analyse_structure(&structure, &policy, &kernel) else {
        panic!("single-frame execution should succeed");
    };
    let trajectory = Trajectory::from_frames(vec![Frame {
        positions: structure.positions().to_vec(),
    }]);
    let Ok(series) = analyse_trajectory(&structure, &trajectory, &policy, &kernel, 1) else {
        panic!("trajectory execution should succeed");
    };
    assert_eq!(single.value.to_bits(), series.value[0].to_bits());
    assert_eq!(single.coverage, series.coverage);
    assert_eq!(
        single.provenance.fingerprint(),
        series.provenance.fingerprint()
    );
}

#[test]
fn serial_and_parallel_series_are_bit_identical_with_the_same_provenance() {
    let structure = structure();
    let trajectory = Trajectory::from_frames(
        (0..1_025)
            .map(|index| Frame {
                positions: vec![[usize_to_f32(index), 0.0, 0.0], [0.25, 0.0, 0.0]],
            })
            .collect(),
    );
    let policy = AnalysisPolicy::default();
    let kernel = x_sum_kernel();
    let outputs: Vec<_> = [1, 2, 4, 16]
        .into_iter()
        .map(|workers| {
            let Ok(result) = analyse_trajectory(&structure, &trajectory, &policy, &kernel, workers)
            else {
                panic!("valid deterministic execution");
            };
            (
                result
                    .value
                    .iter()
                    .map(|value| value.to_bits())
                    .collect::<Vec<_>>(),
                result.provenance.fingerprint(),
                result.coverage,
            )
        })
        .collect();
    assert!(outputs.windows(2).all(|pair| pair[0] == pair[1]));
}

#[test]
fn provenance_contains_the_kernel_identity_parameters_and_frame_count() {
    let structure = structure();
    let policy = AnalysisPolicy::default();
    let kernel = x_sum_kernel();
    let Ok(result) = analyse_structure(&structure, &policy, &kernel) else {
        panic!("analysis should succeed");
    };
    let Some(algorithm) = &result.provenance.algorithm else {
        panic!("algorithm identity must be recorded");
    };
    assert_eq!(algorithm.name(), "coordinate-x-sum");
    assert_eq!(algorithm.version(), "1");
    assert_eq!(result.provenance.input_source.to_string(), "(memory)");
    assert_eq!(
        result.provenance.parameters.get("axis"),
        Some(&ParameterValue::Text("x".into()))
    );
    assert_eq!(
        result.provenance.parameters.get("frame_count"),
        Some(&ParameterValue::Integer(1))
    );
}

#[test]
fn missing_data_policy_is_enforced_without_discarding_coverage() {
    let structure = structure();
    let kernel = structure_kernel(
        AnalysisDescriptor::new("incomplete", "1"),
        |_structure: &Structure, _policy: &AnalysisPolicy| {
            Ok::<_, Infallible>(FrameKernelResult::governed(
                3_u8,
                Status::Complete,
                Coverage {
                    intended: 2,
                    used: 1,
                    missing: 1,
                    ambiguous: 0,
                },
            ))
        },
    );
    let report = AnalysisPolicy::default().with_missing_atoms(MissingPolicy::Report);
    let Ok(result) = analyse_structure(&structure, &report, &kernel) else {
        panic!("report policy computes a partial value");
    };
    assert_eq!(result.status, Status::Partial);
    assert_eq!(result.coverage.missing, 1);

    let fail = AnalysisPolicy::default().with_missing_atoms(MissingPolicy::Fail);
    assert!(analyse_structure(&structure, &fail, &kernel).is_err());
}

#[test]
fn inconsistent_kernel_coverage_is_rejected_before_publication() {
    let structure = structure();
    let kernel = structure_kernel(
        AnalysisDescriptor::new("invalid-coverage", "1"),
        |_structure: &Structure, _policy: &AnalysisPolicy| {
            Ok::<_, Infallible>(FrameKernelResult::governed(
                0_u8,
                Status::Complete,
                Coverage {
                    intended: 1,
                    used: 1,
                    missing: 1,
                    ambiguous: 0,
                },
            ))
        },
    );
    assert!(analyse_structure(&structure, &AnalysisPolicy::default(), &kernel).is_err());
}
