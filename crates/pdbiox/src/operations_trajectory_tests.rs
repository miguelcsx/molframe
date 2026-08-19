use crate::{AnalysisPolicy, FrameInput, Plan, PlanInput, PlanOperation, TrajectoryRequest};

static TRAJECTORY_POSITIONS: [[f32; 3]; 9] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 1.0, 0.0],
    [2.0, 1.0, 0.0],
    [1.0, 2.0, 0.0],
    [2.0, 2.0, 0.0],
    [3.0, 2.0, 0.0],
    [2.0, 3.0, 0.0],
];

#[test]
fn trajectory_plan_rmsd_and_displacement_match_borrowed_kernels() {
    let (frames, input, policy, view) = trajectory_fixture();
    let direct_rmsd = crate::traj::analyse_rmsd_to_reference_view(
        view,
        0,
        crate::traj::FrameAlignment::Rigid,
        &policy,
    )
    .expect("direct RMSD");
    let direct_displacement =
        crate::traj::analyse_mean_squared_displacement_view(view, &[], 2, &policy)
            .expect("direct MSD");
    let mut plan = Plan::new();
    plan.add_trajectory(
        "rmsd",
        TrajectoryRequest::RmsdToReference {
            frames,
            reference: 0,
            alignment: crate::traj::FrameAlignment::Rigid,
            policy: policy.clone(),
        },
    )
    .expect("RMSD operation");
    plan.add(
        "displacement",
        PlanOperation::Trajectory(Box::new(TrajectoryRequest::MeanSquaredDisplacement {
            frames,
            atoms: None,
            maximum_lag: 2,
            policy,
        })),
    )
    .expect("displacement operation");
    let values = execute_trajectory_plan(&plan, &input);
    assert_rmsd(&values, &direct_rmsd);
    assert_displacement(&values, &direct_displacement);
}

#[test]
fn trajectory_plan_pairwise_and_mean_match_borrowed_kernels() {
    let (frames, input, policy, view) = trajectory_fixture();
    let direct_pairwise = crate::traj::analyse_pairwise_fitted_rmsd_view(view, 1_024, &policy)
        .expect("direct pairwise RMSD");
    let direct_mean =
        crate::traj::analyse_generalized_procrustes_mean_view(view, 1.0e-6, 20, &policy)
            .expect("direct Procrustes mean");
    let mut plan = Plan::new();
    plan.add_trajectory(
        "pairwise",
        TrajectoryRequest::PairwiseFittedRmsd {
            frames,
            memory_limit: 1_024,
            policy: policy.clone(),
        },
    )
    .expect("pairwise operation");
    plan.add_trajectory(
        "mean",
        TrajectoryRequest::GeneralizedProcrustesMean {
            frames,
            tolerance: 1.0e-6,
            maximum_iterations: 20,
            policy,
        },
    )
    .expect("Procrustes operation");
    let values = execute_trajectory_plan(&plan, &input);
    assert_pairwise(&values, &direct_pairwise);
    assert_mean(&values, &direct_mean);
}

fn trajectory_fixture() -> (
    usize,
    [FrameInput<'static>; 1],
    AnalysisPolicy,
    crate::traj::FrameView<'static>,
) {
    let input = [FrameInput {
        positions: &TRAJECTORY_POSITIONS,
        frame_count: 3,
        atom_count: 3,
    }];
    let policy = AnalysisPolicy::default();
    let view =
        crate::traj::FrameView::new(&TRAJECTORY_POSITIONS, 3, 3).expect("borrowed frame fixture");
    (0, input, policy, view)
}

fn execute_trajectory_plan(
    plan: &Plan,
    input: &[FrameInput<'_>],
) -> std::collections::BTreeMap<Box<str>, crate::PlanValue> {
    let result = plan
        .execute(PlanInput {
            frames: input,
            ..Default::default()
        })
        .expect("native trajectory plan");
    result
        .entries
        .into_iter()
        .map(|entry| (entry.id, entry.value))
        .collect()
}

fn assert_rmsd(
    values: &std::collections::BTreeMap<Box<str>, crate::PlanValue>,
    expected: &crate::Analysis<Vec<f64>>,
) {
    let Some(crate::PlanValue::Trajectory(value)) = values.get("rmsd") else {
        panic!("missing RMSD result");
    };
    let crate::TrajectoryValue::RmsdToReference(actual) = value.as_ref() else {
        panic!("unexpected RMSD value");
    };
    assert_eq!(actual.value, expected.value);
    assert_eq!(
        actual.provenance.fingerprint(),
        expected.provenance.fingerprint()
    );
}

fn assert_displacement(
    values: &std::collections::BTreeMap<Box<str>, crate::PlanValue>,
    expected: &crate::Analysis<Vec<crate::traj::MeanSquaredDisplacement>>,
) {
    let Some(crate::PlanValue::Trajectory(value)) = values.get("displacement") else {
        panic!("missing displacement result");
    };
    let crate::TrajectoryValue::MeanSquaredDisplacement(actual) = value.as_ref() else {
        panic!("unexpected displacement value");
    };
    assert_eq!(actual.value, expected.value);
    assert_eq!(
        actual.provenance.fingerprint(),
        expected.provenance.fingerprint()
    );
}

fn assert_pairwise(
    values: &std::collections::BTreeMap<Box<str>, crate::PlanValue>,
    expected: &crate::Analysis<crate::traj::EnsembleDistanceMatrix>,
) {
    let Some(crate::PlanValue::Trajectory(value)) = values.get("pairwise") else {
        panic!("missing pairwise result");
    };
    let crate::TrajectoryValue::PairwiseFittedRmsd(actual) = value.as_ref() else {
        panic!("unexpected pairwise value");
    };
    assert_eq!(actual.value, expected.value);
    assert_eq!(
        actual.provenance.fingerprint(),
        expected.provenance.fingerprint()
    );
}

fn assert_mean(
    values: &std::collections::BTreeMap<Box<str>, crate::PlanValue>,
    expected: &crate::Analysis<Vec<[f32; 3]>>,
) {
    let Some(crate::PlanValue::Trajectory(value)) = values.get("mean") else {
        panic!("missing Procrustes result");
    };
    let crate::TrajectoryValue::GeneralizedProcrustesMean(actual) = value.as_ref() else {
        panic!("unexpected Procrustes value");
    };
    assert_eq!(actual.value, expected.value);
    assert_eq!(
        actual.provenance.fingerprint(),
        expected.provenance.fingerprint()
    );
}
