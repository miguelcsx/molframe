use super::*;

fn frames() -> Vec<Vec<[f32; 3]>> {
    vec![
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[4.0, 2.0, 1.0], [5.0, 2.0, 1.0], [4.0, 3.0, 1.0]],
    ]
}

#[test]
fn fitted_distances_remove_rigid_translation() {
    let Ok(matrix) = pairwise_fitted_rmsd(&frames(), 1024) else {
        panic!("valid frames must fit");
    };
    assert!(matrix.get(0, 1).is_some_and(|value| value < 1.0e-6));
}

#[test]
fn pairwise_allocation_is_bounded_before_work() {
    assert_eq!(
        pairwise_fitted_rmsd(&frames(), 1).err(),
        Some(EnsembleGeometryError::MemoryLimit)
    );
}

#[test]
fn procrustes_mean_is_translation_invariant() {
    let Ok(mean) = generalized_procrustes_mean(&frames(), 1.0e-6, 20) else {
        panic!("valid ensemble must converge");
    };
    assert!(rmsd(&mean, &frames()[0]).is_ok_and(|value| value < 1.0e-5));
}

#[test]
fn reference_rmsd_uses_explicit_alignment_policy() {
    let trajectory = frames()
        .into_iter()
        .enumerate()
        .map(|(frame, positions)| crate::Timestep {
            frame,
            positions,
            ..crate::Timestep::default()
        })
        .collect::<Vec<_>>();
    let Ok(raw) = rmsd_to_reference(&trajectory, 0, FrameAlignment::None) else {
        panic!("raw RMSD should be defined");
    };
    let Ok(fitted) = rmsd_to_reference(&trajectory, 0, FrameAlignment::Rigid) else {
        panic!("fitted RMSD should be defined");
    };
    assert!(raw[1] > 1.0);
    assert!(fitted[1] < 1.0e-6);
}

#[test]
fn reference_rmsd_rejects_absent_reference() {
    assert_eq!(
        rmsd_to_reference(&[], 0, FrameAlignment::None),
        Err(EnsembleGeometryError::Empty)
    );
}

#[test]
fn reference_rmsd_rejects_out_of_range_reference() {
    let trajectory = vec![crate::Timestep {
        positions: frames()[0].clone(),
        ..crate::Timestep::default()
    }];
    assert_eq!(
        rmsd_to_reference(&trajectory, 1, FrameAlignment::None),
        Err(EnsembleGeometryError::InvalidParameter)
    );
}
