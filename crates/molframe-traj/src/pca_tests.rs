use super::*;
use std::mem::size_of;

#[test]
fn cartesian_pca_finds_the_only_varying_axis() {
    let frames = vec![
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
    ];
    let Ok(result) = cartesian_pca(&frames, CartesianFit::None, 1, 4096) else {
        panic!("valid observations must decompose");
    };
    assert_eq!(result.components.len(), 1);
    assert!(result.components[0][3] > 0.99);
}

#[test]
fn borrowed_cartesian_pca_matches_owned_frames() {
    let positions = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let view = crate::FrameView::new(&positions, 3, 3).expect("valid borrowed frames");
    let owned = vec![
        positions[..3].to_vec(),
        positions[3..6].to_vec(),
        positions[6..].to_vec(),
    ];
    assert_eq!(
        cartesian_pca_view(view, CartesianFit::None, 1, 4096),
        cartesian_pca(&owned, CartesianFit::None, 1, 4096)
    );
}

#[test]
fn dihedral_pca_is_wrapping_invariant() {
    let Ok(left) = PeriodicAngle::from_radians(-3.0) else {
        panic!("finite angle");
    };
    let Ok(right) = PeriodicAngle::from_radians(-3.0 + std::f64::consts::TAU) else {
        panic!("finite angle");
    };
    let observations = vec![vec![left], vec![right], vec![left]];
    let Ok(result) = dihedral_pca(&observations, 1, 4096) else {
        panic!("valid observations must decompose");
    };
    assert!(result.eigenvalues[0] < 1.0e-20);
}

#[test]
fn pca_memory_is_bounded() {
    let frames = vec![vec![[0.0, 0.0, 0.0]; 3], vec![[1.0, 0.0, 0.0]; 3]];
    assert_eq!(
        cartesian_pca(&frames, CartesianFit::None, 1, 1).err(),
        Some(EnsembleGeometryError::MemoryLimit)
    );
}

#[test]
fn pca_rejects_the_full_workspace_before_attempting_a_fit() {
    let frames = vec![vec![[0.0, 0.0, 0.0]; 128], vec![[0.0, 0.0, 0.0]; 128]];
    assert_eq!(
        cartesian_pca(&frames, CartesianFit::Reference(&frames[0]), 1, 1,),
        Err(EnsembleGeometryError::MemoryLimit)
    );
}

#[test]
fn exact_workspace_ceiling_is_accepted() {
    let frames = vec![
        vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vec![[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    ];
    let observations = frames.len();
    let features = 9;
    let components = 1;
    let data = observations * features;
    let eigensolver = 2 * features.min(observations).pow(2);
    let result = features * components + observations * components + features + components;
    let exact = (data + eigensolver + result) * size_of::<f64>();
    assert_eq!(
        cartesian_pca(&frames, CartesianFit::None, components, exact - 1),
        Err(EnsembleGeometryError::MemoryLimit)
    );
    assert!(cartesian_pca(&frames, CartesianFit::None, components, exact).is_ok());
}
