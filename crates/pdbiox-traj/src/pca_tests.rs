use super::*;

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
