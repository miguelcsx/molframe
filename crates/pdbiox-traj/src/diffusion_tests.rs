use super::*;

fn distances() -> EnsembleDistanceMatrix {
    EnsembleDistanceMatrix {
        size: 3,
        values: vec![0.0, 1.0, 4.0, 1.0, 0.0, 3.0, 4.0, 3.0, 0.0].into_boxed_slice(),
    }
}

#[test]
fn diffusion_excludes_the_constant_mode() {
    let Ok(map) = diffusion_map(&distances(), 2.0, 1, 2) else {
        panic!("valid distance matrix must embed");
    };
    assert_eq!(map.eigenvalues.len(), 2);
    assert_eq!(map.coordinates.len(), 3);
    assert!(map.eigenvalues[0] < 1.0);
}

#[test]
fn epsilon_must_be_positive() {
    assert_eq!(
        diffusion_map(&distances(), 0.0, 1, 2).err(),
        Some(EnsembleGeometryError::InvalidParameter)
    );
}
