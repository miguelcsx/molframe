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

#[test]
fn the_gram_matrix_matches_the_matrix_product_it_replaces() {
    // Deliberately wider than it is tall, the shape the Gram path exists for.
    let observations = 5;
    let features = 23;
    let matrix = DMatrix::from_fn(observations, features, |row, column| {
        let row = f64::from(u32::try_from(row).expect("small row"));
        let column = f64::from(u32::try_from(column).expect("small column"));
        (row + 1.0).mul_add(0.7, column * -0.31) + (row * column).sin()
    });
    let divisor = 4.0;

    let expected = &matrix * matrix.transpose() / divisor;
    let actual = gram_matrix(&matrix, divisor);

    assert_eq!(actual.nrows(), observations);
    assert_eq!(actual.ncols(), observations);
    for row in 0..observations {
        for column in 0..observations {
            let difference = actual[(row, column)] - expected[(row, column)];
            assert!(
                difference.abs() < 1.0e-9,
                "gram[{row}][{column}] was {} not {}",
                actual[(row, column)],
                expected[(row, column)]
            );
        }
    }
}

#[test]
fn the_gram_matrix_is_symmetric_and_survives_a_zero_column() {
    let matrix = DMatrix::from_fn(3, 4, |row, column| {
        if column == 2 {
            0.0
        } else {
            f64::from(u32::try_from(row + column).expect("small index"))
        }
    });
    let gram = gram_matrix(&matrix, 2.0);

    for row in 0..3 {
        for column in 0..3 {
            let difference = gram[(row, column)] - gram[(column, row)];
            assert!(difference.abs() < 1.0e-12, "gram must be symmetric");
        }
    }
}
