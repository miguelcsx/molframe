use super::*;

const TRIANGLE: [[f32; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

#[test]
fn a_set_fits_itself_exactly_and_the_transform_changes_nothing() {
    let Ok(fit) = superpose(&TRIANGLE, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(fit.rmsd < 1e-6, "got {}", fit.rmsd);
    for point in TRIANGLE {
        let moved = fit.transform.apply(point);
        for axis in 0..3 {
            assert!((moved[axis] - point[axis]).abs() < 1e-5);
        }
    }
}

#[test]
fn a_translated_copy_fits_exactly() {
    let moved: Vec<[f32; 3]> = TRIANGLE
        .iter()
        .map(|p| [p[0] + 10.0, p[1] - 5.0, p[2] + 2.0])
        .collect();
    let Ok(fit) = superpose(&moved, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(fit.rmsd < 1e-5, "got {}", fit.rmsd);
}

#[test]
fn a_rotated_copy_fits_exactly_and_the_fit_is_a_rotation_not_a_reflection() {
    // A quarter turn about z.
    let turned: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [-p[1], p[0], p[2]]).collect();
    let Ok(fit) = superpose(&turned, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(fit.rmsd < 1e-5, "got {}", fit.rmsd);
    assert!(
        (fit.transform.determinant() - 1.0).abs() < 1e-9,
        "a fit must never mirror the structure"
    );
}

#[test]
fn a_mirrored_copy_is_not_fitted_by_reflecting_it() {
    let mirrored: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [-p[0], p[1], p[2]]).collect();
    let Ok(fit) = superpose(&mirrored, &TRIANGLE) else {
        panic!("expected a fit")
    };
    assert!(
        (fit.transform.determinant() - 1.0).abs() < 1e-9,
        "the transform must stay a rotation"
    );
    assert!(
        fit.rmsd > 0.1,
        "a mirror image genuinely does not fit, got {}",
        fit.rmsd
    );
}

#[test]
fn deviation_is_symmetric_and_zero_against_itself() {
    assert_eq!(rmsd(&TRIANGLE, &TRIANGLE).ok(), Some(0.0));
    let other: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [p[0] + 1.0, p[1], p[2]]).collect();
    let (Ok(forward), Ok(backward)) = (rmsd(&TRIANGLE, &other), rmsd(&other, &TRIANGLE)) else {
        panic!("expected both")
    };
    assert!((forward - backward).abs() < 1e-12);
    assert!((forward - 1.0).abs() < 1e-6);
}

#[test]
fn flat_coordinates_use_the_same_kernel_without_truncation() {
    let flat = TRIANGLE.into_iter().flatten().collect::<Vec<_>>();
    assert_eq!(rmsd_flat(&flat, &flat), Ok(0.0));
    assert_eq!(
        rmsd_flat(&flat[..flat.len() - 1], &flat[..flat.len() - 1]),
        Err(SuperposeError::LengthMismatch)
    );
}

#[test]
fn sets_of_different_sizes_are_refused_rather_than_truncated() {
    assert_eq!(
        rmsd(&TRIANGLE, &TRIANGLE[..2]),
        Err(SuperposeError::LengthMismatch)
    );
    assert_eq!(
        superpose(&TRIANGLE, &TRIANGLE[..2]).err(),
        Some(SuperposeError::LengthMismatch)
    );
}

#[test]
fn too_few_points_to_fix_a_rotation_are_refused() {
    let pair = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    assert_eq!(
        superpose(&pair, &pair).err(),
        Some(SuperposeError::TooFewPoints)
    );
}

#[test]
fn the_same_pair_of_sets_always_fits_the_same_way() {
    let moved: Vec<[f32; 3]> = TRIANGLE.iter().map(|p| [p[1], p[2], p[0]]).collect();
    let (Ok(first), Ok(second)) = (superpose(&moved, &TRIANGLE), superpose(&moved, &TRIANGLE))
    else {
        panic!("expected both")
    };
    assert!(
        first
            .transform
            .rotation
            .iter()
            .flatten()
            .zip(second.transform.rotation.iter().flatten())
            .all(|(left, right)| left.to_bits() == right.to_bits())
    );
    assert!(
        first
            .transform
            .translation
            .iter()
            .zip(second.transform.translation)
            .all(|(left, right)| left.to_bits() == right.to_bits())
    );
    assert_eq!(first.rmsd.to_bits(), second.rmsd.to_bits());
}

#[test]
fn caller_selected_fit_controls_are_validated_and_enforced() {
    assert_eq!(
        superpose_with_options(
            &TRIANGLE,
            &TRIANGLE,
            SuperposeOptions {
                collinear_relative_tolerance: 0.0,
                eigen: eigen::EigenOptions::standard(),
            },
        ),
        Err(SuperposeError::InvalidOptions)
    );
    assert_eq!(
        superpose_with_options(
            &TRIANGLE,
            &TRIANGLE,
            SuperposeOptions {
                collinear_relative_tolerance: 1e-12,
                eigen: eigen::EigenOptions {
                    relative_tolerance: 1e-14,
                    maximum_sweeps: 0,
                },
            },
        ),
        Err(SuperposeError::Eigen(eigen::EigenError::InvalidOptions))
    );
}

fn cloud(seed: &mut u64, points: usize) -> Vec<[f32; 3]> {
    let mut next = || {
        *seed ^= *seed << 13;
        *seed ^= *seed >> 7;
        *seed ^= *seed << 17;
        // Top 24 bits give an exact f32 in [0, 1).
        let bits = u32::try_from(*seed >> 40).unwrap_or(0);
        #[allow(clippy::cast_possible_truncation)]
        let value = (f64::from(bits) / 16_777_216.0 * 20.0 - 10.0) as f32;
        value
    };
    (0..points).map(|_| [next(), next(), next()]).collect()
}

#[test]
fn qcp_agrees_with_the_jacobi_path_on_random_clouds() {
    let mut seed = 0x9E37_79B9_7F4A_7C15_u64;
    for _ in 0..100 {
        let mobile = cloud(&mut seed, 24);
        let reference = cloud(&mut seed, 24);
        let options = SuperposeOptions::standard();
        let Ok((statistics, count)) = prepared_statistics(&mobile, &reference, options) else {
            panic!("expected statistics")
        };
        let key = key_matrix(&statistics.covariance);
        let Some((lambda, quaternion)) =
            qcp::dominant_eigenpair(&key, initial_correlation(&statistics))
        else {
            panic!("qcp should converge on a generic cloud")
        };
        let Ok(jacobi) = eigen::symmetric_with_options(key, options.eigen) else {
            panic!("jacobi should converge")
        };

        let qcp_rmsd = optimal_rmsd(&statistics, lambda, count);
        let jacobi_rmsd = optimal_rmsd(&statistics, jacobi.values[0], count);
        assert!(
            (qcp_rmsd - jacobi_rmsd).abs() < 1e-9,
            "{qcp_rmsd} vs {jacobi_rmsd}"
        );

        let (Some(a), Some(b)) = (rotation_from(quaternion), rotation_from(jacobi.dominant()))
        else {
            panic!("expected rotations")
        };
        for row in 0..3 {
            for column in 0..3 {
                assert!((a[row][column] - b[row][column]).abs() < 1e-9);
            }
        }
    }
}

#[test]
fn rmsd_after_fit_equals_the_fitted_rmsd() {
    let mut seed = 0x1234_5678_9ABC_DEF1_u64;
    for _ in 0..20 {
        let mobile = cloud(&mut seed, 30);
        let reference = cloud(&mut seed, 30);
        let (Ok(fit), Ok(value)) = (
            superpose(&mobile, &reference),
            rmsd_after_fit(&mobile, &reference),
        ) else {
            panic!("expected fits")
        };
        assert!((fit.rmsd - value).abs() < 1e-9);
    }
    assert_eq!(
        rmsd_after_fit(&TRIANGLE[..2], &TRIANGLE[..2]),
        Err(SuperposeError::TooFewPoints)
    );
}

#[test]
fn a_mirrored_regular_tetrahedron_falls_back_and_still_gives_a_proper_rotation() {
    // Its covariance with its mirror image makes the top eigenvalue triple.
    let tetrahedron: [[f32; 3]; 4] = [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    let mirrored: Vec<[f32; 3]> = tetrahedron.iter().map(|p| [p[0], p[1], -p[2]]).collect();

    let options = SuperposeOptions::standard();
    let Ok((statistics, _)) = prepared_statistics(&mirrored, &tetrahedron, options) else {
        panic!("expected statistics")
    };
    let key = key_matrix(&statistics.covariance);
    assert!(
        qcp::dominant_eigenpair(&key, initial_correlation(&statistics)).is_none(),
        "a repeated top eigenvalue must be left to the fallback"
    );

    let Ok(fit) = superpose(&mirrored, &tetrahedron) else {
        panic!("expected a fit")
    };
    assert!((fit.transform.determinant() - 1.0).abs() < 1e-9);
    assert!(fit.rmsd > 0.1);
    let Ok(value) = rmsd_after_fit(&mirrored, &tetrahedron) else {
        panic!("expected an rmsd")
    };
    // Newton converges only linearly on a repeated root, so the λ-only path
    // is less exact here than the fallback.
    assert!((fit.rmsd - value).abs() < 1e-5);
}
