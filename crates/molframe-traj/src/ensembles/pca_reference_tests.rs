//! Cartesian PCA against a dense `nalgebra` covariance reference.

use super::*;
use nalgebra::{DMatrix, SymmetricEigen};

/// Deterministic low-rank-plus-noise frames with a well-separated spectrum.
fn real(value: usize) -> f64 {
    f64::from(u32::try_from(value).expect("small"))
}

fn synthetic(frames: usize, atoms: usize) -> Vec<Vec<[f32; 3]>> {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f32::from(u16::try_from(state >> 48).expect("sixteen bits")) / 65_536.0 - 0.5
    };
    let base: Vec<[f32; 3]> = (0..atoms).map(|_| [next(), next(), next()]).collect();
    let modes: Vec<Vec<[f32; 3]>> = (0..3)
        .map(|_| (0..atoms).map(|_| [next(), next(), next()]).collect())
        .collect();
    (0..frames)
        .map(|_| {
            let weights = [next() * 6.0, next() * 3.0, next()];
            (0..atoms)
                .map(|atom| {
                    let mut point = base[atom];
                    for (weight, mode) in weights.iter().zip(&modes) {
                        for axis in 0..3 {
                            point[axis] += weight * mode[atom][axis];
                        }
                    }
                    point.map(|value| value + 0.01 * next())
                })
                .collect()
        })
        .collect()
}

struct Reference {
    eigenvalues: Vec<f64>,
    axes: Vec<Vec<f64>>,
    projections: Vec<Vec<f64>>,
}

/// Dense covariance eigen-decomposition, materialising everything.
fn dense_reference(frames: &[Vec<[f32; 3]>], count: usize) -> Reference {
    let rows = frames.len();
    let columns = frames[0].len() * 3;
    let mut matrix = DMatrix::from_fn(rows, columns, |row, column| {
        f64::from(frames[row][column / 3][column % 3])
    });
    for column in 0..columns {
        let mean = matrix.column(column).sum() / real(rows);
        matrix.column_mut(column).add_scalar_mut(-mean);
    }
    let covariance = matrix.tr_mul(&matrix) / (real(rows) - 1.0);
    let decomposition = SymmetricEigen::new(covariance);
    let mut order: Vec<usize> = (0..columns).collect();
    order.sort_by(|&a, &b| decomposition.eigenvalues[b].total_cmp(&decomposition.eigenvalues[a]));
    let mut eigenvalues = Vec::new();
    let mut axes = Vec::new();
    let mut projections = vec![Vec::new(); rows];
    for &index in order.iter().take(count) {
        let mut direction: Vec<f64> = decomposition
            .eigenvectors
            .column(index)
            .iter()
            .copied()
            .collect();
        let pivot = direction.iter().copied().fold(0.0_f64, |best, value| {
            if value.abs() > best.abs() {
                value
            } else {
                best
            }
        });
        if pivot < 0.0 {
            for value in &mut direction {
                *value = -*value;
            }
        }
        let projected = &matrix * DMatrix::from_column_slice(columns, 1, &direction);
        for (row, target) in projections.iter_mut().enumerate() {
            target.push(projected[(row, 0)]);
        }
        eigenvalues.push(decomposition.eigenvalues[index]);
        axes.push(direction);
    }
    Reference {
        eigenvalues,
        axes,
        projections,
    }
}

fn assert_matches_dense(frames: &[Vec<[f32; 3]>], count: usize) {
    let reference = dense_reference(frames, count);
    let result = cartesian_pca(frames, CartesianFit::None, count, usize::MAX / 2)
        .expect("synthetic frames decompose");
    assert_eq!(result.eigenvalues.len(), count);
    for component in 0..count {
        let expected = reference.eigenvalues[component];
        assert!(
            (result.eigenvalues[component] - expected).abs() <= 1.0e-9 * expected.max(1.0),
            "eigenvalue {component}: {} vs {expected}",
            result.eigenvalues[component]
        );
        for (got, want) in result.components[component]
            .iter()
            .zip(&reference.axes[component])
        {
            assert!(
                (got - want).abs() < 1.0e-7,
                "axis {component}: {got} vs {want}"
            );
        }
        for (row, projection) in result.projections.iter().enumerate() {
            let want = reference.projections[row][component];
            assert!((projection[component] - want).abs() < 1.0e-6);
        }
    }
}

#[test]
fn covariance_route_matches_the_dense_reference() {
    // frames > 3N: the covariance route.
    assert_matches_dense(&synthetic(120, 12), 3);
}

#[test]
fn snapshot_route_matches_the_dense_reference() {
    // frames < 3N: the snapshot route; the reference still builds p x p.
    assert_matches_dense(&synthetic(40, 60), 3);
}

#[test]
fn results_are_independent_of_the_worker_budget() {
    let frames = synthetic(90, 40);
    let first = cartesian_pca(&frames, CartesianFit::None, 4, usize::MAX / 2).expect("fit");
    let again = cartesian_pca(&frames, CartesianFit::None, 4, usize::MAX / 2).expect("fit");
    assert_eq!(first, again);
}

#[test]
fn fitted_pca_matches_a_pre_fitted_unfitted_run() {
    let frames = synthetic(30, 25);
    let reference = frames[0].clone();
    let fitted: Vec<Vec<[f32; 3]>> = frames
        .iter()
        .map(|frame| {
            let fit = molframe_geom::superpose(frame, &reference).expect("fit");
            frame
                .iter()
                .map(|&point| fit.transform.apply(point))
                .collect()
        })
        .collect();
    let direct = cartesian_pca(
        &frames,
        CartesianFit::Reference(&reference),
        3,
        usize::MAX / 2,
    )
    .expect("fit");
    let manual = cartesian_pca(&fitted, CartesianFit::None, 3, usize::MAX / 2).expect("fit");
    assert_eq!(direct, manual);
}
