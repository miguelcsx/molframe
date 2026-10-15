use super::*;
use crate::shared_tests::worker_context as context;
use nalgebra::DMatrix;

fn real(value: usize) -> f64 {
    f64::from(u32::try_from(value).expect("small"))
}

fn data(rows: usize, columns: usize) -> Vec<f64> {
    (0..rows * columns)
        .map(|index| {
            let (row, column) = (index / columns, index % columns);
            let seed = (row * 31 + column * 17) % 23;
            f64::from(u8::try_from(seed).expect("small")) * 0.37 - 4.0
        })
        .collect()
}

fn means(values: &[f64], rows: usize, columns: usize) -> Vec<f64> {
    (0..columns)
        .map(|column| {
            (0..rows)
                .map(|row| values[row * columns + column])
                .sum::<f64>()
                / real(rows)
        })
        .collect()
}

fn centred_dense(values: &[f64], rows: usize, columns: usize, mean: &[f64]) -> DMatrix<f64> {
    DMatrix::from_fn(rows, columns, |row, column| {
        values[row * columns + column] - mean[column]
    })
}

fn to_dense(matrix: &Mat<f64>) -> DMatrix<f64> {
    DMatrix::from_fn(matrix.nrows(), matrix.ncols(), |row, column| {
        matrix[(row, column)]
    })
}

#[test]
fn covariance_and_gram_match_centred_nalgebra_products() {
    let ctx = context(4);
    let (rows, columns) = (700, 40);
    let values = data(rows, columns);
    let mean = means(&values, rows, columns);
    let observations = Observations::new(&values, rows, columns, &mean).expect("shape");
    let dense = centred_dense(&values, rows, columns, &mean);
    let expected = dense.tr_mul(&dense) / 699.0;
    let got = to_dense(&covariance(observations, 699.0, &ctx).expect("covariance"));
    assert!((&got - &expected).abs().max() < 1e-9);

    let (rows, columns) = (30, 600);
    let values = data(rows, columns);
    let mean = means(&values, rows, columns);
    let observations = Observations::new(&values, rows, columns, &mean).expect("shape");
    let dense = centred_dense(&values, rows, columns, &mean);
    let expected = &dense * dense.transpose() / 29.0;
    let got = to_dense(&gram(observations, 29.0, &ctx).expect("gram"));
    assert!((&got - &expected).abs().max() < 1e-9);
}

#[test]
fn products_are_bit_identical_at_every_worker_count() {
    let (rows, columns) = (1500, 33);
    let values = data(rows, columns);
    let mean = means(&values, rows, columns);
    let observations = Observations::new(&values, rows, columns, &mean).expect("shape");
    let weights = Mat::from_fn(columns, 2, |row, column| real(row * (column + 1)) * 0.5);
    let reference = covariance(observations, 7.0, &context(1)).expect("serial");
    let projected = project(observations, weights.as_ref(), &context(1)).expect("serial");
    for workers in [2, 5, 16] {
        let ctx = context(workers);
        assert_eq!(
            covariance(observations, 7.0, &ctx).expect("parallel"),
            reference
        );
        assert_eq!(
            project(observations, weights.as_ref(), &ctx).expect("parallel"),
            projected
        );
    }
}

#[test]
fn projection_and_back_projection_match_dense_products() {
    let (rows, columns) = (300, 700);
    let values = data(rows, columns);
    let mean = means(&values, rows, columns);
    let observations = Observations::new(&values, rows, columns, &mean).expect("shape");
    let dense = centred_dense(&values, rows, columns, &mean);
    let weights = Mat::from_fn(columns, 3, |row, column| real((row + column) % 7) - 3.0);
    let got = project(observations, weights.as_ref(), &context(3)).expect("project");
    let expected = &dense * to_dense(&weights);
    for row in [0, 150, 299] {
        for column in 0..3 {
            assert!((got[row][column] - expected[(row, column)]).abs() < 1e-8);
        }
    }
    let snapshots = Mat::from_fn(rows, 2, |row, column| real((row * 3 + column) % 5) - 2.0);
    let back = back_project(observations, snapshots.as_ref(), &context(3)).expect("back");
    let expected = dense.tr_mul(&to_dense(&snapshots));
    assert!((to_dense(&back) - expected).abs().max() < 1e-8);
}

#[test]
fn shape_mismatch_is_rejected() {
    assert!(Observations::new(&[1.0; 5], 2, 3, &[0.0; 3]).is_err());
    assert!(Observations::new(&[1.0; 6], 2, 3, &[0.0; 2]).is_err());
}
