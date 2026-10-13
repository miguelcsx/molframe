use super::*;

fn matrix(rows: usize, columns: usize) -> DMatrix<f64> {
    DMatrix::from_fn(rows, columns, |row, column| {
        let seed = (row * 31 + column * 17) % 23;
        f64::from(u8::try_from(seed).expect("small")) * 0.37 - 4.0
    })
}

fn context(workers: usize) -> ExecutionContext {
    ExecutionContext::builder()
        .worker_budget(workers)
        .build()
        .expect("context")
}

#[test]
fn covariance_and_gram_match_nalgebra_products() {
    let m = matrix(700, 40);
    let ctx = context(4);
    let expected = m.tr_mul(&m) / 699.0;
    let got = covariance(&m, 699.0, &ctx).expect("covariance");
    assert!((&got - &expected).abs().max() < 1e-9);
    let wide = matrix(30, 600);
    let expected = &wide * wide.transpose() / 29.0;
    let got = gram(&wide, 29.0, &ctx).expect("gram");
    assert!((&got - &expected).abs().max() < 1e-9);
}

#[test]
fn products_are_bit_identical_at_every_worker_count() {
    let m = matrix(1500, 33);
    let reference = covariance(&m, 7.0, &context(1)).expect("serial");
    let components: Vec<Box<[f64]>> = vec![(0..33).map(f64::from).collect(), vec![0.5; 33].into()];
    let projected = project(&m, &components, &context(1)).expect("serial");
    for workers in [2, 5, 16] {
        let ctx = context(workers);
        assert_eq!(covariance(&m, 7.0, &ctx).expect("parallel"), reference);
        assert_eq!(project(&m, &components, &ctx).expect("parallel"), projected);
    }
}

#[test]
fn projection_matches_the_row_dot_products() {
    let m = matrix(300, 12);
    let components: Vec<Box<[f64]>> = vec![(0..12).map(f64::from).collect()];
    let got = project(&m, &components, &context(3)).expect("project");
    for row in [0, 150, 299] {
        let expected: f64 = (0..12)
            .map(|c| m[(row, c)] * f64::from(u8::try_from(c).expect("small")))
            .sum();
        assert!((got[row][0] - expected).abs() < 1e-9);
    }
}
