use super::*;

#[test]
fn a_square_matrix_is_symmetric_with_an_exact_zero_diagonal() {
    let points = [[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 12.0]];
    let Ok(matrix) = distance_matrix(&points) else {
        panic!("small square matrix must fit");
    };
    assert_eq!((matrix.rows(), matrix.columns()), (3, 3));
    assert_eq!(matrix.get(0, 0), Some(0.0));
    assert_eq!(matrix.get(0, 1), Some(5.0));
    assert_eq!(matrix.get(1, 0), matrix.get(0, 1));
    assert_eq!(matrix.get(0, 2), Some(12.0));
    assert_eq!(matrix.get(3, 0), None);
}

#[test]
fn a_rectangular_matrix_keeps_row_major_shape_and_order() {
    let left = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let right = [[0.0, 1.0, 0.0]];
    let Ok(matrix) = distance_matrix_between(&left, &right) else {
        panic!("small rectangular matrix must fit");
    };
    assert_eq!((matrix.rows(), matrix.columns()), (2, 1));
    let expected = [1.0, 2.0f64.sqrt()];
    assert!(
        matrix
            .as_slice()
            .iter()
            .zip(expected)
            .all(|(actual, expected)| (actual - expected).abs() < f64::EPSILON)
    );
}

#[test]
fn values_can_be_moved_without_changing_row_major_order() {
    let left = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
    let right = [[0.0, 1.0, 0.0]];
    let Ok(matrix) = distance_matrix_between(&left, &right) else {
        panic!("small rectangular matrix must fit");
    };
    assert_eq!(matrix.into_values(), vec![1.0, 2.0f64.sqrt()]);
}

#[test]
fn caller_owned_output_matches_the_allocating_api() {
    let points = [[0.0, 0.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 12.0]];
    let mut output = vec![f64::NAN; points.len() * points.len()];
    if let Err(error) = distance_matrix_into(&points, &mut output) {
        panic!("into calculation failed: {error}");
    }
    let matrix = match distance_matrix(&points) {
        Ok(value) => value,
        Err(error) => panic!("matrix calculation failed: {error}"),
    };
    assert_eq!(output, matrix.as_slice());
    assert_eq!(
        distance_matrix_into(&points, &mut output[..8]),
        Err(MatrixError::OutputLength)
    );
}

#[test]
fn parallel_row_blocks_are_bitwise_worker_count_independent() {
    let points = (0_u16..300)
        .map(|index| {
            let value = f32::from(index);
            [value, value.mul_add(0.25, 1.0), value.mul_add(-0.5, 2.0)]
        })
        .collect::<Vec<_>>();
    let serial_context = match ExecutionContext::builder().worker_budget(1).build() {
        Ok(value) => value,
        Err(error) => panic!("serial context failed: {error}"),
    };
    let parallel_context = ExecutionContext::default();
    let serial = match distance_matrix_with_context(&points, &serial_context) {
        Ok(value) => value,
        Err(error) => panic!("serial matrix failed: {error}"),
    };
    let parallel = match distance_matrix_with_context(&points, &parallel_context) {
        Ok(value) => value,
        Err(error) => panic!("parallel matrix failed: {error}"),
    };
    assert_eq!(serial.as_slice(), parallel.as_slice());
    for workers in [2, 4, 8] {
        let context = ExecutionContext::builder()
            .worker_budget(workers)
            .build()
            .expect("a positive worker budget is valid");
        let matrix = distance_matrix_with_context(&points, &context).expect("fits the budget");
        assert_eq!(matrix.as_slice(), serial.as_slice(), "{workers} workers");
    }
}

#[test]
fn a_budget_for_the_matrix_alone_suffices_and_a_smaller_one_is_refused() {
    let points = (0_u16..300)
        .map(|index| [f32::from(index), 0.0, 0.0])
        .collect::<Vec<_>>();
    let matrix_bytes = 300 * 300 * 8;
    let budgeted = |bytes| {
        ExecutionContext::builder()
            .worker_budget(4)
            .memory_budget(molframe_core::MemoryBudget::new(bytes).expect("a positive budget"))
            .scratch_policy(molframe_core::ScratchPolicy::new(0))
            .build()
            .expect("a valid context")
    };
    // One matrix plus a few row blocks in flight; the old path needed two matrices.
    let roomy = budgeted(matrix_bytes + matrix_bytes / 2);
    assert!(distance_matrix_with_context(&points, &roomy).is_ok());
    assert_eq!(roomy.reserved_bytes(), 0);
    let tight = budgeted(matrix_bytes - 1);
    assert_eq!(
        distance_matrix_with_context(&points, &tight).map(|_| ()),
        Err(MatrixError::MemoryBudget)
    );
}
